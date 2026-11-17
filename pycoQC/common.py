# -*- coding: utf-8 -*-

# Standard library imports
from datetime import datetime as _datetime
from os import access, R_OK, W_OK, listdir, path, makedirs, walk
from pathlib import Path
import gzip
import inspect
from glob import glob
import logging
import re
import traceback
from collections import Counter, OrderedDict

# Third party imports
from dateutil.parser import parse as _dateutil_parse
import numpy as np
import pandas as pd
import pysam as ps


#~~~~~~~~~~~~~~CUSTOM EXCEPTION AND WARN CLASSES~~~~~~~~~~~~~~#
# Anything that is not a letter or a digit separates words in a path.
_TOKEN_SPLIT = re.compile(r"[^a-z0-9]+")


class pycoQCError(Exception):
    """ Basic exception class for pycoQC package """
    pass


class pycoQCWarning(Warning):
    """ Basic Warning class for pycoQC package """
    pass


##~~~~~~~ READ QUALITY ~~~~~~~#
# Error probability for every possible per-base quality value, precomputed once
# so that scoring a read is a single table lookup plus a mean.
#   _ERROR_PROB[q]        q is a raw Phred score (pysam query_qualities)
#   _ERROR_PROB_PHRED33[b] b is a byte of an ASCII Phred+33 string (FASTQ/SAM)
_ERROR_PROB = np.power(10.0, -np.arange(256, dtype=np.float64) / 10.0)
_ERROR_PROB_PHRED33 = np.power(10.0, -(np.arange(256, dtype=np.float64) - 33.0) / 10.0)


def mean_qscore(qualities, phred33=False):
    """Mean quality of a read, expressed as a Phred score

    Phred scores are logarithmic, so averaging them directly overstates the
    quality. The per-base error probabilities are averaged instead and the
    result converted back to a Phred score.

    * qualities
        Per-base qualities, either as raw Phred values (any buffer of bytes,
        such as the array returned by pysam's query_qualities) or as an ASCII
        Phred+33 encoded string when phred33 is True.
    * phred33
        Treat qualities as ASCII Phred+33 encoded rather than as raw scores.
    """
    if not qualities:
        return 0.0
    if isinstance(qualities, str):
        qualities = qualities.encode("ascii")
    codes = np.frombuffer(qualities, dtype=np.uint8)
    if codes.size == 0:
        return 0.0
    table = _ERROR_PROB_PHRED33 if phred33 else _ERROR_PROB
    # Counting the 256 possible quality values and taking a single dot product
    # beats gathering one probability per base, since reads are much longer
    # than the table.
    counts = np.bincount(codes, minlength=256)
    return float(-10.0 * np.log10(counts @ table / codes.size))


def path_has_token(path, token):
    """Is `token` a whole word anywhere in this path?

    The basecallers separate their output into fastq_pass and fastq_fail
    directories, and name files like FAW88668_pass_barcode06_0.pod5, so where
    a file sits says whether its reads passed. Testing the whole path for the
    substring "fail" misread any run stored under a directory such as
    failsafe_backup or /mnt/failover, marking every read as failing. Only
    whole words between the separators count.

    * path
        A file path.
    * token
        The word to look for, compared without regard to case.
    """
    words = _TOKEN_SPLIT.split(str(path).lower())
    return token.lower() in words


def path_below(path, root):
    """The part of `path` that lies under `root`.

    Pass and fail are read from where the basecaller put a file, and the
    basecaller only chose the part of the path below the directory it was
    given. Reading the whole absolute path meant that a run stored under a
    directory someone had named for a failed experiment came out as a run
    of failing reads. Where the two are unrelated the file name is used,
    which is never worse than the full path.

    * path
        A file path.
    * root
        The directory the converter was pointed at.
    """
    try:
        return str(Path(path).resolve().relative_to(Path(root).resolve()))
    except (ValueError, OSError):
        return str(Path(path).name)


def run_health_fields(tags):
    """Duration, template duration, end reason and mux from Dorado's tags.

    Dorado writes these tags into a BAM, and samtools fastq -T copies them into
    a FASTQ header:
    du, the read's duration in seconds; ns, its number of signal samples; ts,
    the samples trimmed from its start, which is the adapter; er, the reason
    it ended; and mx, the mux it was read on. The strand's own time is
    du * (ns - ts) / ns, which is what a translocation speed should divide by.
    Missing tags give blanks.

    * tags
        Mapping from tag name to value, as strings or numbers.
    """
    def number(key):
        try:
            return float(tags[key])
        except (KeyError, TypeError, ValueError):
            return None

    du, ns, ts = number("du"), number("ns"), number("ts")
    template = ""
    if du is not None and ns and ts is not None and 0 <= ts < ns:
        template = round(du * (ns - ts) / ns, 5)
    mux = number("mx")
    return {
        "duration": round(du, 5) if du is not None else "",
        "template_duration": template,
        "end_reason": tags.get("er", ""),
        "mux": int(mux) if mux is not None else "",
    }


def parse_timestamp(value):
    """Parse a basecaller timestamp

    Basecallers write ISO 8601 timestamps, which the standard library parses
    around 160 times faster than the general purpose dateutil parser. Anything
    it cannot read is handed to dateutil so unusual formats still work.

    * value
        Timestamp string, e.g. "2022-06-27T23:14:28Z".
    """
    try:
        # fromisoformat only accepts a numeric UTC offset before Python 3.11.
        return _datetime.fromisoformat(value.replace("Z", "+00:00"))
    except (ValueError, AttributeError, TypeError):
        return _dateutil_parse(value)


def start_time_to_elapsed_seconds(start_time):
    """Convert a column of read start times to seconds since the first read

    The result is returned as whole seconds where possible, falling back to
    floats when some times are missing.

    * start_time
        Series of parseable timestamps.
    """
    parsed = pd.to_datetime(start_time)
    elapsed = (parsed - parsed.min()).dt.total_seconds()
    try:
        return elapsed.astype("int64")
    except (pd.errors.IntCastingNaNError, ValueError):
        return elapsed


def elapsed_seconds_per_run(df, time_col="start_time", run_col="run_id"):
    """Convert read start times to seconds since the start of their own run.

    A basecaller's own summary counts from the moment its run started, and
    pycoQC relies on that when it lays several runs out in order. Counting
    from the first read in the directory instead puts the whole wall clock
    gap between two runs inside the second one: two one hour runs recorded
    two days apart were reported as a single 50 hour run.

    * df
        Frame holding the times, modified in place.
    * time_col, run_col
        Column names. Without the run column the whole frame is one run.
    """
    if time_col not in df:
        return df

    if run_col in df and df[run_col].notna().any():
        # groupby drops rows whose key is missing, so unlabelled reads are
        # gathered under one key of their own rather than being lost.
        keys = df[run_col].fillna("")
        df[time_col] = (df.groupby(keys, sort=False)[time_col]
                          .transform(start_time_to_elapsed_seconds))
    else:
        df[time_col] = start_time_to_elapsed_seconds(df[time_col])
    return df


##~~~~~~~ FUNCTIONS ~~~~~~~#
def as_labels(values):
    """Barcodes as names, without turning 6 into "6.0".

    A column of bare numbers, as Deepbinner writes, is typed as a number,
    and once a merge introduces a missing value it becomes a float, so the
    plain string of it carries a decimal point into the report, the plot
    and the name of the file Barcode_split writes.
    """
    if pd.api.types.is_float_dtype(values):
        whole = values.dropna() % 1 == 0
        if whole.all():
            values = values.astype("Int64")
    return values.astype("object").where(
        values.isna(), values.astype(str).str.strip())


def is_readable_file(fn):
    """Verify the readability of a file or list of file"""
    return path.isfile(fn) and access(fn, R_OK)


def check_arg(arg_name, arg_val, required_type, allow_none=True, min=None, max=None, choices=()):
    """Check an argument's type and range, coercing it where possible.

    Failures raise pycoQCError, the one exception the command line turns
    into a plain message. A bare Exception, which this used to raise, came
    out as a traceback, and `except:` around the coercion also swallowed
    KeyboardInterrupt.
    """
    if allow_none and arg_val is None:
        return arg_val

    if not isinstance(arg_val, required_type):
        try:
            arg_val = required_type(arg_val)
        except (TypeError, ValueError):
            raise pycoQCError("Argument `{}` value `{}` is not in correct type: `{}` and cannot be coerced".format(
                arg_name, arg_val, required_type.__name__))

    if required_type in [float, int]:
        # A bound of 0 is falsy, so testing truthiness skipped it entirely
        # and a negative min_pass_len or min_barcode_percent was accepted.
        if min is not None and arg_val < min:
            raise pycoQCError("Argument `{}` value `{}` is too low. Minimal value: {}".format(
                arg_name, arg_val, min))
        if max is not None and arg_val > max:
            raise pycoQCError("Argument `{}` value `{}` is too high. Maximal value: {}".format(
                arg_name, arg_val, max))

    if choices and arg_val not in choices:
        raise pycoQCError("Argument `{}` value `{}` is not in the list of possible choices. Choices: {}".format(
            arg_name, arg_val, ", ".join(choices)))

    return arg_val


def sequencing_summary_file_sample(infile, outfile=None, n_seq=10000):
    """
    Sample a number read lines in infile and write the output_over_time in output_file
    If the file contains several runids the function will sample proportionally to the
    * infile: STR
    Path to a sequencing_summary input file
    * outfile: STR(default None)
        Path to a sequencing_summary output file. If not given, will return a dataframe instead
    * n_seq: STR(default 10000)
        Overall number of sequence lines to sample
    """
    df = pd.read_csv(infile, sep ="\t")
    df = df.dropna()
    total = len(df)
    print("{} sequences".format(total))

    l = []
    for runid, runid_df in df.groupby("run_id", sort=False):
        n_to_sample = int(round(len(runid_df)/total*n_seq, 0))
        if n_to_sample == 0:
            n_to_sample=1
        print("{} = {} seq, to sample = {}".format(runid, len(runid_df), n_to_sample))
        sdf = runid_df.sample(n_to_sample).sort_values("start_time")
        l.append(sdf)

    df = pd.concat(l).reset_index(drop=True)
    if outfile:
        if outfile.endswith("gz"):
            df.to_csv(outfile, index=False, sep="\t", compression="gzip")
        else:
            df.to_csv(outfile, index=False, sep="\t", compression=None)
    else:
        return df


def dict_to_str(c, prefix="\t", suffix="\n"):
    """ Transform a dict to a tabulated str """
    m = ""
    if type(c) == Counter:
        for i, j in c.most_common():
            m += "{}{}: {:,}{}".format(prefix, i, j, suffix)
    else:
        for i, j in c.items():
            m += "{}{}: {}{}".format(prefix, i, j, suffix)
    return m


def recursive_file_gen(dir, ext):
    """List every file with the given extension under a directory.

    * dir
        Directory to walk. A path that is not a directory yields nothing.
    * ext
        Extension, without the dot, or a list of them.

    This used to stop descending as soon as one matching file was found in a
    directory, which quietly lost reads from a very ordinary layout: a
    fastq_pass folder holding unclassified.fastq beside one folder per
    barcode returned the unclassified reads and nothing else, with no
    warning and a read count that looked plausible.
    """
    extensions = ext if isinstance(ext, list) else [ext]
    suffixes = tuple("." + e for e in extensions)

    if not path.isdir(dir):
        return

    for root, dirnames, filenames in walk(dir):
        # Sorted, so converting a run twice gives the same summary twice.
        dirnames.sort()
        for fn in sorted(filenames):
            if fn.endswith(suffixes):
                yield path.join(root, fn)


# Converter verbosity, from --verbose_level.
VERBOSE_LEVELS = {2: logging.DEBUG, 1: logging.INFO, 0: logging.WARNING}


def format_worker_error(exc):
    """What a converter's worker process reports when it fails.

    A pycoQCError already says what is wrong with the input, so its message
    is sent alone; wrapped in the worker's traceback it reached the user as
    thirty lines of somebody else's stack with the reason at the bottom.
    Anything else is a bug, and its traceback is how it gets found.

    * exc
        The exception being handled.
    """
    if isinstance(exc, pycoQCError):
        return str(exc)
    return traceback.format_exc()


def get_logger(name=None, verbose=False, quiet=False, level=None):
    """A logger for one module, set to the requested verbosity.

    Progress is logged at INFO, so `quiet` leaves only warnings and errors.
    Output goes through the root logger, which is given a plain handler only
    if nothing has configured logging yet: importing a pycoQC module used to
    configure the root logger at INFO as a side effect, which turned on the
    debug chatter of every other library in the process.

    * name
        Logger name, normally the calling module's __name__.
    * verbose, quiet
        Debug output, or warnings only. Ignored when level is given.
    * level
        An explicit logging level.
    """
    logging.basicConfig(format='%(message)s')
    logger = logging.getLogger(name)

    if level is None:
        if verbose:
            level = logging.DEBUG
        elif quiet:
            level = logging.WARNING
        else:
            level = logging.INFO
    logger.setLevel(level)

    return logger


def check_converter_io(in_dir, out_fn, kind):
    """The checks every converter makes before it starts.

    * in_dir
        Directory to read, which has to exist: os.access alone accepted a
        plain file, which was then walked as a directory and reported as
        holding no reads.
    * out_fn
        Summary to write. Its directory has to be writable; a bare file
        name means the current directory.
    * kind
        Name of the input format, for the messages.
    """
    if not path.isdir(in_dir) or not access(in_dir, R_OK):
        raise pycoQCError("Cannot read the indicated {} directory: {}".format(kind, in_dir))
    if not access(path.dirname(path.abspath(out_fn)), W_OK):
        raise pycoQCError("Cannot write the indicated seq_summary_fn: {}".format(out_fn))


def doc_func(func):
    """Parse the function description string"""

    docstr_list = []
    for l in inspect.getdoc(func).split("\n"):
        l = l.strip()
        if l:
            if l.startswith("*"):
                break
            else:
                docstr_list.append(l)

    return " ".join(docstr_list)


def _help_defaults(cls):
    """The defaults a class shows where its signature has a stand in, as a
    converter's fields=None does for its DEFAULT_FIELDS. Only a classmethod
    called help_defaults is asked: anything else of that name, on a class
    that is no converter, is left alone."""
    hook = inspect.getattr_static(cls, "help_defaults", None)
    return cls.help_defaults() if isinstance(hook, classmethod) else {}


def _is_stand_in(default):
    """None or an empty list or tuple: a default a converter replaces with its
    own, as opposed to a real one."""
    return default is None or (isinstance(default, (list, tuple)) and len(default) == 0)


def make_arg_dict(func):
    """Parse the arguments default value, type and doc"""

    # Init method for classes. A converter's fields=None stands for its
    # DEFAULT_FIELDS, which the help shows for the class. Given an __init__
    # alone, the help shows its signature as it is: the function does not
    # say which class it is used by, and a subclass inheriting it has
    # defaults of its own.
    shown = {}
    if inspect.isclass(func):
        shown = _help_defaults(func)
        func = func.__init__

    if inspect.isfunction(func) or inspect.ismethod(func):
        # Parse arguments default values and annotations
        d = OrderedDict()
        for name, p in inspect.signature(func).parameters.items():
            if p.name not in ["self","cls"]: # Object stuff. Does not make sense to include in doc
                d[name] = OrderedDict()
                if name not in ["kwargs","args"]: # Include but skip default required and type
                    # Get Annotation
                    if p.annotation != inspect._empty:
                        d[name]["type"] = p.annotation
                    # Get default value if available
                    if p.default == inspect._empty:
                        d[name]["required"] = True
                    else:
                        # Only in place of the stand in, never of a real default.
                        d[name]["default"] = (shown[name] if name in shown and _is_stand_in(p.default)
                                              else p.default)

        # Parse the docstring in a dict
        docstr_dict = OrderedDict()
        lab=None
        for l in inspect.getdoc(func).split("\n"):
            l = l.strip()
            if l:
                if l.startswith("*"):
                    lab = l[1:].strip()
                    docstr_dict[lab] = []
                elif lab:
                    docstr_dict[lab].append(l)

        # Concatenate and copy doc in main dict
        for name in d.keys():
            if name in docstr_dict:
                d[name]["help"] = " ".join(docstr_dict[name])
        return d


def arg_opt(func, arg, **kwargs):
    """Get options corresponding to argument name and deal with special cases"""
    arg_dict = make_arg_dict(func)[arg]

    if "default" in arg_dict and "help" in arg_dict:
        arg_dict["help"] += "(default: %(default)s)"

    if "type" in arg_dict and "help" in arg_dict:
        arg_dict["help"] += " [%(type)s]"

    # Special case for boolean args
    if arg_dict["type"] == bool:
        if arg_dict["default"] is False:
            arg_dict["action"] = 'store_true'
            del arg_dict["type"]
        elif arg_dict["default"] is True:
            arg_dict["action"] = 'store_false'
            del arg_dict["type"]

    # Special case for lists args
    elif arg_dict["type"] == list:
        arg_dict["nargs"]='*'

    return arg_dict


def jhelp(f):
    """
    Display a Markdown pretty help message for functions and class methods(default __init__ is a class is passed)

    * f
        A python function or method.
    jhelp also display default values and type annotations if available.
    The docstring synthax should follow the same synthax as the one used for this function
    * f
        Function or method to display the help message for
    """
    # Private import as this is only needed if using jupyter
    from IPython.display import display, Markdown

    f_doc = doc_func(f)
    arg_doc = make_arg_dict(f)

    # Signature and function documentation
    s = "**{}**({})\n\n{}\n\n---\n\n".format(f.__name__, ", ".join(arg_doc.keys()), f_doc)

    # Args doc
    for arg_name, arg_val in arg_doc.items():
        # Arg signature section
        s+= "* **{}**".format(arg_name)
        if "default" in arg_val:
            if arg_val["default"] == "":
                  arg_val["default"] = "\"\""
            s+= "(default: {})".format(arg_val["default"])
        if "required" in arg_val:
            s+= "(required)"
        if "type" in arg_val:
            # An annotation need not be a class, and only classes have a
            # __name__. The test guarding this compared type(list) with type,
            # which is always true, so a string annotation raised.
            s+= " [{}]".format(getattr(arg_val["type"], "__name__", arg_val["type"]))
        s+="\n\n"
        # Arg doc section
        if "help" in arg_val:
            s+= "{}\n\n".format(arg_val["help"])

    # Display in Jupyter
    display(Markdown(s))


def is_gzip_name(fn):
    """Whether a file is read as gzip, which is decided by its name."""
    return str(fn).endswith(".gz")


def head(fp, n=10, sep="\t", comment=None):
    """
    Emulate linux head cmd, for plain or gzipped text files
    * fp
        Path to the file to be parsed.
    * n
        Number of lines to print starting from the beginning of the file(Default 10)
    """
    line_list = []

    # Get lines. Summaries are routinely kept gzipped, and opening one as
    # text raised UnicodeDecodeError rather than showing it.
    opener = gzip.open if is_gzip_name(fp) else open
    try:
        with opener(fp, "rt") as fh:
            line_num = 0
            while(line_num < n):
                l= next(fh).strip()
                if comment and l.startswith(comment):
                    continue
                if sep:
                    line_list.append(l.split(sep))
                else:
                    line_list.append(l)
                line_num+=1

    except StopIteration:
        pass

    # Add padding if sep given
    if sep:
        try:
            # Find longest elem per col
            col_len_list = [0 for _ in range(len(line_list[0]))]
            for ls in line_list:
                for i in range(len(ls)):
                    len_col = len(ls[i])
                    if len_col > col_len_list[i]:
                        col_len_list[i] = len_col

            # Add padding
            line_list_tab = []
            for ls in line_list:
                s = ""
                for i in range(len(ls)):
                    len_col = col_len_list[i]
                    len_cur_col = len(ls[i])
                    s += ls[i][0:len_col] + " "*(len_col-len_cur_col)+" "
                line_list_tab.append(s)
            line_list = line_list_tab

        # Fall back to non tabulated display
        except IndexError:
            return head(fp=fp, n=n, sep=None)

    for l in line_list:
        print(l)
    print()


def ls(dir_path):
    for f in listdir(dir_path):
        print(f)


def expand_file_names(fn, bam_check=False):
    """Turn file names and glob patterns into a list of readable files.

    * fn
        A path or glob pattern, or a list of them.
    * bam_check
        Also require each file to be a sorted, indexed BAM.

    A name that exists is taken literally. glob reads square brackets as a
    character class, so a summary called run[1].txt could not be named at
    all. Matches are sorted, since glob returns them in directory order and
    the order of the files decides the order of runs with no start times.
    A file named twice, or matched by two patterns, is read once: reading it
    twice doubled every count in the report.
    """
    if isinstance(fn, str):
        patterns = [fn]
    elif isinstance(fn, (list, tuple)):
        patterns = list(fn)
    else:
        raise pycoQCError("{} has to be either a file or a regular expression or a list of files".format(fn))

    fn_list = []
    seen = set()
    for pattern in patterns:
        matches = [pattern] if path.isfile(pattern) else sorted(glob(pattern))
        for match in matches:
            key = path.abspath(match)
            if key not in seen:
                seen.add(key)
                fn_list.append(match)

    # Verify that files are readable
    if not fn_list:
        raise pycoQCError("No files found in {}".format(fn))
    for f in fn_list:
        if not is_readable_file(f):
            raise pycoQCError("Cannot read file {}".format(f))
        # Extra checks for bam files
        if bam_check:
            with ps.AlignmentFile(f, "rb") as bam:
                if not bam.has_index():
                    raise pycoQCError("No index found for bam file: {}. Please index with samtools index".format(f))
                # A BAM may carry no @HD line at all, which raised a bare
                # KeyError instead of saying what was wrong with the file.
                if bam.header.to_dict().get("HD", {}).get("SO") != "coordinate":
                    raise pycoQCError("Bam file not sorted: {}. Please sort with samtools sort".format(f))
    return fn_list


def merge_files_to_df(fn_list, logger=None):
    """Read tab separated files and stack them into one frame.

    * fn_list
        Paths to read.
    * logger
        Told which columns the files disagree on, if any.
    """
    try:
        if len(fn_list) == 1:
            df = pd.read_csv(fn_list[0], sep="\t")

        else:
            df_list = []
            for fn in fn_list:
                df_list.append(pd.read_csv(fn, sep="\t"))

            # Only the columns every file has survive the stack. That is the
            # safe choice, but it used to happen in silence: two summaries
            # where one lacked start_time produced a report with no time
            # plots at all and nothing to say why.
            shared = set(df_list[0].columns).intersection(
                *[set(d.columns) for d in df_list[1:]])
            dropped = sorted(set().union(*[set(d.columns) for d in df_list])
                             - shared)
            if dropped and logger:
                logger.warning(
                    "WARNING: {} not present in every input file and left "
                    "out: {}".format(
                        "columns are" if len(dropped) > 1 else "column is",
                        ", ".join(dropped)))

            df = pd.concat(df_list, ignore_index=True, sort=False, join="inner")
    except pd.errors.EmptyDataError:
        raise pycoQCError(
            "No data found in {}. The file is empty, or holds a header with "
            "no reads under it.".format(
                fn_list[0] if len(fn_list) == 1 else ", ".join(fn_list)))

    if len(df) == 0:
        raise pycoQCError("No valid read found in input file")

    return df


def mkdir(fn, exist_ok=False):
    """ Create directory recursivelly. Raise IO error if path exist or if error at creation """
    try:
        makedirs(fn, exist_ok=exist_ok)
    except OSError as E:
        # This used to raise pycoMethError, a name from another project that
        # does not exist here, so the real failure surfaced as a NameError.
        raise pycoQCError("Error creating output folder `{}`: {}".format(fn, E))


def mkbasedir(fn, exist_ok=False):
    """ Create directory for a given file recursively. Raise IO error if path exist or if error at creation """
    dir_fn = path.dirname(fn)
    if dir_fn:
        mkdir(dir_fn, exist_ok=True)
