# Disable multithreading for MKL and openBlas
import os
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["MKL_THREADING_LAYER"] = "sequential"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["OMP_NUM_THREADS"] = "1"
os.environ['OPENBLAS_NUM_THREADS'] = '1'

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import multiprocessing as mp
from time import time
from collections import Counter, OrderedDict
import traceback
import logging
from itertools import islice
import gzip

# Third party imports
import pandas as pd
from tqdm import tqdm

# Local imports
from pycoQC.common import (
    dict_to_str,
    path_below,
    path_has_token,
    mean_qscore,
    parse_timestamp,
    pycoQCError,
    recursive_file_gen,
    stop_processes,
    elapsed_seconds_per_run)


# Logger setup
logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)
logLevel_dict = {2: logging.DEBUG, 1: logging.INFO, 0: logging.WARNING}


#~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~#
class FastqParser(object):

    # A fastq record is four lines: header, sequence, separator, quality.
    LINES_PER_READ = 4

    @staticmethod
    def make_chunks(file_handle, size):
        while True:
            chunk = list(islice(file_handle, size))
            if not chunk:
                break
            yield chunk

    @staticmethod
    def open_fastq(fastq_fn):
        """Open a fastq file, transparently handling gzip compression"""
        if fastq_fn.endswith('.gz'):
            return gzip.open(fastq_fn, "rt")
        return open(fastq_fn, "r")

    @staticmethod
    def pass_fail_flag(fastq_fn, root=None):
        """Basecallers write passing and failing reads to separate paths

        Matched on whole words, and only below the directory the converter
        was given, since that is the only part the basecaller chose. A run
        kept under a directory such as failsafe_backup, or one called after
        a failed experiment, is not a directory of failing reads.
        """
        location = path_below(fastq_fn, root) if root else fastq_fn
        return 'fail' if path_has_token(location, 'fail') else 'pass'

    @staticmethod
    def read_entry(entry, flag):
        """Parse a block of fastq lines into {read_id: fields}

        Blocks are cut on read boundaries, so leftover lines mean the file
        itself ended mid-record. That read is skipped, with a warning, rather
        than silently dropped.
        """
        # A list, not a dict keyed by read id. Keying by id silently merged
        # repeated ids, but only when they fell in the same block, so the
        # read count depended on where the block boundaries happened to
        # land. Every record is kept and pycoQC's --filter_duplicated
        # decides what to do about repeats.
        my_dict = list()
        lines = list()
        for line in entry:
            line = line.rstrip()
            lines.append(line)
            if len(lines) == FastqParser.LINES_PER_READ:
                my_dict.extend(
                    FastqParser.single_fastq_entry_to_dict(lines, flag).values())
                lines = list()
        if lines:
            logger.warning(
                "Truncated fastq record at end of block: dropping {} trailing "
                "line(s). The file may be incomplete.".format(len(lines)))
        return my_dict

    @staticmethod
    def parse_header_fields(items):
        """Read the metadata attached to a fastq header into a dict

        Two dialects are in use and both appear in the wild:

        * Guppy and older MinKNOW write `key=value`::

              @20a271c3-8f4b runid=acfd15dc sampleid=run3 read=2743 ch=207
              start_time=2022-06-27T23:14:28Z barcode=barcode17

        * Dorado and current MinKNOW write SAM style `attr:type:value`, which
          may be separated by spaces or tabs::

              @bd8655fb-383c RG:Z:e4994c62_dna_r10.4.1_sup@v5.0.0 ch:i:123
              st:Z:2025-01-06T10:06:36.778368+00:00 du:f:1.23 qs:f:12.5

        Unrecognised items are ignored rather than raising, so an unexpected
        field cannot stop a file being parsed.

        * items
            The whitespace-separated pieces of the header, excluding the read id.
        """
        fields = dict()
        for item in items:
            # A SAM style tag is two letters, a one letter type and a value.
            # Check it before `=`, since a value may itself contain one.
            tag, sep, rest = item.partition(':')
            if sep and len(tag) == 2:
                tag_type, sep, value = rest.partition(':')
                if sep and len(tag_type) == 1:
                    fields[tag] = value
                    continue
            key, sep, value = item.partition('=')
            if sep:
                fields[key] = value
        return fields

    @staticmethod
    def single_fastq_entry_to_dict(fastq_entry, flag):
        """Turn one fastq record into {read_id: fields}"""

        header, seq, extra, qual = fastq_entry  # get each component of list in a variable

        # Split header line to get items and store in dictionary
        items = header.split()
        read_dict = FastqParser.parse_header_fields(items[1:])

        # Fields are read under either dialect's name, and anything missing is
        # left blank rather than failing: a plain fastq carries none of them.
        # Dorado has no run id tag, but its read group starts with one.
        run_id = read_dict.get('runid', '')
        if not run_id and 'RG' in read_dict:
            run_id = read_dict['RG'].split('_')[0]

        start_time = read_dict.get('start_time') or read_dict.get('st')
        time_string = parse_timestamp(start_time) if start_time else ''

        channel = read_dict.get('ch', '')

        # Dorado writes the flowcell as PU and the sample library as LB.
        flowcell = read_dict.get('PU', '')
        sample_id = read_dict.get('sampleid') or read_dict.get('LB') or ''

        # Guppy names the barcode, Dorado uses SM with al as the alias, and
        # some versions write BC.
        barcode = (read_dict.get('barcode') or read_dict.get('SM')
                   or read_dict.get('BC') or read_dict.get('al') or '')

        # The following values are always present in any given fastq
        read_id = items[0][1:]  # remove the leading "@"
        length = len(seq)

        # Average phred score
        average_phred = round(mean_qscore(qual, phred33=True), 2)

        # %GC. A zero length record is valid fastq, and dividing by its
        # length killed the conversion of every remaining file.
        if length:
            gc = round((seq.count('G') + seq.count('C')) / float(length) * 100, 2)
        else:
            gc = 0.0

        # Flag
        if flag == 'pass':
            flag = 'TRUE'
        else:
            flag = 'FALSE'

        # Output everyting in a dictionary
        entry_dict = dict()
        # entry_dict[read_id] = (run_id, time_string, channel, barcode, average_phred, length, gc)
        entry_dict[read_id] = dict()
        entry_dict[read_id]['read_id'] = read_id
        entry_dict[read_id]['run_id'] = run_id
        entry_dict[read_id]['channel'] = channel
        entry_dict[read_id]['barcode_arrangement'] = barcode
        entry_dict[read_id]['sequence_length_template'] = length
        entry_dict[read_id]['mean_qscore_template'] = average_phred
        entry_dict[read_id]['gc_percent'] = gc
        entry_dict[read_id]['start_time'] = time_string
        entry_dict[read_id]['passes_filtering'] = flag
        entry_dict[read_id]['flow_cell_id'] = flowcell
        entry_dict[read_id]['sample_id'] = sample_id

        return entry_dict


class Fastq_to_seq_summary():
    """
    Create a summary file akin the one generated by Albacore or Guppy from a directory containing
    multiple fastq files. The script will attempt to extract all the required fields but will not
    raise an error if not found.
    """

    # Number of reads handed to a worker at a time. Counting in reads rather
    # than lines means a block can never end in the middle of a record, which
    # would drop that read.
    READS_PER_CHUNK = 1000

    
    def __init__(self,
                 fastq_dir: str,
                 seq_summary_fn: str,
                 max_fastq: int = 0,
                 threads: int = 4,
                 basecall_id: int = 0,
                 verbose_level: int = 0,
                 include_path: bool = False,
                 # Matches the --fields default of the command line. The two
                 # disagreed, so calling this class directly silently dropped
                 # the pass or fail column that the command produced.
                 fields: list = ("read_id", "run_id", "channel", "start_time", "gc_percent",
                                 "sequence_length_template", "mean_qscore_template",
                                 "calibration_strand_genome_template", "barcode_arrangement",
                                 "passes_filtering")):
        """
        * fastq_dir
            Directory containing fastqfiles. Can contain multiple subdirectories
        * seq_summary_fn
            path of the summary sequencing file where to write the data extracted from the fastq files
        * max_fastq
            Maximum number of file to try to parse. 0 to deactivate
        * threads
            Total number of threads to use. 1 thread is used for the reader and 1 for the writer. Minimum 3(default = 4)
        * fields
            list of field names corresponding to attributes to try to fetch from the fastq files. List a valid field names:
            mean_qscore_template, sequence_length_template, called_events, skip_prob, stay_prob, step_prob, strand_score, read_id, start_time,
            duration, start_mux, read_number, channel, channel_digitisation, channel_offset, channel_range, channel_sampling,
            run_id, sample_id, device_id, protocol_run, flow_cell, calibration_strand, calibration_strand, calibration_strand,
            calibration_strand, barcode_arrangement, barcode_full, barcode_score, gc_percent
        * basecall_id
            id of the basecalling group. By default, leave to 0, but if you perform multiple basecalling on the same fastq files,
            this can be used to indicate the corresponding group(1, 2 ...)
        * include_path
            If True the absolute path to the corresponding file is added in an extra column
        * verbose_level
            Level of verbosity, from 2(Chatty) to 0(Nothing)
        """
        # Set logging level
        logger.setLevel(logLevel_dict.get(verbose_level, logging.WARNING))

        # Perform checks
        logger.info("Check input data and options")
        if not os.access(fastq_dir, os.R_OK):
            raise pycoQCError("Cannot read the indicated fastq directory")
        # dirname of a bare name is "", which is never writable, so
        # a plain output filename was rejected.
        if not os.access(os.path.dirname(os.path.abspath(seq_summary_fn)), os.W_OK):
            raise pycoQCError("Cannot write the indicated seq_summary_fn")
        if threads < 3:
            raise pycoQCError("At least 3 threads required")

        # Save self args
        self.fastq_dir = fastq_dir
        self.seq_summary_fn = seq_summary_fn
        self.threads = threads-2
        self.max_fastq = max_fastq
        self.fields = fields
        self.basecall_id = basecall_id
        self.include_path = include_path
        self.verbose_level = verbose_level

        # Init Multiprocessing variables. in_q carries blocks of fastq lines
        # rather than file paths, so it is kept short to bound how much of a
        # file is held in memory at once.
        in_q = mp.Queue(maxsize=max(4, self.threads * 2))
        out_q = mp.Queue(maxsize=1000)
        error_q = mp.Queue()
        counter_q = mp.Queue()

        # Define processes
        ps_list = [mp.Process(target=self._list_fastq, args=(in_q, error_q))]
        for i in range(self.threads):
            ps_list.append(mp.Process(target=self._read_fastq, args=(in_q, out_q, error_q, counter_q, i)))
        ps_list.append(mp.Process(target=self._write_seq_summary, args=(out_q, error_q, counter_q)))

        logger.info("Start processing fastq files")
        try:
            # Start all processes
            for proc in ps_list:
                proc.start()
            # Monitor error queue
            for tb in iter(error_q.get, None):
                raise pycoQCError(tb)
            # Join processes
            for proc in ps_list:
                proc.join()

        except BaseException:
            logger.info("\nAn error occurred. All processes were killed\n")
            raise

        finally:
            # Every exit from here, including a clean one, leaves nothing
            # running. Listing the exceptions to clean up after meant any
            # other kind left the whole pool alive.
            stop_processes(ps_list)

    def _list_fastq(self, in_q, error_q):
        """
        Mono-threaded worker walking the fastq files found in a directory tree
        recursively and feeding blocks of lines to the parsing workers

        Blocks of lines are queued rather than whole file paths. Each worker
        used to open a process pool of its own to parse one file, so the
        process count was the thread count squared; queueing chunks spreads a
        single file across all workers without that nesting.
        """
        logger.debug("[READER] Start reading fastq files")
        try:
            # Load an input queue with blocks of fastq lines
            my_ext = ['fq', 'fq.gz', 'fastq', 'fastq.gz']
            i = 0
            for fastq_fn in recursive_file_gen(dir=self.fastq_dir, ext=my_ext):
                flag = FastqParser.pass_fail_flag(fastq_fn, self.fastq_dir)
                with FastqParser.open_fastq(fastq_fn) as fastq_fh:
                    chunk_lines = self.READS_PER_CHUNK * FastqParser.LINES_PER_READ
                    for chunk in FastqParser.make_chunks(fastq_fh, chunk_lines):
                        in_q.put((chunk, fastq_fn, flag))
                i += 1
                if self.max_fastq and i >= self.max_fastq:
                    break

            # Raise error is no file found
            if i == 0:
                raise pycoQCError("No valid fastq files found in indicated folder")

            logger.debug("[READER] Read a total of {} files".format(i))

        # Manage exceptions and deal poison pills
        except Exception:
            error_q.put(traceback.format_exc())
        finally:
            for i in range(self.threads):
                in_q.put(None)

    def _read_fastq(self, in_q, out_q, error_q, counter_q, worker_id):
        """
        Multi-threaded workers in charge of parsing blocks of fastq lines.
        """
        logger.debug("[WORKER_{:02}] Start parsing fastq chunks".format(worker_id))
        try:
            c = {"overall": Counter(), "fields_found": Counter(), "fields_not_found": Counter()}

            # Parse each block of lines handed over by the reader
            for chunk, fastq_fn, flag in iter(in_q.get, None):
                # Fetch required fields is available
                """
                "read_id", "run_id", "channel", "start_time",
                "sequence_length_template", "mean_qscore_template", "gc_percent"
                "calibration_strand_genome_template", "barcode_arrangement"
                """
                for read_info_dict in FastqParser.read_entry(chunk, flag):
                    # A fresh dict per read: the queue serialises entries in a
                    # background thread, so a shared dict would be mutated
                    # before it is pickled.
                    d = OrderedDict()
                    for field in self.fields:
                        if field in read_info_dict:
                            d[field] = read_info_dict[field]
                            c["fields_found"][field] += 1
                        else:
                            c["fields_not_found"][field] += 1

                    if self.include_path:
                        d["path"] = os.path.abspath(fastq_fn)

                    # Put read data in queue
                    if d:
                        out_q.put(d)
                        c["overall"]["valid reads"] += 1
                    else:
                        c["overall"]["invalid reads"] += 1

            # Put counter in counter queue
            counter_q.put(c)

        # Manage exceptions and deal poison pills
        except Exception:
            error_q.put(traceback.format_exc())
        finally:
            out_q.put(None)

    def _write_seq_summary(self, out_q, error_q, counter_q):
        """
        Mono-threaded Worker writing the sequencing summary file
        """
        logger.debug("[WRITER] Start collecting summary data")

        t = time()
        try:
            l = []
            with tqdm(unit=" reads", mininterval=0.1, smoothing=0.1, disable=self.verbose_level == 2) as pbar:
                # Collect line data
                for _ in range(self.threads):
                    for d in iter(out_q.get, None):
                        l.append(d)
                        pbar.update(1)

            # Transform collected data to dataframe and write to file
            logger.debug("[WRITER] Write data to file")
            if not l:
                raise pycoQCError(
                    "No reads found in the fastq file(s). The files were read "
                    "but hold no records.")
            df = pd.DataFrame(l)

            # start_time is only present if it was asked for, and --fields
            # lets a caller leave it out, which used to raise KeyError once
            # every file had already been parsed.
            if 'start_time' in df:
                # Keep the absolute start of the run before the per read
                # times become seconds elapsed, so the report can say when
                # it ran.
                absolute = pd.to_datetime(df['start_time'], errors="coerce")
                if absolute.notna().any():
                    if 'run_id' in df and df['run_id'].notna().any():
                        # One directory can hold more than one run, and each
                        # one started when it started.
                        df['run_start_time'] = (
                            absolute.groupby(df['run_id'].fillna(""), sort=False)
                                    .transform("min").dt.strftime("%Y-%m-%dT%H:%M:%S"))
                    else:
                        df['run_start_time'] = absolute.min().isoformat()

                # Express read start times as seconds since their own run
                # began, which is what the basecaller's summary does.
                elapsed_seconds_per_run(df)


            df.to_csv(self.seq_summary_fn, sep="\t", index=False)

            # Collapse data from counters coming from workers
            logger.debug("[WRITER] Summarize counters")
            c = {"overall": Counter(), "fields_found": Counter(), "fields_not_found": Counter()}
            for _ in range(self.threads):
                worker_c = counter_q.get()
                for k, v in worker_c.items():
                    for k2, v2 in v.items():
                        c[k][k2] += v2

            logger.info("Overall counts {}".format(dict_to_str(c["overall"])))
            logger.info("fields found {}".format(dict_to_str(c["fields_found"])))
            logger.info("fields not found {}".format(dict_to_str(c["fields_not_found"])))

            # Print final
            logger.warning("Total reads: {} / Average speed: {} reads/s\n".format(len(df),
                                                                                  round(len(df) / (time() - t), 2)))

        # Manage exceptions and deal poison pills
        except Exception:
            error_q.put(traceback.format_exc())
        finally:
            error_q.put(None)
