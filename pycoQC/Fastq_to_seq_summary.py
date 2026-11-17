#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import os
from collections import OrderedDict
import logging
from itertools import islice
import gzip

# Third party imports
import pandas as pd

# Local imports
from pycoQC._pipeline import ConverterPipeline
from pycoQC.common import (
    is_gzip_name,
    path_below,
    path_has_token,
    mean_qscore,
    parse_timestamp,
    run_health_fields,
    elapsed_seconds_per_run)


# Configured when a conversion starts, not on import: setting up the root
# logger here switched on INFO output for every library in the process.
logger = logging.getLogger(__name__)


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
    def is_gzipped(fastq_fn):
        """Whether a fastq file is read as gzip, by its name."""
        return is_gzip_name(fastq_fn)

    @staticmethod
    def open_fastq(fastq_fn):
        """Open a fastq file, transparently handling gzip compression"""
        if FastqParser.is_gzipped(fastq_fn):
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
        # Present when the header carries Dorado's tags, as samtools fastq -T writes them.
        entry_dict[read_id].update(run_health_fields(read_dict))

        return entry_dict


class Fastq_to_seq_summary(ConverterPipeline):
    """
    Create a summary file akin the one generated by Albacore or Guppy from a directory containing
    multiple fastq files. The script will attempt to extract all the required fields but will not
    raise an error if not found.
    """

    kind = "fastq"
    extensions = ("fq", "fq.gz", "fastq", "fastq.gz")

    # Number of reads parsed as one block. Counting in reads rather than
    # lines means a block can never end in the middle of a record, which
    # would drop that read.
    READS_PER_CHUNK = 1000

    # How far over an even share of the input a file may be and still be
    # handed to one parser whole.
    WHOLE_FILE_SLACK = 1.25

    # How many times larger a gzipped FASTQ is once decompressed, to weigh
    # it against plain files: 2.2 to 2.9 on Dorado output, whose quality
    # strings compress poorly.
    GZIP_RATIO = 2.5

    # The fields written when none are named, here and on the command line,
    # which reads its --fields default from this. A subclass changes its
    # default here. The command line once kept a list of its own, and calling
    # this class directly silently dropped the pass or fail column that the
    # command produced.
    DEFAULT_FIELDS = ("read_id", "run_id", "channel", "start_time", "gc_percent",
                      "sequence_length_template", "mean_qscore_template",
                      "calibration_strand_genome_template", "barcode_arrangement",
                      "passes_filtering", "duration", "template_duration",
                      "end_reason", "mux")

    def __init__(self,
                 fastq_dir: str,
                 seq_summary_fn: str,
                 max_fastq: int = 0,
                 threads: int = 4,
                 basecall_id: int = 0,
                 verbose_level: int = 0,
                 include_path: bool = False,
                 fields: list = None):
        """
        * fastq_dir
            Directory containing fastqfiles. Can contain multiple subdirectories
        * seq_summary_fn
            path of the summary sequencing file where to write the data extracted from the fastq files
        * max_fastq
            Maximum number of file to try to parse. 0 to deactivate
        * threads
            Total number of threads to use. 1 thread is used for the reader and 1 for the writer. Minimum 3, and at most 31 on Windows before
            Python 3.13 (default = 4)
        * fields
            Fields to fetch from the fastq files. None or an empty list, the
            default, is DEFAULT_FIELDS. Valid field names:
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
        self._configure(fastq_dir, seq_summary_fn, max_fastq, threads,
                        basecall_id, include_path, verbose_level)
        self.fields = self._requested_fields(fields)

        # Blocks of lines are large, so few are kept waiting for a parser.
        self.in_queue_size = max(4, self.threads * 2)

        # Captured on the instance so the workers, which get a pickled copy
        # of it rather than the parent's class, see a value changed on the
        # class before the call.
        self.reads_per_chunk = self.READS_PER_CHUNK

        self._run()

    #~~~~~~~~~~~~~~PIPELINE STAGES~~~~~~~~~~~~~~#
    def _work_items(self):
        """Files small enough to be one parser's share whole, blocks of reads if not.

        A parser handed a file path reads and decompresses it itself, so that
        work is spread over every parser, and only the parsed records cross
        between processes. A file larger than an even share of the input
        would keep one parser busy while the others ran out of work, as one
        large pass file among small fail files did, so the lister splits
        such a file into blocks of reads for every parser to take from.
        Each item is (path, flag, lines); a whole file has None for its lines.
        """
        files = list(self._input_files())
        # The work in a file is its reads, so a gzipped file is weighed by
        # roughly what it holds rather than what it takes on disk.
        sizes = [os.path.getsize(fn) * (self.GZIP_RATIO if FastqParser.is_gzipped(fn) else 1)
                 for fn in files]
        # A quarter over an even share is still taken whole, so files of
        # about the same size are not split for being a little over average.
        # That leaves a parser idle when there are a few files fewer than
        # parsers, and it is the cheaper choice: a split file is read and
        # decompressed by the lister alone, and on a 180,000 read benchmark
        # that one process took 8.9 s for the whole input, against about
        # 4.7 s for a parser to take a quarter of it whole.
        limit = self.WHOLE_FILE_SLACK * sum(sizes) / max(1, self.threads)
        for fastq_fn, size in zip(files, sizes):
            flag = FastqParser.pass_fail_flag(fastq_fn, self.in_dir)
            if size <= limit:
                yield fastq_fn, flag, None
                continue
            for chunk in self._chunks(fastq_fn):
                yield fastq_fn, flag, chunk

    def _item_file(self, item):
        return item[0]

    def _chunks(self, fastq_fn):
        # Counted in reads rather than lines, so a block never ends in the
        # middle of a record.
        with FastqParser.open_fastq(fastq_fn) as fastq_fh:
            yield from FastqParser.make_chunks(
                fastq_fh, self.reads_per_chunk * FastqParser.LINES_PER_READ)

    def _records(self, item, counts):
        fastq_fn, flag, chunk = item
        chunks = self._chunks(fastq_fn) if chunk is None else (chunk,)
        path = os.path.abspath(fastq_fn) if self.include_path else None
        for chunk in chunks:
            for read_info_dict in FastqParser.read_entry(chunk, flag):
                d = OrderedDict()
                for field in self.fields:
                    if field in read_info_dict:
                        d[field] = read_info_dict[field]
                    # A blank is written as one, but counted as not found:
                    # the run health fields are always present, blank when
                    # the header carries no Dorado tags, and counting those
                    # as found hid that they were missing.
                    if read_info_dict.get(field, "") not in ("", None):
                        counts["fields_found"][field] += 1
                    else:
                        counts["fields_not_found"][field] += 1
                if path is not None:
                    d["path"] = path
                yield d

    def _finalise(self, df, notes):
        # start_time is only present if it was asked for, and --fields lets
        # a caller leave it out, which used to raise KeyError once every
        # file had already been parsed.
        if 'start_time' in df:
            # Keep the absolute start of the run before the per read times
            # become seconds elapsed, so the report can say when it ran.
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

            # Express read start times as seconds since their own run began,
            # which is what the basecaller's summary does.
            elapsed_seconds_per_run(df)
        return df
