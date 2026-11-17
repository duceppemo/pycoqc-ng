# Changelog

### Unreleased

Every entry point was run against real data: the original Albacore and
Guppy summaries, a MinKNOW R10.4.1 summary, a 1.25 million read Dorado
summary with a minimap2 BAM, Dorado BAM and FASTQ demultiplexed by barcode,
a PromethION POD5 and MinKNOW multi-read FAST5 files. What broke is below,
along with a restructuring of the plotting and report code that most of the
fixes fell out of.

**Added**

* Three plots on how the run went, rather than what it produced:
  * `end_reason_over_time`: the share of reads ending each way over the run.
    A rising share of `unblock_mux_change` means pores are clogging; on a 72
    hour Dorado run it went from 4.6% to 14.8%. On an adaptive sampling run,
    `data_service_unblock_mux_change` is the rejection rate.
  * `read_speed_1D` and `read_speed_over_time`: translocation speed in bases
    per second, measured on the strand's own time so short reads are not
    dragged down by the adapter. The JSON report carries its median and
    percentiles.
  * `active_channels_over_time`: channels producing reads per 1.5 hour
    interval, MinKNOW's mux scan cycle, in total and per mux, against the
    flow cell's channel count.

  All three are in the default report, and the two that need no basecalls in
  the report for POD5 and raw FAST5. The JSON report counts the end reasons.
* The converters carry what these need. `Bam_to_seq_summary` and
  `Fastq_to_seq_summary` write `duration`, `template_duration`, `end_reason`
  and `mux` from Dorado's `du`, `ns`, `ts`, `er` and `mx` tags (for FASTQ,
  when its headers carry those tags). `Pod5_to_seq_summary` writes `end_reason` and `mux`,
  and `Fast5_to_seq_summary` adds `start_mux` to its default fields. Summaries
  from Guppy, Dorado and MinKNOW already carry these columns, and pycoQC now
  keeps them rather than discarding them on reading.

**Fixed**

* Multi-read FAST5 files written by MinKNOW before it added the `file_type`
  attribute (`file_version` 1.0 and 2.0) were taken for single read files,
  and `Fast5_to_seq_summary` died looking for `/Raw/Reads`. They are now
  recognised by their layout.
* `pycoQC -v` on a run that has not been basecalled, a POD5 or MinKNOW
  summary, printed a traceback in the middle of the run.
* On a run with no failing reads, or no passing ones, the *Fail Reads* button
  of the density, over time and channel plots changed nothing, so the plot
  went on showing the previous read set under the new label.
* `--quiet` printed every progress step, which were logged as warnings.
* `runid_list` was ignored for a summary without start times, and is now
  also available on the command line as `--runid_list`, as it was in the
  original pycoQC. A single run id given as a string in Python is no longer
  read as a list of its characters.
* Barcode names with numbers in them are sorted as numbers everywhere,
  so `barcode2` comes before `barcode10`; the per barcode boxplots also
  sorted `barcode10` before `barcode2` on their own.
* `Barcode_split` named its files `sequencing_summary_6.0.txt` for a summary
  of numeric barcodes with a blank among them, and a barcode name holding a
  `/` wrote outside the output directory.
* The coverage plot drew a rectangle and a label for every reference, so an
  assembly with 233,000 contigs made a 31 MB plot that took a minute and a
  half. Past 100 references they are left off. With an even number of
  references the last one was never shaded.
* The channel activity map raised `IndexError` for a channel numbered past
  the platform's count, and drew channel 0's reads on the last channel.
* A summary named twice on the command line, or matched by two patterns,
  was read twice and every count doubled. A file whose name contains
  square brackets could not be named at all. Matches of a pattern are
  sorted, so runs without start times come in a stable order.
* A barcode file in none of the layouts pycoQC reads (Guppy, Deepbinner, or
  plain `read_id` and `barcode` columns) raised `KeyError`; it now says
  which layouts are understood.
* Supplementary alignments were labelled `Suplementary` in the alignment
  table and pie, and the categories came in whatever order the BAM showed
  them, with colours handed out by position, so a category's colour changed
  from run to run. Each category now has its own.
* An `N` in an MD tag was not counted as a mismatch.
* Run details such as a sample name holding `&` or `<` went into the HTML
  page as markup: the template's name ends in `.j2`, which the escaping rule
  did not match. Source file paths are escaped too, and the page is written
  as UTF-8 whatever the locale.
* A flat 2D density, as from a handful of identical reads, gave plotly a
  contour step of zero.
* A configuration naming something other than a plot, such as
  `summary_stats_dict`, failed after calling it.
* POD5 read ids repeated across files were dropped in silence, and POD5 or
  FAST5 files holding no reads gave an empty summary rather than an error.
* pandas 3.1 deprecates the `inplace` keyword. The four remaining
  `inplace=True` calls now assign the result instead.
* kaleido 1.5.1 often runs its script in the blank page it is leaving rather
  than the one it loads, failing with "kaleido_scopes is not defined" and
  aborting the PDF report. The `pdf` extra now excludes that version, and
  each page gets three attempts in case of other transient browser failures.
* The FASTQ and FAST5 converters forked their workers from a process already
  running threads (OpenBLAS's pool, pyarrow's jemalloc, or whatever pytest
  or a notebook had started), which can deadlock a child and which Python
  3.12 warns about. On Linux, workers now start from a fork server, as is
  the default from Python 3.14; macOS and Windows already used spawn. The server preloads pycoQC with those thread
  pools held off, so it stays single threaded and the workers share its
  imports: eight workers take about 100 MB rather than the 470 MB they
  would importing pandas each. Starting the server adds under a second to
  a conversion.
* A converter worker that died without reporting, killed for memory or
  failing as it started, left the parent waiting for ever. The parent now
  checks its workers' exit codes and stops with an error, which for a
  script missing its `if __name__ == "__main__":` guard says so.
* The buttons switching a plot between read sets ran down its left hand
  side, where the legend and the y axis title also were. They covered one
  or the other, and on the output and channel plots, with six buttons,
  both. Every plot's buttons now sit in a row above the plot area, under
  the title.
* The percentile labels of the distribution plots printed over one another
  wherever the percentiles were close together, which on a translocation
  speed plot left the middle three unreadable. Each label now sits a step
  lower than the one before.
* The wiki's example figures are re-rendered with these changes. New ones
  show the translocation speed, the end reasons, the active channels, the
  time in pore of a run without basecalls, and the read length per barcode.

**Changed**

* The original pycoQC 2.5.2 documentation site (`mkdocs.yml` and the pages,
  notebooks, example data and reports under `docs/`) is removed. Nothing
  built or linked to it, and it described the tool as it was before this
  fork. The wiki, kept in `docs/wiki/`, is the documentation. The three
  summaries the tests read moved to `tests/data/`.
* All four converters share one worker pipeline (lister, parsers, writer) in
  `pycoQC/_pipeline.py`, rather than each carrying its own copy of the
  queues, error handling, file listing, argument checks and cleanup. A
  converter now supplies only how a file becomes records and any last
  change to the table. Rows come out in the order the files are listed,
  whatever the number of threads, and what the user is told is logged by
  the calling process rather than a worker. The summary is written only
  once every worker has finished: a conversion that fails leaves none
  behind, and an earlier summary at the same path is kept. Read ids
  repeated across files and within one file are reported apart, since only
  the first suggests the directory holds a run twice.
* `Bam_to_seq_summary` and `Pod5_to_seq_summary` parse in processes rather
  than threads. Their per-read work holds Python's GIL, so the threads only
  took turns: on 54,000 BAM records, 8 threads took 14.1 s against 5.9 s
  for 3. They now take 2.7 s with 8. They print the same progress and
  totals as the FASTQ and FAST5 converters in place of their own lines.
  Workers on every platform start with numpy's and pyarrow's thread pools
  capped at one thread, which the converter modules used to set for the
  whole process on import. On Windows before Python 3.13 a conversion uses
  at most 31 threads, the most the calling process can wait on there, and
  says so when it lowers `--threads`.
* The figures in the run, basecall and alignment summary tables and in the
  JSON report are computed once per read set and shared, in place of a
  method per figure called again by each table and again by the report.
  The per barcode table groups the reads once rather than scanning the
  whole run for each barcode, twice: half the time with 3 barcodes, and
  the saving grows with the barcode count.
* `pycoQC_plot.py`, 3,200 lines holding every plot, is split by plot family
  into `pycoQC/plots/`: summary tables, distributions, over time, run
  health, barcodes and alignment, each a mixin, with the read sets and
  shared helpers in `base.py`. `pycoQC_plot` is composed from them and
  keeps every method, so nothing that uses it changes, and the reports are
  identical. An unused method, `_barcodes_count`, was dropped.
* The four converter commands are built by one function, which reads each
  `--fields` default from the converter class rather than keeping a second
  copy of the list. The copies had drifted apart once already. The BAM
  command's `--basecall_id` help described an option BAM ignores; it now
  says so.
* `--basecall_id` warns that it is ignored by the FASTQ converter too, as
  it already did for BAM and POD5; only a FAST5 has basecall groups.
* `Fastq_to_seq_summary` is about 1.8 times faster with six parsers. One
  process used to decompress every file and pass all of its text, 1.2 GB
  for 180,000 reads, to the parsers, and adding parsers barely helped. Each
  parser now reads its own files, and only a file well over an even share
  of the input is split into blocks of reads, so that one large file among
  small ones is still parsed by every parser. Parsers send records to the
  writer in batches rather than one read at a time.
* The command line prints a one line `ERROR:` message and exits with status
  1 for problems with the input, instead of a traceback. A converter's
  worker sends back the message alone for such problems, and its traceback
  only for a bug.
* A plot that fails unexpectedly is reported as a warning and the rest of
  the report is still written.
* The per barcode boxplots are drawn from statistics computed by pycoQC
  rather than from every sampled read, which went into the page four times:
  two boxplots were 12 MB of a 100,000 read report. Outlier points are no
  longer drawn. The channel map and the 2D densities are stored in single
  precision. The Dorado run with a BAM went from a 51 MB report in 105 s to
  12.7 MB in 17 s.
* Importing pycoQC no longer changes the process it is imported into:
  `FutureWarning` was silenced for every library, the converters configured
  the root logger at INFO, switching on other libraries' debug output, and
  set BLAS thread counts in the environment. The test suite now fails on a
  pandas `FutureWarning` raised from pycoQC code.
* The plotting module builds every three read set plot through the same
  helpers, rather than five copies of the code choosing a read set and
  building its buttons. Private helpers are no longer name mangled.
* The HTML and PDF reports share the loop that runs the configured plots.
* `check_arg` raises `pycoQCError` rather than a bare `Exception`. `head`
  reads gzipped files, as its documentation said.

### 2026-10-20 — pycoqc-ng 3.1.0

A report for runs that have not been basecalled, and a long list of
fixes found by reading the code against real Albacore, Guppy and Dorado
output rather than against the test fixtures.

**Added**

* A report for runs that have not been basecalled. POD5, and a FAST5 holding
  signal only, carry no sequence, so `pycoQC` used to stop with an error.
  It now reports on the run itself: reads over time, time in pore, activity
  per channel and reads per barcode, with a signal summary table in place of
  the basecall one. There is deliberately no estimated base count, since
  bases can only be guessed from the signal length.
* Pass and fail come from the `passes_filtering` column when there is no
  quality score to threshold on.

**Fixed**

* `--filter_calibration` did nothing on any Guppy 2.2+ run, which names the
  column `calibration_strand_genome` rather than `..._genome_template`.
  Recognising the name alone was not enough: the filter kept a list of the
  markers meaning "no control found", and Guppy's markers are not Albacore's,
  so the whole run would have been discarded as calibration strands.
* Barcodes written as bare numbers, as Deepbinner writes them, were typed as
  numbers and so split across chunks, became `6.0` after a merge, and could
  not be sorted against `unclassified`.
* The FASTQ and BAM converters counted elapsed time from the first read in
  the directory rather than from the start of each run, so two one hour runs
  recorded two days apart were reported as a single 50 hour run.
* Reads missing an optional column were discarded from the read count, the
  yield and the N50.
* Histogram bins: the shortest read fell outside the log bins, the axis ran a
  quarter of a decade past the longest read, the 1D density was drawn half a
  bin to the right, and the 2D density half a bin left and down.
* The yield on the over time plot disagreed with the summary table, because
  bases were estimated from a sample scaled by the ratio of read counts.
* The two mean coverage figures in one report measured different things.
* The barcode figure's buttons took two clicks and only updated the table;
  at 96 barcodes the legend, the pie and the table were unreadable or cut
  off, in the PDF as well.
* The converters stopped descending into a directory as soon as they found
  one matching file, so a `fastq_pass` holding `unclassified.fastq` beside
  one folder per barcode converted the unclassified reads and nothing else.
* `--fields` and `--include_path` were accepted and ignored by
  `Bam_to_seq_summary`, as was `--include_path` by `Pod5_to_seq_summary`;
  `--basecall_id` now says it has no meaning for those formats rather than
  being silently dropped. `Fast5_to_seq_summary`'s command line default
  left out four run level fields its Python default includes; both now come
  from one `DEFAULT_FIELDS` per converter. A field named twice in `--fields`
  is read and counted once; the BAM converter wrote it as two columns. In
  Python, `fields=None` or `fields=[]` gives the converter's default, where
  FASTQ wrote empty rows and FAST5 failed, and asking for `path` alone
  without `--include_path` gives the read ids. The `fields` argument of
  `Fastq_to_seq_summary` and `Fast5_to_seq_summary` now defaults to `None`,
  standing for the class's `DEFAULT_FIELDS`, which the notebook help for the
  class (`jhelp(Fastq_to_seq_summary)`) shows as the default; a subclass
  changes its default there. `fields` takes a single name as a string, or
  names in any ordered iterable: a list, a tuple, a numpy or pandas array, a
  pandas Index, `dict.keys()`, a generator. Blank names alone, as
  `--fields ""` gives, are the default; a blank among names is refused.
  Bytes, sets, dicts and their values, Series, files and numbers are
  refused with an error naming `fields`, as each gave wrong names, or
  columns in a different order on each run.
* FASTQ records were collected into a dict keyed by read id, so repeated
  ids merged, but only when they fell in the same 1,000 read block.
* A bam or fastq directory holding no reads wrote an empty summary or
  raised `KeyError`.
* `--config_file`'s help has always promised a lookup in `~/.pycoQC` and
  `~/.config/pycoQC/config`; nothing looked there.
* Pass and fail were split by dropping the passing reads' index labels,
  which takes every row carrying a repeated read id: a run basecalled
  twice reported no failing reads at all.
* Numeric barcodes came out as `6.0`, and sorting their names put barcode
  10 between 1 and 2.
* A summary with a read length but no quality score, or the reverse, was
  reported as if nothing had been basecalled instead of saying which
  column was missing.
* `Fast5_to_seq_summary` left `duration` and `passes_filtering` out of its
  defaults, so a raw FAST5 report had no duration plots and nothing
  marked as failing.
* A zero duration made every log bin edge NaN; a blank `passes_filtering`
  counted as a failure; a barcode too rare to be sampled showed a read
  count beside a blank N50.
* Pass and fail are read only from the part of the path below the
  directory given to the converter, so a run kept under a directory named
  for a failed experiment is no longer a run of failing reads.
* The FASTQ and FAST5 converters sent their workers SIGTERM on the error
  path and never waited for them, so a child that outlived the signal was
  left for Python to join at interpreter exit: the conversion raised the
  error it meant to, and the process then hung on the way out with nothing
  to say why. Cleanup now waits, and SIGKILLs whatever ignored SIGTERM, on
  every exit path rather than for three hand-listed exception types.

**Packaging**

* The vendored stylesheet directory was not declared as a package, which
  setuptools warned about on every build. The files were shipped anyway,
  by the package-data globs, so no release was missing them.

### 2026-10-19 — pycoqc-ng 3.0.0

First release on PyPI: `pip install pycoqc-ng`

Renamed the distribution to `pycoqc-ng`, continuing from pycoQC 2.5.2. The
import package is still `pycoQC`, so existing scripts keep working.

**Packaging**

* Metadata moved to `pyproject.toml`; the project is now pure Python and the
  conda recipe is `noarch`.
* Removed the Cython extension. The two callers used `import functions`, an
  absolute import of a module built at the repository root, so an installed
  copy crashed on use.
* Replaced `pkg_resources` with `importlib.resources`.
* Declared the missing `pod5` dependency and dropped `cython` at runtime.

**Fixed**

* `Fastq_to_seq_summary` reused one dict for every read in a file and queued
  that same object repeatedly; because `mp.Queue` pickles asynchronously, a
  200-read file produced 200 rows containing only 2 distinct reads.
* `max_fastq` processed one file fewer than requested.
* N50 sorted a pandas-owned array in place, which fails on read-only views.
* Writing elapsed timedeltas into a datetime column, rejected by pandas 3.
* `__over_time_data` dropped missing values from the value column only while
  reading times from the undropped frame, pairing values with the wrong times.
* `plotly` 7 compatibility (the removed `show_link` argument).
* `python -m pycoQC` did nothing, as `__main__.py` had no entry guard.
* `mkdir` raised `pycoMethError`, a name from another project, so any failure
  to create an output directory surfaced as a `NameError`.
* FASTQ chunk sizes were counted in lines. A fastq record is four lines, so a
  size that was not a multiple of four split a record and silently dropped
  it. Chunks are now counted in reads, and a block ending mid-record warns.

**Performance**

* Timestamps are parsed with `datetime.fromisoformat` instead of dateutil,
  160-220x faster per read.
* Mean q-scores are computed with a 256-bin count and one dot product.
* Together: 3.3x faster FASTQ parsing, with byte-identical output.
* Vectorised the over-time percentile binning and the coverage offsets.
* Flattened the nested multiprocessing in `Fastq_to_seq_summary`. Each worker
  process used to open a process pool of its own, so the process count was
  roughly the square of the thread count: 202 processes at `-t 16`, now 17.
  Removing the oversubscription also cut 48,000 reads over 16 files from
  4.16s to 3.33s, and a single 50,000-read file from 5.01s to 3.97s.

**Refactoring**

* Replaced `from pycoQC.common import *` and `from collections import *` with
  explicit imports throughout. Two modules had been using `logging` and `pd`
  without importing them, relying on them leaking through the star import,
  which also made the namespace fragile. With that done pyflakes resolves
  names for the first time and found several unused imports and variables,
  now removed.

**Project**

* Added a test suite (73 tests) and GitHub Actions CI across Python 3.9-3.13;
  replaced the dead Travis configuration.
* pyflakes runs as a CI lint job, so undefined names and dead imports fail
  the build.


### 18/11/2023 v-3.0.0.0

* Added support for Fastq input.
* If Fastq provided, %GC graphs are added
* Added support for BAM input if Doraro is used to basecall
* Added support for pod5 input
* Improved support for barcodes

### 11/10/2019 v-2.5.0.17

* Style update to comply with codacity
* Update summary table functions and split into 3 Run, Basecall and align

### 11/10/2019 v-2.5.0.14

* Add min_pass_len option to specify a minimal read length to be "pass"
* Update pycoQC usage doc and fix notebooks broken links

### 10/10/2019 v-2.5.0.12

* Start documenting changelog
* Change alignment score for identity frequency
* Reorganise documentation directory
* File names in HTML reports are now absolute paths
* Fix errors in Fast5_to_seq_summary
