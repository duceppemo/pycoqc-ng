# Changelog

### Unreleased

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
  left out four run level fields its Python default includes.
* FASTQ records were collected into a dict keyed by read id, so repeated
  ids merged, but only when they fell in the same 1,000 read block.
* A bam or fastq directory holding no reads wrote an empty summary or
  raised `KeyError`.
* `--config_file`'s help has always promised a lookup in `~/.pycoQC` and
  `~/.config/pycoQC/config`; nothing looked there.

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
