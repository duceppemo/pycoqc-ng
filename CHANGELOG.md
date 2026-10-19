# Changelog

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
