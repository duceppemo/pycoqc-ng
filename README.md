# pycoqc-ng

**Interactive QC metrics and plots for Oxford Nanopore sequencing data**

[![CI](https://github.com/duceppemo/pycoqc-ng/actions/workflows/ci.yml/badge.svg)](https://github.com/duceppemo/pycoqc-ng/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/pycoqc-ng.svg)](https://pypi.org/project/pycoqc-ng/)
[![License](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
[![Python](https://img.shields.io/pypi/pyversions/pycoqc-ng.svg)](https://pypi.org/project/pycoqc-ng/)

`pycoqc-ng` is a maintained continuation of [pycoQC](https://github.com/a-slide/pycoQC)
by Adrien Leger and Tommaso Leonardi, which is no longer developed. It keeps the
original reports and adds support for the current Oxford Nanopore toolchain.

## What this fork adds

* **Dorado / unaligned BAM input** via `Bam_to_seq_summary`
* **POD5 input** via `Pod5_to_seq_summary`
* **FASTQ input** via `Fastq_to_seq_summary`, including **%GC metrics and plots**
* Per-barcode **length and %GC boxplots**, and a basecall summary broken down
  by barcode
* Runs on current dependencies: numpy 2.x, pandas 2.x and 3.x, plotly 5 to 7,
  pysam 0.22+, Python 3.9 to 3.13
* No compiled extension, so it installs as a pure Python package

## Installation

```sh
pip install pycoqc-ng
```

Or the current development version:

```sh
pip install git+https://github.com/duceppemo/pycoqc-ng.git
```

A bioconda recipe is planned; see
[issues](https://github.com/duceppemo/pycoqc-ng/issues) for progress.

To work on the code:

```sh
git clone https://github.com/duceppemo/pycoqc-ng.git
cd pycoqc-ng
pip install -e ".[dev]"
pytest
```

## Usage

Generate a report from a basecaller summary:

```sh
pycoQC -f sequencing_summary.txt -o pycoQC_report.html
```

With a Guppy barcoding file, plus a machine-readable summary:

```sh
pycoQC -f sequencing_summary.txt -b barcoding_summary.txt \
       -o pycoQC_report.html -j pycoQC_report.json
```

If you have no summary file, build one from the basecaller output first:

```sh
# Dorado unaligned BAM
Bam_to_seq_summary -f bam_dir/ -s sequencing_summary.txt

# POD5
Pod5_to_seq_summary -f pod5_dir/ -s sequencing_summary.txt

# FASTQ (adds %GC)
Fastq_to_seq_summary -f fastq_dir/ -s sequencing_summary.txt

# legacy FAST5
Fast5_to_seq_summary -f fast5_dir/ -s sequencing_summary.txt
```

Split a summary by barcode:

```sh
Barcode_split -f sequencing_summary.txt -o split_dir/
```

Every command takes `--help`. `pycoQC --default_config` prints the plot
configuration, which you can edit and pass back with `--config_file`.

The API mirrors the CLI:

```python
from pycoQC.pycoQC import pycoQC

pycoQC(summary_file="sequencing_summary.txt",
       html_outfile="pycoQC_report.html",
       json_outfile="pycoQC_report.json")
```

## Documentation

The original documentation covers the plots and options in detail, and still
applies to the shared functionality: https://a-slide.github.io/pycoQC

Example reports built from the bundled test data live in `docs/pycoQC/results/`.

## Citing

The method and the original implementation are published in JOSS. If you use
this tool, please cite the original authors:

> Leger A. and Leonardi T. (2019). pycoQC, interactive quality control for
> Oxford Nanopore Sequencing. *Journal of Open Source Software*, 4(34), 1236.
> https://doi.org/10.21105/joss.01236

[![JOSS](http://joss.theoj.org/papers/ea8e08dc950622bdd5d16a65649954aa/status.svg)](http://joss.theoj.org/papers/ea8e08dc950622bdd5d16a65649954aa)
[![DOI](https://zenodo.org/badge/94531811.svg)](https://zenodo.org/badge/latestdoi/94531811)

## Disclaimer

This is research software under active development. It is tested on Linux and
macOS.

## Licence

GPLv3, unchanged from the original project. See [LICENSE](LICENSE).

Copyright © 2019-2024 Adrien Leger & Tommaso Leonardi
Copyright © 2023-2026 Marc-Olivier Duceppe

## Authors

* Adrien Leger and Tommaso Leonardi — original pycoQC
* Marc-Olivier Duceppe — this fork
