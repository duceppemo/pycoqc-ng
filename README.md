<p align="center">
  <img src="https://raw.githubusercontent.com/duceppemo/pycoqc-ng/main/docs/images/pycoqc-ng-logo.svg" alt="pycoqc-ng" width="680">
</p>

<p align="center">
  <a href="https://github.com/duceppemo/pycoqc-ng/actions/workflows/ci.yml"><img src="https://github.com/duceppemo/pycoqc-ng/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="https://codecov.io/gh/duceppemo/pycoqc-ng"><img src="https://codecov.io/gh/duceppemo/pycoqc-ng/graph/badge.svg" alt="Coverage"></a>
  <a href="https://github.com/duceppemo/pycoqc-ng/releases/latest"><img src="https://img.shields.io/github/v/release/duceppemo/pycoqc-ng?label=release&cacheSeconds=3600" alt="Latest release"></a>
  <a href="https://anaconda.org/bioconda/pycoqc-ng"><img src="https://img.shields.io/conda/vn/bioconda/pycoqc-ng?label=bioconda" alt="Bioconda"></a>
  <a href="https://pypi.org/project/pycoqc-ng/"><img src="https://img.shields.io/pypi/v/pycoqc-ng?label=pypi" alt="PyPI"></a>
  <img src="https://img.shields.io/pypi/pyversions/pycoqc-ng" alt="Python versions">
  <a href="https://github.com/duceppemo/pycoqc-ng/blob/main/LICENSE"><img src="https://img.shields.io/github/license/duceppemo/pycoqc-ng" alt="License: GPLv3"></a>
  <a href="https://github.com/duceppemo/pycoqc-ng/wiki"><img src="https://img.shields.io/badge/docs-wiki-informational" alt="Documentation"></a>
  <a href="https://doi.org/10.5281/zenodo.23484147"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.23484147.svg" alt="DOI"></a>
</p>

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
* A lightweight static **PDF report** alongside the interactive HTML one
* No compiled extension, so it installs as a pure Python package

## Installation

```sh
pip install pycoqc-ng
```

Or the current development version:

```sh
pip install git+https://github.com/duceppemo/pycoqc-ng.git
```

Conda, once the recipe is merged:

```sh
conda create -n pycoqc-ng -c conda-forge -c bioconda pycoqc-ng
```

The bioconda recipe is [in review](https://github.com/bioconda/bioconda-recipes/pull/70469),
so the bioconda badge above reads "not found" until it lands.

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

Full documentation is in the [wiki](https://github.com/duceppemo/pycoqc-ng/wiki):
[Installation](https://github.com/duceppemo/pycoqc-ng/wiki/Installation) ·
[Usage](https://github.com/duceppemo/pycoqc-ng/wiki/Usage) ·
[Input formats](https://github.com/duceppemo/pycoqc-ng/wiki/Input-formats) ·
[Outputs](https://github.com/duceppemo/pycoqc-ng/wiki/Outputs) ·
[Configuration](https://github.com/duceppemo/pycoqc-ng/wiki/Configuration) ·
[FAQ](https://github.com/duceppemo/pycoqc-ng/wiki/FAQ)

The wiki is maintained in [`docs/wiki/`](docs/wiki) and published by a workflow, so
send documentation changes as pull requests against those files rather than editing
the wiki directly.

Example reports built from the bundled test data live in `docs/pycoQC/results/`.

## Citing

The method and the original implementation are published in JOSS. If you use
this tool, please cite the original authors:

> Leger A. and Leonardi T. (2019). pycoQC, interactive quality control for
> Oxford Nanopore Sequencing. *Journal of Open Source Software*, 4(34), 1236.
> https://doi.org/10.21105/joss.01236

[![JOSS](http://joss.theoj.org/papers/ea8e08dc950622bdd5d16a65649954aa/status.svg)](http://joss.theoj.org/papers/ea8e08dc950622bdd5d16a65649954aa)

To reference a specific version of this software, the Zenodo concept DOI
[10.5281/zenodo.23484147](https://doi.org/10.5281/zenodo.23484147) always resolves to
the latest release. See [CITATION.cff](CITATION.cff).

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
