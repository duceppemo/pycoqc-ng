# Installation

`pycoqc-ng` is a pure Python package and needs Python 3.9 or newer. It is tested on
Python 3.9 through 3.13, on Linux and macOS.

## conda (recommended)

A conda environment is the easiest route, because it also brings in `pysam` and `pod5`
without needing a compiler:

```bash
conda create -n pycoqc-ng -c conda-forge -c bioconda pycoqc-ng
conda activate pycoqc-ng
```

> The bioconda recipe is still in review. Until it is merged, use pip.

## pip

```bash
pip install pycoqc-ng
```

To work inside a virtual environment:

```bash
python -m venv venv
source venv/bin/activate
pip install pycoqc-ng
```

## From source

```bash
git clone https://github.com/duceppemo/pycoqc-ng.git
cd pycoqc-ng
pip install .
```

For development, install the extras and run the tests:

```bash
pip install -e ".[dev,notebook]"
pytest
```

See [Development](Development) for more.

## Checking the installation

```bash
pycoQC --version
```

should print `pycoqc-ng v3.0.0`. The package installs six commands:

| Command | Purpose |
|---|---|
| `pycoQC` | build the QC report |
| `Fastq_to_seq_summary` | summary from FASTQ |
| `Bam_to_seq_summary` | summary from unaligned BAM (Dorado) |
| `Pod5_to_seq_summary` | summary from POD5 |
| `Fast5_to_seq_summary` | summary from legacy FAST5 |
| `Barcode_split` | split a summary into one file per barcode |

Each takes `--help`. The package can also be run as a module, `python -m pycoQC`.

## Dependencies

Installed automatically:

`numpy` ≥1.24, `scipy` ≥1.10, `pandas` ≥2.0, `plotly` ≥5.16, `jinja2` ≥3.1,
`h5py` ≥3.10, `tqdm` ≥4.66, `pysam` ≥0.22, `python-dateutil` ≥2.8, `pod5` ≥0.2.4

The `notebook` extra adds `ipython`, needed only for the `jhelp()` helper in notebooks.
The `pdf` extra adds `kaleido` and `pypdf`, needed only for `--pdf_outfile`; kaleido
drives a headless browser, which is why it is not installed by default:

```sh
pip install 'pycoqc-ng[pdf]'
```

There is no compiled extension. Earlier versions of pycoQC shipped a Cython module;
it was replaced by a vectorised NumPy implementation that is both faster and exact, so
no build toolchain is required.
