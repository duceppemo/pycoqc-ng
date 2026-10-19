# Usage

## pycoQC

Builds the QC report from a basecaller summary file.

```bash
pycoQC -f sequencing_summary.txt -o pycoQC_report.html
```

With a Guppy barcoding file and a JSON summary as well:

```bash
pycoQC -f sequencing_summary.txt \
       -b barcoding_summary.txt \
       -o pycoQC_report.html \
       -j pycoQC_report.json
```

With alignment information from a BAM file:

```bash
pycoQC -f sequencing_summary.txt -a alignment.bam -o pycoQC_report.html
```

`-f`, `-b` and `-a` each accept several space-separated paths or a shell glob, so a run
split across files can be passed in one go:

```bash
pycoQC -f run1/sequencing_summary.txt run2/sequencing_summary.txt -o report.html
pycoQC -f "sequencing_summary_*.txt" -o report.html
```

At least one of `-o` or `-j` is required.

### Input and output

| Option | Description |
|---|---|
| `-f`, `--summary_file` | Summary from Albacore 1.0.0+, Guppy 2.1.3+ or Dorado 0.4.1+ (required) |
| `-b`, `--barcode_file` | Barcoding summary from Guppy 2.1.3+ or Deepbinner 0.2.0+ |
| `-a`, `--bam_file` | BAM of the same reads, preferably aligned with minimap2 |
| `-o`, `--html_outfile` | HTML report (required unless `-j` is given) |
| `-j`, `--json_outfile` | JSON summary (required unless `-o` is given) |

### Filtering

| Option | Default | Description |
|---|---|---|
| `--min_pass_qual` | 7 | Minimum mean quality for a read to count as *pass* |
| `--min_pass_len` | 0 | Minimum length for a read to count as *pass* |
| `--filter_calibration` | off | Drop reads flagged as calibration strand |
| `--filter_duplicated` | off | Drop duplicate `read_id`s, keeping the first (Guppy sometimes emits a read twice) |
| `--min_barcode_percent` | 0.1 | Barcodes below this share of reads are relabelled `unclassified` |

### Report

| Option | Default | Description |
|---|---|---|
| `--report_title` | `PycoQC report` | Title shown in the HTML |
| `--config_file` | — | JSON file choosing and tuning the plots, see [Configuration](Configuration) |
| `--template_file` | — | Custom Jinja2 template |
| `--skip_coverage_plot` | off | Skip the coverage plot, useful for transcriptome references with many sequences |
| `--sample` | 100000 | Reads sampled for plotting, deterministically. `0` uses everything |
| `-d`, `--default_config` | — | Print the default configuration and exit |

`--sample` only affects the plots; the summary statistics are computed on all reads.
Lowering it makes a large run render faster and produces a smaller HTML file.

### Verbosity

`-v` increases it, `-q` reduces it. They are mutually exclusive.

### Exit codes

`0` on success, `2` on a usage error such as a missing `--summary_file`. Worth knowing
if you drive the tool from Nextflow or Snakemake, which key off the exit status.

## Converters

If the basecaller did not write a summary, build one. Each converter walks a directory
tree recursively and writes a summary file that `pycoQC` can read. See
[Input formats](Input-formats) for which to use.

```bash
Fastq_to_seq_summary  -f fastq_dir/ -s sequencing_summary.txt
Bam_to_seq_summary    -f bam_dir/   -s sequencing_summary.txt
Pod5_to_seq_summary   -f pod5_dir/  -s sequencing_summary.txt
Fast5_to_seq_summary  -f fast5_dir/ -s sequencing_summary.txt
```

Shared options:

| Option | Default | Description |
|---|---|---|
| `-f`, `--*_dir` | — | Input directory, searched recursively (required) |
| `-s`, `--seq_summary_fn` | — | Output summary file (required) |
| `-t`, `--threads` | 4 | Total processes. One reads, one writes, the rest parse, so the minimum is 3 |
| `--max_*` | 0 | Stop after this many input files. `0` means no limit |
| `--fields` | see `--help` | Fields to extract |
| `--include_path` | off | Add a column with the absolute path of the source file |
| `--verbose_level` | 0 | `2` chatty, `1` normal, `0` quiet |

The output path must be in a directory that already exists.

## Barcode_split

Splits a summary into one file per barcode, writing
`sequencing_summary_<barcode>.txt` into the output directory, which is created if
needed. Barcode columns are dropped from the output, since the filename carries that.

```bash
Barcode_split -f sequencing_summary.txt -o split_dir/
```

| Option | Default | Description |
|---|---|---|
| `-f`, `--summary_file` | — | Summary file(s) (required) |
| `-b`, `--barcode_file` | — | Separate barcoding summary, if barcodes are not in the summary |
| `-o`, `--output_dir` | — | Where to write the per-barcode files |
| `-u`, `--output_unclassified` | off | Also write unclassified reads |
| `-p`, `--min_barcode_percent` | 0.1 | Skip barcodes below this share of reads |

## Python API

The API mirrors the command line:

```python
from pycoQC.pycoQC import pycoQC

pycoQC(summary_file="sequencing_summary.txt",
       html_outfile="pycoQC_report.html",
       json_outfile="pycoQC_report.json",
       min_pass_qual=7,
       sample=100000)
```

The converters are classes with the same arguments as their commands:

```python
from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary

Fastq_to_seq_summary(fastq_dir="fastq_dir/",
                     seq_summary_fn="sequencing_summary.txt",
                     threads=4)
```

In a notebook, `jhelp()` renders the signature and docstring as formatted Markdown
(requires the `notebook` extra):

```python
from pycoQC.common import jhelp
from pycoQC.pycoQC import pycoQC

jhelp(pycoQC)
```

To work with the parsed data directly:

```python
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC.pycoQC_plot import pycoQC_plot

parser = pycoQC_parse(summary_file="sequencing_summary.txt")
plotter = pycoQC_plot(parser, min_pass_qual=7)

fig = plotter.read_len_1D()       # a plotly figure
df = plotter.all_df               # the full dataframe
```
