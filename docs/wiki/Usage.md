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

At least one of `-o`, `-j` or `-p` is required.

A PDF alongside the HTML, for emailing or a lab notebook:

```bash
pycoQC -f sequencing_summary.txt -o report.html -p report.pdf
```

### Input and output

| Option | Description |
|---|---|
| `-f`, `--summary_file` | Summary from Albacore 1.0.0+, Guppy 2.1.3+ or Dorado 0.4.1+ (required) |
| `-b`, `--barcode_file` | Barcoding summary from Guppy 2.1.3+ or Deepbinner 0.2.0+ |
| `-a`, `--bam_file` | BAM of the same reads, preferably aligned with minimap2 |
| `-o`, `--html_outfile` | HTML report (at least one of `-o`, `-j` and `-p` is required) |
| `-j`, `--json_outfile` | JSON summary |
| `-p`, `--pdf_outfile` | Static PDF report, far smaller than the HTML. Needs the `pdf` extra |

### Filtering

| Option | Default | Description |
|---|---|---|
| `--min_pass_qual` | 7 | Minimum mean quality for a read to count as *pass* |
| `--min_pass_len` | 0 | Minimum length for a read to count as *pass* |
| `--filter_calibration` | off | Drop reads flagged as calibration strand |
| `--filter_duplicated` | off | Drop duplicate `read_id`s, keeping the first (Guppy sometimes emits a read twice) |
| `--min_barcode_percent` | 0.1 | Barcodes below this share of reads are relabelled `unclassified` |
| `--runid_list` | all | Report only these run ids, laid out in time in the order given. By default every run is kept and runs are ordered by decreasing throughput |

### Report

| Option | Default | Description |
|---|---|---|
| `--report_title` | `PycoQC report` | Title shown in the HTML |
| `--config_file` | — | JSON file choosing and tuning the plots, see [Configuration](Configuration) |
| `--template_file` | — | Custom Jinja2 template |
| `--skip_coverage_plot` | off | Leave the coverage plot out |
| `--sample` | 100000 | Reads sampled for plotting, deterministically. `0` uses everything |
| `-d`, `--default_config` | — | Print the default configuration and exit |

`--sample` only affects the plots; the summary statistics are computed on all reads.
Lowering it makes a large run render faster and produces a smaller HTML file.

### Verbosity

`-v` adds debug output, `-q` leaves only warnings and errors. They are mutually exclusive.
Progress goes to stderr, so it does not mix with anything you redirect from stdout.

### Exit codes

`0` on success, `2` on a usage error such as a missing `--summary_file`, and `1` when
the input cannot be used, with a one line `ERROR:` message saying why: a file that is
not there, a summary with no reads left after filtering, a barcode file in neither
known layout. The converters exit the same way. Worth knowing if you drive the tool
from Nextflow or Snakemake, which key off the exit status.

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

POD5, and a FAST5 holding signal only, carry no sequence, so a summary built from
either gives a shorter report: the run and the pores, without read lengths, quality
scores or yield.

Shared options:

| Option | Default | Description |
|---|---|---|
| `-f`, `--*_dir` | — | Input directory, searched recursively (required) |
| `-s`, `--seq_summary_fn` | — | Output summary file (required). A name ending in `.gz`, `.bz2`, `.xz` or `.zip` is written compressed |
| `-t`, `--threads` | 4 | Processes: one lists the input files, one writes the summary and the rest parse, so the minimum is 3 and each thread above it adds a parser. At most 31 on Windows before Python 3.13 |
| `--max_*` | 0 | Stop after this many input files. `0` means no limit |
| `--fields` | see `--help` | Columns to keep (FASTQ, BAM and FAST5; POD5 has no such option) |
| `--include_path` | off | Add a column with the absolute path of the source file |
| `--verbose_level` | 0 | `2` chatty, `1` normal, `0` quiet |

The output path must be in a directory that already exists. The summary is written only
once every input file has been read: a conversion that fails leaves no summary behind,
and an earlier one at the same path is kept.

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

if __name__ == "__main__":
    Fastq_to_seq_summary(fastq_dir="fastq_dir/",
                         seq_summary_fn="sequencing_summary.txt",
                         threads=4)
```

The FASTQ and FAST5 converters start worker processes from a fork server,
and each worker imports the calling script again as it starts. In a script,
keep the call under `if __name__ == "__main__":`, or every worker runs the
conversion again and the converter stops with an error saying so. Notebooks
need no guard.

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
