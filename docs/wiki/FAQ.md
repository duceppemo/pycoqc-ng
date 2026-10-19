# FAQ

## How does this relate to pycoQC?

`pycoqc-ng` continues [pycoQC](https://github.com/a-slide/pycoQC), whose last release
was 2.5.2 in 2020. The method, the report and most of the code are the original authors'
work. This fork keeps it running on current dependencies and adds Dorado, POD5 and
FASTQ input.

The name had to change because `pycoQC` is taken on both PyPI and bioconda. The
**import package is still `pycoQC`**, so existing scripts keep working:

```python
from pycoQC.pycoQC import pycoQC   # unchanged
```

The commands kept their names too, so `pycoQC -f ... -o ...` is unchanged.

## Can I install it alongside the original pycoqc?

Not in the same environment. Both provide a `pycoQC` command and a `pycoQC` import
package, so they would overwrite each other. Use separate conda environments.

## `ERROR: --summary_file is a required argument`

`pycoQC` needs a summary file. If the basecaller did not write one, build one from your
reads first, see [Input formats](Input-formats).

## `No valid fastq files found in indicated folder`

The converter searches recursively for `.fastq`, `.fq`, `.fastq.gz` and `.fq.gz`. Check
you pointed at the directory rather than a file, and that the extensions match. The
converters take a **directory**, not a filename.

## `Cannot write the indicated seq_summary_fn`

The directory for the output file must already exist. The converters do not create it,
unlike `Barcode_split`, which does create its output directory.

## `At least 3 threads required`

The converters run one reader, one writer and at least one parser, so `-t` must be 3 or
more. Higher values add parsers.

## Why is my report so large?

Every plot and the plotly runtime are embedded so the file stands alone. Tens of MB is
normal for a barcoded run. To shrink it, lower `--sample`, or trim plots from the
configuration, see [Configuration](Configuration).

## The coverage plot is enormous or very slow

It has one position per reference base, so a transcriptome or a draft assembly with many
contigs is punishing. Pass `--skip_coverage_plot`.

## Some plots are missing

Plots are skipped when their data is not present, with a note in the log. Run with `-v`
to see which and why.

* alignment plots need `-a` with a BAM
* barcode plots need barcodes, in the summary or via `-b`
* %GC plots need a summary built from FASTQ or BAM; the basecaller's own summary carries
  no per-read %GC

## All my reads are unclassified

Barcodes below `--min_barcode_percent` (0.1% by default) are relabelled `unclassified`.
On a small or very uneven run, lower it:

```bash
pycoQC -f sequencing_summary.txt -o report.html --min_barcode_percent 0
```

## Pass and fail counts do not match the basecaller's

`pycoQC` applies its own threshold, `--min_pass_qual` 7 by default, rather than reading
the basecaller's flag. Set it to match your basecaller's setting if you want the numbers
to agree.

## Can I use it on a run that is still going?

Yes. The summary file grows as the run proceeds and a report built from it reflects the
reads written so far. Expect the final figures to differ.

## Does it work with RNA, or 1D²?

Yes for direct RNA, and the plots are the same. 1D² data from the older chemistries is
read, since the summary columns are detected rather than assumed.

## What Python versions are supported?

3.9 through 3.13, tested on Linux and macOS in CI.

## It used to need a compiler. Does it still?

No. The Cython extension is gone, replaced by a vectorised NumPy implementation that is
both faster and exact. The package is pure Python and the conda build is `noarch`.

## How should I cite it?

Cite the original paper for the method, and the Zenodo DOI if you need to reference a
specific software version. See [Home](Home#citing).
