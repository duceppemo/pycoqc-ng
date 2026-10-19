# Input formats

`pycoQC` itself reads one thing: a **sequencing summary** file. Everything else is a
converter that produces one.

```
                        ┌─ Albacore / Guppy / Dorado ─► sequencing_summary.txt ─┐
                        │                                                       │
FASTQ  ─► Fastq_to_seq_summary ────────────────────────────────────────────────►├─► pycoQC
BAM    ─► Bam_to_seq_summary  ─────────────────────────────────────────────────►│
POD5   ─► Pod5_to_seq_summary ─────────────────────────────────────────────────►│
FAST5  ─► Fast5_to_seq_summary ────────────────────────────────────────────────►┘
```

## Sequencing summary

Written by the basecaller: Albacore 1.0.0+, Guppy 2.1.3+ or Dorado 0.4.1+. This is the
preferred input, because it is what the basecaller actually recorded rather than
something recomputed afterwards.

Gzipped files are read transparently, so `sequencing_summary.txt.gz` works.

Columns are detected rather than required, and the report adapts: plots whose columns
are missing are skipped with a note in the log instead of failing. The ones that matter
most are `read_id`, `run_id`, `channel`, `start_time`, `sequence_length_template` and
`mean_qscore_template`.

### Barcodes

Barcodes can arrive two ways:

* already in the summary, as a `barcode_arrangement` column
* in a separate barcoding summary from Guppy's `guppy_barcoder` or Deepbinner, passed
  with `-b`

Barcodes representing less than `--min_barcode_percent` of reads (0.1% by default) are
relabelled `unclassified`, which keeps the per-barcode plots readable when a run has a
long tail of spurious assignments.

### Alignment

Pass a BAM with `-a` to unlock the alignment plots: aligned length, identity frequency,
coverage and alignment rate. Align against your reference with minimap2 first. The BAM
must contain the same reads as the summary.

For a transcriptome, or any reference with very many sequences, add
`--skip_coverage_plot`; the coverage plot has one position per reference base and gets
unwieldy.

## FASTQ

```bash
Fastq_to_seq_summary -f fastq_dir/ -s sequencing_summary.txt -t 8
```

Reads `.fastq`, `.fq`, `.fastq.gz` and `.fq.gz`, recursively. Metadata is taken from the
Guppy-style header:

```
@20a271c3-8f4b-466c-92b8-50a469aaa3e4 runid=acfd15dc sampleid=run3 read=2743 ch=207 start_time=2022-06-27T23:14:28Z barcode=barcode17
```

Missing fields are tolerated and left blank, so plain FASTQ works, just with fewer
plots available.

FASTQ input additionally yields **%GC per read**, which the basecaller summary does not
carry, enabling the %GC plots.

Reads under a path containing `fail` are marked as failing, which matches the
`fastq_pass` and `fastq_fail` layout the basecallers write.

## BAM (Dorado)

```bash
Bam_to_seq_summary -f bam_dir/ -s sequencing_summary.txt -t 8
```

For **unaligned** BAM as produced by `dorado basecaller`. Metadata comes from the tags
Dorado writes: `st` (start time), `ch` (channel), `BC` (barcode), `qs` (quality),
`rn` (read number), `fn` (source file), plus the read group for the run id, flowcell,
device and basecalling model. Like FASTQ, this gives %GC.

Not to be confused with `-a/--bam_file` on `pycoQC`, which takes an *aligned* BAM to add
alignment plots to a report.

## POD5

```bash
Pod5_to_seq_summary -f pod5_dir/ -s sequencing_summary.txt -t 8
```

POD5 holds the raw signal and run metadata, not basecalls, so the summary it produces
describes the run rather than read quality. Use it when you want channel and timing
information for data you have not basecalled.

## FAST5

```bash
Fast5_to_seq_summary -f fast5_dir/ -s sequencing_summary.txt -t 8
```

The legacy format, for data basecalled with Albacore or older Guppy where basecalls
were written back into the FAST5. Slower than the alternatives, because every file is
opened and walked.

## Which should I use?

| You have | Use |
|---|---|
| `sequencing_summary.txt` | it directly, no conversion |
| Dorado output | `Bam_to_seq_summary` |
| FASTQ, and want %GC | `Fastq_to_seq_summary` |
| POD5 and no basecalls yet | `Pod5_to_seq_summary` |
| Old FAST5 with basecalls | `Fast5_to_seq_summary` |

Prefer the basecaller's own summary where you have it. The converters recompute what
they can from the reads, which is close but not identical, and they are slower.
