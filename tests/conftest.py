# -*- coding: utf-8 -*-
"""Shared fixtures: synthetic ONT reads in the formats pycoQC consumes."""

import gzip
import math
import random
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
N_READS = 200


def reference_mean_qscore(qual_string):
    """Plain-Python mean q-score, used to check the vectorised implementation.

    Averages per-base error probabilities and converts back to a Phred score.
    """
    probs = [10 ** (-(ord(c) - 33) / 10) for c in qual_string]
    return -10 * math.log10(sum(probs) / len(probs))


def _random_reads(n_reads, seed=42):
    """Generate (read_id, sequence, quality string) tuples deterministically."""
    rng = random.Random(seed)
    reads = []
    for i in range(n_reads):
        length = rng.randint(200, 3000)
        seq = "".join(rng.choice("ACGT") for _ in range(length))
        qual = "".join(chr(rng.randint(2, 40) + 33) for _ in range(length))
        reads.append((f"read{i:05d}", seq, qual))
    return reads


@pytest.fixture(scope="session")
def reads():
    return _random_reads(N_READS)


@pytest.fixture(scope="session")
def fastq_dir(tmp_path_factory, reads):
    """A directory holding one gzipped FASTQ with Guppy-style headers."""
    d = tmp_path_factory.mktemp("fastq")
    with gzip.open(d / "reads.fastq.gz", "wt") as fh:
        for i, (read_id, seq, qual) in enumerate(reads):
            fh.write(
                f"@{read_id} runid=acfd15dcc73bb7bafd19e1fdbe72a535d7558bf1 "
                f"sampleid=testrun read={i} ch={(i * 7) % 512 + 1} "
                f"start_time=2022-06-27T23:{i % 60:02d}:28Z "
                f"barcode=barcode{i % 4 + 1:02d}\n"
                f"{seq}\n+\n{qual}\n")
    return d


@pytest.fixture(scope="session")
def multi_fastq_dir(tmp_path_factory, reads):
    """Three small FASTQ files, for exercising the max_fastq limit."""
    d = tmp_path_factory.mktemp("multi_fastq")
    per_file = 10
    for n in range(3):
        with gzip.open(d / f"reads_{n}.fastq.gz", "wt") as fh:
            for read_id, seq, qual in reads[n * per_file:(n + 1) * per_file]:
                fh.write(
                    f"@{read_id}_f{n} runid=run{n} ch=1 "
                    f"start_time=2022-06-27T23:00:28Z\n{seq}\n+\n{qual}\n")
    return d


@pytest.fixture(scope="session")
def bam_dir(tmp_path_factory, reads):
    """A directory holding one unaligned BAM carrying Dorado-style tags."""
    pysam = pytest.importorskip("pysam")
    d = tmp_path_factory.mktemp("bam")
    run_id = "e2419f1214664b2b48a1f01eebf8dfa6"
    read_group = f"{run_id}_dna_r10.4.1_e8.2_400bps_sup@v4.2.0"
    header = {
        "HD": {"VN": "1.6", "SO": "unknown"},
        "RG": [{
            "ID": read_group,
            "PU": "FAW88668",
            "PM": "promethion",
            "DS": f"basecall_model=dna_r10.4.1_e8.2_400bps_sup@v4.2.0 runid={run_id}",
        }],
    }
    with pysam.AlignmentFile(str(d / "reads.bam"), "wb", header=header) as bam:
        for i, (read_id, seq, qual) in enumerate(reads):
            rec = pysam.AlignedSegment()
            rec.query_name = read_id
            rec.query_sequence = seq
            rec.query_qualities = pysam.qualitystring_to_array(qual)
            rec.flag = 4
            rec.set_tags([
                ("BC", f"barcode{i % 4 + 1:02d}", "Z"),
                ("qs", 11, "i"),
                ("du", 31.6362, "f"),
                ("ch", (i * 7) % 512 + 1, "i"),
                ("st", f"2021-03-29T15:{i % 60:02d}:10.066+00:00", "Z"),
                ("rn", i, "i"),
                ("fn", f"FAW88668_pass_barcode06_1ec22c6e_{i}.pod5", "Z"),
                ("RG", read_group, "Z"),
            ])
            bam.write(rec)
    return d


@pytest.fixture(scope="session")
def guppy_summary():
    """A real (small) Guppy sequencing summary shipped with the docs."""
    fn = REPO_ROOT / "docs/pycoQC/data/Guppy-basecall-1D-DNA_sequencing_summary.txt.gz"
    if not fn.exists():
        pytest.skip(f"test data not available: {fn}")
    return fn


@pytest.fixture(scope="session")
def barcoded_summary(tmp_path_factory, guppy_summary):
    """A small sequencing summary carrying barcode assignments.

    Derived from the real summary shipped with the docs so the column schema
    stays authentic, with barcodes added and some reads left unclassified.
    """
    pd = pytest.importorskip("pandas")
    df = pd.read_csv(guppy_summary, sep="\t").head(400).copy()
    # Uneven barcode sizes, so the low-frequency cutoff has something to act on.
    assignments = (["barcode01"] * 200 + ["barcode02"] * 150
                   + ["barcode03"] * 45 + ["unclassified"] * 5)
    df["barcode_arrangement"] = assignments[:len(df)]
    out = tmp_path_factory.mktemp("barcoded") / "sequencing_summary.txt"
    df.to_csv(out, sep="\t", index=False)
    return out
