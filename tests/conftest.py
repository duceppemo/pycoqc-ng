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


def _write_basecalled_fast5(path, n_reads=30, basecalled=True, seed=7):
    """A multi-read FAST5 in the layout Guppy wrote with --fast5_out.

    When basecalled is False the Analyses groups are omitted, which is what
    MinKNOW wrote on its own and what current Guppy and Dorado produce.
    """
    h5py = pytest.importorskip("h5py")
    import numpy as np
    rng = np.random.default_rng(seed)
    run_id = "acfd15dcc73bb7bafd19e1fdbe72a535d7558bf1"

    with h5py.File(path, "w") as fh:
        fh.attrs["file_type"] = np.bytes_("multi-read")
        fh.attrs["file_version"] = np.bytes_("2.0")
        for i in range(n_reads):
            read_id = f"{i:08x}-b178-468d-ab62-af921129c079"
            grp = fh.create_group(f"read_{read_id}")

            raw = grp.create_group("Raw")
            raw.attrs["read_id"] = np.bytes_(read_id)
            raw.attrs["start_time"] = np.uint64(i * 400000)
            raw.attrs["duration"] = np.uint32(rng.integers(5000, 50000))
            raw.attrs["read_number"] = np.int32(i)
            raw.attrs["start_mux"] = np.uint8(1)

            channel = grp.create_group("channel_id")
            channel.attrs["channel_number"] = np.bytes_(str((i * 7) % 512 + 1))
            channel.attrs["digitisation"] = np.float64(8192.0)
            channel.attrs["offset"] = np.float64(-243.0)
            channel.attrs["range"] = np.float64(1500.0)
            channel.attrs["sampling_rate"] = np.float64(4000.0)

            tracking = grp.create_group("tracking_id")
            tracking.attrs["run_id"] = np.bytes_(run_id)
            tracking.attrs["sample_id"] = np.bytes_("test_sample")
            tracking.attrs["device_id"] = np.bytes_("MN23126")
            tracking.attrs["flow_cell_id"] = np.bytes_("FAH77625")

            if basecalled:
                summary = grp.create_group(
                    "Analyses/Basecall_1D_000/Summary/basecall_1d_template")
                summary.attrs["mean_qscore"] = np.float64(
                    round(float(rng.uniform(6, 16)), 3))
                summary.attrs["sequence_length"] = np.int64(
                    int(rng.integers(300, 9000)))
                barcoding = grp.create_group(
                    "Analyses/Barcoding_000/Summary/barcoding")
                barcoding.attrs["barcode_arrangement"] = np.bytes_(
                    f"barcode{i % 4 + 1:02d}")
    return path


@pytest.fixture(scope="session")
def basecalled_fast5_dir(tmp_path_factory):
    """One FAST5 carrying basecalls, as Guppy --fast5_out produced."""
    d = tmp_path_factory.mktemp("fast5_basecalled")
    _write_basecalled_fast5(d / "batch0.fast5")
    return d


@pytest.fixture(scope="session")
def raw_fast5_dir(tmp_path_factory):
    """A FAST5 with raw signal only, with no basecalls written back."""
    d = tmp_path_factory.mktemp("fast5_raw")
    _write_basecalled_fast5(d / "batch0.fast5", basecalled=False)
    return d


def _write_pod5(path, n_reads=25, seed=0, id_offset=0):
    """A POD5 file written with the real pod5 writer, carrying MinKNOW-like
    run metadata. POD5 holds raw signal only, so there are no basecalls."""
    pod5 = pytest.importorskip("pod5")
    import datetime
    import uuid
    import numpy as np

    run_info = pod5.RunInfo(
        acquisition_id="e2419f1214664b2b48a1f01eebf8dfa6f28fb57a",
        acquisition_start_time=datetime.datetime(
            2026, 3, 29, 15, 18, 10, tzinfo=datetime.timezone.utc),
        adc_max=4095, adc_min=-4096,
        context_tags={"sample_frequency": "5000"},
        experiment_name="test_experiment",
        flow_cell_id="FAW88668",
        flow_cell_product_code="FLO-MIN114",
        protocol_name="sequencing/sequencing_MIN114_DNA_e8_2_400K",
        protocol_run_id=str(uuid.uuid4()),
        protocol_start_time=datetime.datetime(
            2026, 3, 29, 15, 0, 0, tzinfo=datetime.timezone.utc),
        sample_id="test_sample",
        sample_rate=5000,
        sequencing_kit="sqk-nbd114-24",
        sequencer_position="MN23126",
        sequencer_position_type="MinION",
        software="MinKNOW 23.04.5",
        system_name="minit",
        system_type="MinION Mk1C",
        tracking_id={"device_id": "MN23126"},
    )

    rng = np.random.default_rng(seed)
    with pod5.Writer(str(path)) as writer:
        for i in range(n_reads):
            writer.add_read(pod5.Read(
                read_id=uuid.UUID(int=id_offset + i),
                pore=pod5.Pore(channel=(i * 7) % 512 + 1, well=1,
                               pore_type="not_set"),
                calibration=pod5.Calibration(offset=-243.0, scale=0.1875),
                read_number=i,
                start_sample=i * 100000,
                median_before=220.0,
                end_reason=pod5.EndReason.from_reason_with_default_forced(
                    pod5.EndReasonEnum.SIGNAL_POSITIVE),
                run_info=run_info,
                signal=rng.integers(-200, 2000, size=5000, dtype=np.int16),
            ))
    return path


@pytest.fixture(scope="session")
def pod5_dir(tmp_path_factory):
    """A directory holding one POD5 file whose name carries no barcode."""
    d = tmp_path_factory.mktemp("pod5")
    _write_pod5(d / "reads.pod5")
    return d


@pytest.fixture(scope="session")
def pod5_barcoded_dir(tmp_path_factory):
    """POD5 files named the way MinKNOW writes them, with the pass or fail
    state and the barcode in the file name."""
    d = tmp_path_factory.mktemp("pod5_barcoded")
    _write_pod5(d / "FAW88668_pass_barcode06_1ec22c6e_0.pod5",
                n_reads=10, id_offset=0)
    _write_pod5(d / "FAW88668_fail_barcode11_1ec22c6e_1.pod5",
                n_reads=10, seed=1, id_offset=1000)
    _write_pod5(d / "FAW88668_pass_unclassified_1ec22c6e_2.pod5",
                n_reads=10, seed=2, id_offset=2000)
    return d
