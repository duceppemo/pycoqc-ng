# -*- coding: utf-8 -*-
"""Tests for the FASTQ and BAM to sequencing-summary converters."""

import pandas as pd
import pytest

from pycoQC.Bam_to_seq_summary import Bam_to_seq_summary
from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary
from conftest import N_READS, reference_mean_qscore


@pytest.fixture(scope="module")
def fastq_summary(tmp_path_factory, fastq_dir):
    out = tmp_path_factory.mktemp("fastq_out") / "sequencing_summary.txt"
    Fastq_to_seq_summary(fastq_dir=str(fastq_dir), seq_summary_fn=str(out), threads=3)
    return pd.read_csv(out, sep="\t")


@pytest.fixture(scope="module")
def bam_summary(tmp_path_factory, bam_dir):
    out = tmp_path_factory.mktemp("bam_out") / "sequencing_summary.txt"
    Bam_to_seq_summary(bam_dir=str(bam_dir), seq_summary_fn=str(out), threads=3)
    return pd.read_csv(out, sep="\t")


class TestFastqToSeqSummary:

    def test_every_read_appears_exactly_once(self, fastq_summary):
        """Regression: a dict shared across reads was queued repeatedly, so
        the summary held the right number of rows but only a couple of
        distinct reads."""
        assert len(fastq_summary) == N_READS
        assert fastq_summary["read_id"].nunique() == N_READS
        assert not fastq_summary["read_id"].duplicated().any()

    def test_per_read_values_are_correct(self, fastq_summary, reads):
        by_id = fastq_summary.set_index("read_id")
        for read_id, seq, qual in reads:
            row = by_id.loc[read_id]
            assert row["sequence_length_template"] == len(seq)
            assert row["mean_qscore_template"] == pytest.approx(
                round(reference_mean_qscore(qual), 2), abs=0.011)
            expected_gc = round(
                (seq.count("G") + seq.count("C")) / len(seq) * 100, 2)
            assert row["gc_percent"] == pytest.approx(expected_gc, abs=0.011)

    def test_expected_columns_are_present(self, fastq_summary):
        for column in ("read_id", "run_id", "channel", "start_time",
                       "gc_percent", "sequence_length_template",
                       "mean_qscore_template", "barcode_arrangement"):
            assert column in fastq_summary.columns

    def test_start_time_is_elapsed_seconds(self, fastq_summary):
        assert fastq_summary["start_time"].min() == 0
        assert (fastq_summary["start_time"] >= 0).all()

    def test_barcodes_are_preserved(self, fastq_summary):
        assert set(fastq_summary["barcode_arrangement"]) == {
            "barcode01", "barcode02", "barcode03", "barcode04"}

    @pytest.mark.parametrize("max_fastq, expected_files", [(1, 1), (2, 2), (3, 3)])
    def test_max_fastq_reads_exactly_that_many_files(
            self, tmp_path, multi_fastq_dir, max_fastq, expected_files):
        """Regression: the limit broke out of the loop before enqueueing, so
        one file fewer than requested was parsed."""
        out = tmp_path / f"summary_{max_fastq}.txt"
        Fastq_to_seq_summary(
            fastq_dir=str(multi_fastq_dir), seq_summary_fn=str(out),
            max_fastq=max_fastq, threads=3)
        df = pd.read_csv(out, sep="\t")
        assert df["run_id"].nunique() == expected_files

    def test_rejects_too_few_threads(self, tmp_path, fastq_dir):
        from pycoQC.common import pycoQCError
        with pytest.raises(pycoQCError):
            Fastq_to_seq_summary(
                fastq_dir=str(fastq_dir),
                seq_summary_fn=str(tmp_path / "x.txt"), threads=2)


class TestBamToSeqSummary:

    def test_every_read_appears_exactly_once(self, bam_summary):
        assert len(bam_summary) == N_READS
        assert bam_summary["read_id"].nunique() == N_READS

    def test_per_read_values_are_correct(self, bam_summary, reads):
        by_id = bam_summary.set_index("read_id")
        for read_id, seq, qual in reads:
            row = by_id.loc[read_id]
            assert row["sequence_length_template"] == len(seq)
            assert row["mean_qscore_template"] == pytest.approx(
                round(reference_mean_qscore(qual), 2), abs=0.011)

    def test_dorado_metadata_is_extracted(self, bam_summary):
        assert set(bam_summary["flowcell_id"]) == {"FAW88668"}
        assert set(bam_summary["device_id"]) == {"promethion"}
        assert bam_summary["dorado_model"].iloc[0].startswith("dna_r10.4.1")

    def test_start_time_is_elapsed_seconds(self, bam_summary):
        assert bam_summary["start_time"].min() == 0
