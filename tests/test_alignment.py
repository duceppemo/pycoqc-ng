# -*- coding: utf-8 -*-
"""Adding alignment information from a BAM file."""

import pandas as pd
import pytest

from pycoQC.common import pycoQCError
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC.pycoQC_plot import pycoQC_plot


def write_summary(path, n=600):
    pd.DataFrame({
        "read_id": [f"r{i:05d}" for i in range(n)],
        "run_id": "run1",
        "channel": [i % 512 + 1 for i in range(n)],
        "start_time": [float(i) for i in range(n)],
        "sequence_length_template": [200 + (i % 2000) for i in range(n)],
        "mean_qscore_template": [5.0 + (i % 10) for i in range(n)],
    }).to_csv(path, sep="\t", index=False)
    return path


@pytest.fixture(scope="module")
def aligned_parser(tmp_path_factory, aligned_bam):
    summary = write_summary(tmp_path_factory.mktemp("sum") / "s.txt")
    return pycoQC_parse(summary_file=str(summary), bam_file=str(aligned_bam),
                        quiet=True)


class TestAlignmentColumns:

    def test_alignment_columns_are_added(self, aligned_parser):
        df = aligned_parser.reads_df
        for column in ("align_len", "ref_id", "ref_start", "ref_end",
                       "identity_freq"):
            assert column in df, column

    def test_reference_lengths_are_recorded(self, aligned_parser):
        assert aligned_parser.ref_len_dict["chr1"] == 50000

    def test_identity_is_a_fraction(self, aligned_parser):
        identity = aligned_parser.reads_df["identity_freq"].dropna()
        assert len(identity) > 0
        assert identity.between(0, 1).all()

    def test_alignment_coordinates_are_within_the_reference(self, aligned_parser):
        df = aligned_parser.reads_df.dropna(subset=["ref_start", "ref_end"])
        assert (df["ref_start"] >= 0).all()
        assert (df["ref_end"] <= 50000).all()
        assert (df["ref_end"] > df["ref_start"]).all()

    def test_non_primary_alignments_are_counted_separately(self, aligned_parser):
        """Unmapped, secondary and supplementary records are tallied rather
        than contributing read statistics."""
        counts = aligned_parser.alignments_df.set_index("Alignments")["Counts"]
        assert {"Primary", "Unmapped", "Secondary"} <= set(counts.index)
        assert counts["Primary"] == 600
        assert counts["Unmapped"] == 1
        assert counts["Secondary"] == 1
        # Percentages are derived from the same tally.
        percents = aligned_parser.alignments_df["Percents"]
        assert percents.sum() == pytest.approx(100, abs=0.1)


@pytest.fixture(scope="module")
def plotter(aligned_parser):
    return pycoQC_plot(aligned_parser, quiet=True)


class TestAlignmentPlots:

    def test_alignment_is_detected(self, plotter):
        assert plotter.has_alignment

    @pytest.mark.parametrize("method", [
        "alignment_summary", "align_len_1D", "identity_freq_1D",
        "alignment_reads_status", "alignment_rate", "alignment_coverage",
    ])
    def test_alignment_plots_render(self, plotter, method):
        """These all raise pycoQCError without a BAM, so this is the first
        exercise they get."""
        assert getattr(plotter, method)() is not None

    def test_coverage_uses_the_reference_length(self, plotter):
        assert plotter.total_ref_len == 50000


class TestBamValidation:

    def test_unsorted_bam_is_rejected(self, tmp_path, aligned_bam):
        """pycoQC needs coordinate order to compute coverage."""
        pysam = pytest.importorskip("pysam")
        unsorted = tmp_path / "unsorted.bam"
        with pysam.AlignmentFile(str(aligned_bam), "rb") as src:
            header = src.header.to_dict()
            header["HD"]["SO"] = "unknown"
            with pysam.AlignmentFile(str(unsorted), "wb",
                                     header=header) as out:
                for rec in src:
                    out.write(rec)
        pysam.index(str(unsorted))

        summary = write_summary(tmp_path / "s.txt")
        with pytest.raises(pycoQCError, match="[Nn]ot sorted"):
            pycoQC_parse(summary_file=str(summary),
                         bam_file=str(unsorted), quiet=True)

    def test_unindexed_bam_is_rejected(self, tmp_path, aligned_bam):
        import shutil
        copied = tmp_path / "noindex.bam"
        shutil.copy(aligned_bam, copied)
        summary = write_summary(tmp_path / "s.txt")
        with pytest.raises(pycoQCError, match="[Nn]o index"):
            pycoQC_parse(summary_file=str(summary),
                         bam_file=str(copied), quiet=True)
