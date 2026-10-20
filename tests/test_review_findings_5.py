# -*- coding: utf-8 -*-
"""Defects found reviewing the fourth pass's own changes."""

import os

import numpy as np
import pandas as pd
import pytest

from pycoQC.common import pycoQCError
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC.pycoQC_plot import pycoQC_plot


def base_columns(n):
    return {
        "read_id": [f"r{i:05d}" for i in range(n)],
        "run_id": "run1",
        "channel": [i % 512 + 1 for i in range(n)],
        "start_time": [float(i) for i in range(n)],
    }


class TestFailingReadsAreNotLost:
    """pass and fail were split by dropping the passing index labels. The
    index is the read id, duplicates are not filtered by default, and
    dropping a repeated label takes every row carrying it."""

    def test_duplicate_ids_do_not_empty_the_fail_set(self, tmp_path):
        n = 1000
        good = pd.DataFrame({**base_columns(n),
                             "sequence_length_template": 2000,
                             "mean_qscore_template": 12.0})
        # The same reads basecalled a second time, worse: two summaries of
        # one run in the same directory.
        repeat = good.head(50).copy()
        repeat["mean_qscore_template"] = 3.0
        path = tmp_path / "s.txt"
        pd.concat([good, repeat]).to_csv(path, sep="\t", index=False)

        plotter = pycoQC_plot(pycoQC_parse(summary_file=str(path), quiet=True),
                              quiet=True)
        assert len(plotter.all_df) == 1050
        assert len(plotter.pass_df) == 1000
        assert len(plotter.fail_df) == 50
        assert len(plotter.pass_df) + len(plotter.fail_df) == len(plotter.all_df)


class TestNumericBarcodeNames:

    @staticmethod
    def _summary(path, n=1200):
        pd.DataFrame({
            **base_columns(n),
            "sequence_length_template": [500 + i % 900 for i in range(n)],
            "mean_qscore_template": 10.0,
            # Deepbinner style numbers, with gaps so the column is a float.
            "barcode_arrangement": [(i % 12) + 1 if i % 20 else None
                                    for i in range(n)],
        }).to_csv(path, sep="\t", index=False)
        return path

    def test_a_whole_number_barcode_keeps_no_decimal_point(self, tmp_path):
        """Stringifying a float column carried ".0" into the report, the
        plot and the name of the file Barcode_split writes."""
        path = self._summary(tmp_path / "s.txt")
        parser = pycoQC_parse(summary_file=str(path), min_barcode_percent=0,
                              quiet=True)
        labels = set(parser.reads_df["barcode"])
        assert "6" in labels
        assert not any("." in label for label in labels)

    def test_barcodes_are_ordered_by_number(self, tmp_path):
        """Sorted as text, barcode 10 came between 1 and 2."""
        path = self._summary(tmp_path / "s.txt")
        plotter = pycoQC_plot(
            pycoQC_parse(summary_file=str(path), min_barcode_percent=0,
                         quiet=True), quiet=True)
        shown = list(plotter.basecall_summary_per_barcode()
                     ["data"][0]["cells"]["values"][0])
        assert shown == [str(i) for i in range(1, 13)] + ["unclassified"]


class TestHalfABasecallIsAnError:
    """has_basecalls needs both columns, so a summary with a length but no
    quality score quietly became a raw signal report, with the length
    ignored and the yield and N50 gone."""

    @pytest.mark.parametrize("present,missing", [
        ("sequence_length_template", "mean_qscore"),
        ("mean_qscore_template", "read_len"),
    ])
    def test_one_column_without_the_other_is_reported(self, tmp_path, present,
                                                      missing):
        n = 600
        pd.DataFrame({**base_columns(n), present: 900}).to_csv(
            tmp_path / "s.txt", sep="\t", index=False)
        with pytest.raises(pycoQCError, match=missing):
            pycoQC_parse(summary_file=str(tmp_path / "s.txt"), quiet=True)

    def test_neither_column_is_still_a_signal_run(self, tmp_path):
        n = 600
        pd.DataFrame({**base_columns(n), "duration": 3.0}).to_csv(
            tmp_path / "s.txt", sep="\t", index=False)
        parser = pycoQC_parse(summary_file=str(tmp_path / "s.txt"), quiet=True)
        assert pycoQC_plot(parser, quiet=True).has_basecalls is False


class TestZeroIsNotAMissingValue:

    def test_a_zero_duration_does_not_empty_the_histogram(self, tmp_path):
        """log10(0) is -inf, so every bin edge became NaN. Read lengths are
        filtered for this; durations are not, and a zero sample read is a
        real thing to find in a pod5."""
        n = 600
        duration = np.linspace(0.5, 9.0, n)
        duration[0] = 0.0
        pd.DataFrame({**base_columns(n), "duration": duration}).to_csv(
            tmp_path / "s.txt", sep="\t", index=False)

        plotter = pycoQC_plot(
            pycoQC_parse(summary_file=str(tmp_path / "s.txt"), quiet=True),
            quiet=True)
        histogram = plotter.summary_stats_dict()["All Reads"]["signal"][
            "duration_hist"]
        assert np.isfinite(np.asarray(histogram["x"], dtype=float)).all()
        assert np.nansum(np.asarray(histogram["y"], dtype=float)) > 0


class TestAMissingVerdictIsNotAFailure:

    def test_blank_passes_filtering_counts_as_unjudged(self, tmp_path):
        n = 600
        flags = ["True"] * (n - 10) + [""] * 10
        pd.DataFrame({**base_columns(n), "duration": 3.0,
                      "passes_filtering": flags}).to_csv(
            tmp_path / "s.txt", sep="\t", index=False)

        plotter = pycoQC_plot(
            pycoQC_parse(summary_file=str(tmp_path / "s.txt"), quiet=True),
            quiet=True)
        assert len(plotter.fail_df) == 0
        assert len(plotter.pass_df) == 600


class TestARareBarcodeKeepsItsStatistics:

    def test_a_barcode_missing_from_the_sample_is_measured_in_full(
            self, tmp_path):
        """Its row showed a real read count beside a blank N50 and median."""
        n = 20000
        barcodes = ["barcode01"] * (n - 3) + ["barcode02"] * 3
        pd.DataFrame({
            **base_columns(n),
            "sequence_length_template": [500 + i % 900 for i in range(n)],
            "mean_qscore_template": 10.0,
            "barcode_arrangement": barcodes,
        }).to_csv(tmp_path / "s.txt", sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(tmp_path / "s.txt"),
                              min_barcode_percent=0, quiet=True)
        # Small enough that the rare barcode is very unlikely to be sampled.
        plotter = pycoQC_plot(parser, sample=200, quiet=True)
        cells = plotter.basecall_summary_per_barcode()["data"][0]["cells"]["values"]
        rows = dict(zip(cells[0], zip(cells[1], cells[3])))
        reads, n50 = rows["barcode02"]
        assert reads == 3
        assert not np.isnan(n50)


class TestRawFast5GivesTheReportTheWikiPromises:
    """The FAST5 converter knew how to read duration and could tell pass
    from fail by the path, but neither was in its default fields, so a raw
    FAST5 summary produced a report with no duration plots and no failures."""

    @pytest.fixture
    def raw_fast5_summary(self, tmp_path):
        pytest.importorskip("h5py")
        import sys
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        from conftest import _write_basecalled_fast5
        from pycoQC.Fast5_to_seq_summary import Fast5_to_seq_summary

        (tmp_path / "fast5_pass").mkdir()
        (tmp_path / "fast5_fail").mkdir()
        _write_basecalled_fast5(tmp_path / "fast5_pass" / "a.fast5",
                                n_reads=20, basecalled=False, seed=1)
        _write_basecalled_fast5(tmp_path / "fast5_fail" / "b.fast5",
                                n_reads=10, basecalled=False, seed=2)
        out = tmp_path / "s.txt"
        Fast5_to_seq_summary(fast5_dir=str(tmp_path), seq_summary_fn=str(out),
                             threads=3)
        return out

    def test_duration_is_written_in_seconds(self, raw_fast5_summary):
        df = pd.read_csv(raw_fast5_summary, sep="\t")
        assert "duration" in df.columns
        # Recorded in samples; a read lasts seconds, not tens of thousands.
        assert df["duration"].max() < 1000

    def test_pass_and_fail_come_from_the_path(self, raw_fast5_summary):
        df = pd.read_csv(raw_fast5_summary, sep="\t")
        assert df["passes_filtering"].sum() == 20
        assert (~df["passes_filtering"]).sum() == 10

    def test_the_duration_plots_are_available(self, raw_fast5_summary):
        plotter = pycoQC_plot(
            pycoQC_parse(summary_file=str(raw_fast5_summary), quiet=True),
            quiet=True)
        assert plotter.has_duration
        assert plotter.read_duration_1D()["data"]
        assert plotter.read_duration_over_time()["data"]
        assert len(plotter.fail_df) == 10


class TestPassAndFailComeFromTheBasecallersPath:
    """Pass and fail are read from where the basecaller put a file, and the
    basecaller only chose the part of the path below the directory it was
    given. Reading the whole absolute path meant a run stored under a
    directory named for a failed experiment came out as failing reads."""

    def test_an_unrelated_parent_directory_is_not_a_verdict(self, tmp_path):
        from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary

        root = tmp_path / "failed_experiment_2024" / "run1"
        (root / "fastq_pass").mkdir(parents=True)
        with open(root / "fastq_pass" / "reads.fastq", "w") as fh:
            for i in range(20):
                fh.write("@r{} runid=r1 ch={} "
                         "start_time=2024-01-01T00:00:0{}Z\n{}\n+\n{}\n".format(
                             i, i + 1, i % 10, "ACGT" * 25, "I" * 100))

        out = tmp_path / "s.txt"
        Fastq_to_seq_summary(fastq_dir=str(root), seq_summary_fn=str(out),
                             threads=3)

        assert pd.read_csv(out, sep="\t")["passes_filtering"].all()

    def test_the_basecallers_own_fail_directory_still_counts(self, tmp_path):
        from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary

        root = tmp_path / "run1"
        (root / "fastq_fail").mkdir(parents=True)
        with open(root / "fastq_fail" / "reads.fastq", "w") as fh:
            for i in range(20):
                fh.write("@r{} runid=r1 ch={} "
                         "start_time=2024-01-01T00:00:0{}Z\n{}\n+\n{}\n".format(
                             i, i + 1, i % 10, "ACGT" * 25, "I" * 100))

        out = tmp_path / "s.txt"
        Fastq_to_seq_summary(fastq_dir=str(root), seq_summary_fn=str(out),
                             threads=3)

        assert not pd.read_csv(out, sep="\t")["passes_filtering"].any()
