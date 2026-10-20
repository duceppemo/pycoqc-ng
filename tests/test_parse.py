# -*- coding: utf-8 -*-
"""Reading summary, barcode and alignment files together."""

import pandas as pd
import pytest

from pycoQC.common import pycoQCError
from pycoQC.pycoQC_parse import pycoQC_parse


def write_summary(path, n=600, run_id="run1", **overrides):
    df = pd.DataFrame({
        "read_id": [f"r{i:05d}" for i in range(n)],
        "run_id": run_id,
        "channel": [i % 512 + 1 for i in range(n)],
        "start_time": [float(i) for i in range(n)],
        "sequence_length_template": [200 + (i % 2000) for i in range(n)],
        "mean_qscore_template": [5.0 + (i % 10) for i in range(n)],
    })
    for key, value in overrides.items():
        df[key] = value
    df.to_csv(path, sep="\t", index=False)
    return path


class TestBarcodeFile:

    def test_guppy_barcode_file(self, tmp_path):
        """Guppy names the columns read_id and barcode_arrangement."""
        summary = write_summary(tmp_path / "s.txt", n=600)
        barcodes = tmp_path / "bc.txt"
        pd.DataFrame({
            "read_id": [f"r{i:05d}" for i in range(600)],
            "barcode_arrangement": [f"barcode{i % 3 + 1:02d}" for i in range(600)],
        }).to_csv(barcodes, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(summary),
                              barcode_file=str(barcodes), quiet=True)
        assert "barcode" in parser.reads_df
        assert set(parser.reads_df["barcode"]) >= {"barcode01", "barcode02"}

    def test_deepbinner_barcode_file(self, tmp_path):
        """Deepbinner names them read_ID and barcode_call."""
        summary = write_summary(tmp_path / "s.txt", n=600)
        barcodes = tmp_path / "bc.txt"
        pd.DataFrame({
            "read_ID": [f"r{i:05d}" for i in range(600)],
            "barcode_call": [f"barcode{i % 3 + 1:02d}" for i in range(600)],
        }).to_csv(barcodes, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(summary),
                              barcode_file=str(barcodes), quiet=True)
        assert "barcode" in parser.reads_df

    def test_none_becomes_unclassified(self, tmp_path):
        summary = write_summary(tmp_path / "s.txt", n=600)
        barcodes = tmp_path / "bc.txt"
        pd.DataFrame({
            "read_id": [f"r{i:05d}" for i in range(600)],
            "barcode_arrangement": ["none"] * 300 + ["barcode01"] * 300,
        }).to_csv(barcodes, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(summary),
                              barcode_file=str(barcodes), quiet=True)
        assert "unclassified" in set(parser.reads_df["barcode"])

    def test_reads_absent_from_the_barcode_file_are_unclassified(self, tmp_path):
        summary = write_summary(tmp_path / "s.txt", n=600)
        barcodes = tmp_path / "bc.txt"
        pd.DataFrame({
            "read_id": [f"r{i:05d}" for i in range(300)],
            "barcode_arrangement": ["barcode01"] * 300,
        }).to_csv(barcodes, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(summary),
                              barcode_file=str(barcodes), quiet=True)
        assert (parser.reads_df["barcode"] == "unclassified").sum() >= 300

    def test_low_frequency_barcodes_are_unset(self, tmp_path):
        summary = write_summary(tmp_path / "s.txt", n=600)
        barcodes = tmp_path / "bc.txt"
        # One barcode holds a single read, far below any sensible cutoff.
        calls = ["barcode01"] * 599 + ["barcode99"]
        pd.DataFrame({
            "read_id": [f"r{i:05d}" for i in range(600)],
            "barcode_arrangement": calls,
        }).to_csv(barcodes, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(summary),
                              barcode_file=str(barcodes),
                              min_barcode_percent=5, quiet=True)
        assert "barcode99" not in set(parser.reads_df["barcode"])


class TestFilters:

    def test_duplicated_reads_are_dropped(self, tmp_path):
        summary = tmp_path / "s.txt"
        write_summary(summary, n=600)
        df = pd.read_csv(summary, sep="\t")
        doubled = pd.concat([df, df.head(100)], ignore_index=True)
        doubled.to_csv(summary, sep="\t", index=False)

        kept = pycoQC_parse(summary_file=str(summary),
                            filter_duplicated=True, quiet=True)
        assert len(kept.reads_df) == 600

    def test_duplicates_are_kept_by_default(self, tmp_path):
        summary = tmp_path / "s.txt"
        write_summary(summary, n=600)
        df = pd.read_csv(summary, sep="\t")
        pd.concat([df, df.head(100)], ignore_index=True).to_csv(
            summary, sep="\t", index=False)
        parser = pycoQC_parse(summary_file=str(summary), quiet=True)
        assert len(parser.reads_df) == 700

    def test_calibration_reads_are_filtered(self, tmp_path):
        summary = write_summary(
            tmp_path / "s.txt", n=600,
            calibration_strand_genome_template=["no_match"] * 600)
        parser = pycoQC_parse(summary_file=str(summary),
                              filter_calibration=True, quiet=True)
        assert len(parser.reads_df) == 600

    def test_zero_length_reads_are_dropped(self, tmp_path):
        summary = tmp_path / "s.txt"
        write_summary(summary, n=600)
        df = pd.read_csv(summary, sep="\t")
        df.loc[:49, "sequence_length_template"] = 0
        df.to_csv(summary, sep="\t", index=False)
        parser = pycoQC_parse(summary_file=str(summary), quiet=True)
        assert len(parser.reads_df) == 550

    def test_selecting_an_unknown_runid_is_reported(self, tmp_path):
        summary = write_summary(tmp_path / "s.txt", n=600)
        with pytest.raises(pycoQCError, match="run ids"):
            pycoQC_parse(summary_file=str(summary),
                         runid_list=["nope"], quiet=True)

    def test_selecting_a_known_runid_keeps_it(self, tmp_path):
        summary = write_summary(tmp_path / "s.txt", n=600, run_id="runA")
        parser = pycoQC_parse(summary_file=str(summary),
                              runid_list=["runA"], quiet=True)
        assert len(parser.reads_df) == 600


class TestSeveralRuns:

    def test_runs_are_offset_in_time(self, tmp_path):
        """Two runs in one file must not overlap on the time axis."""
        first = write_summary(tmp_path / "a.txt", n=400, run_id="r1")
        second = write_summary(tmp_path / "b.txt", n=400, run_id="r2")
        parser = pycoQC_parse(summary_file=[str(first), str(second)],
                              quiet=True)
        df = parser.reads_df
        assert set(df["run_id"]) == {"r1", "r2"}
        spans = df.groupby("run_id")["start_time"].agg(["min", "max"])
        # The later run starts after the earlier one ends.
        ordered = spans.sort_values("min")
        assert ordered["min"].iloc[1] > ordered["max"].iloc[0]

    def test_missing_summary_file_is_reported(self, tmp_path):
        with pytest.raises(pycoQCError):
            pycoQC_parse(summary_file=str(tmp_path / "nope.txt"), quiet=True)

    def test_summary_without_a_read_id_is_reported(self, tmp_path):
        bad = tmp_path / "bad.txt"
        pd.DataFrame({"name": ["a", "b"], "something": [1, 2]}).to_csv(
            bad, sep="\t", index=False)
        with pytest.raises(pycoQCError, match="not found"):
            pycoQC_parse(summary_file=str(bad), quiet=True)
