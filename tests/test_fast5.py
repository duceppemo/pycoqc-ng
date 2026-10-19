# -*- coding: utf-8 -*-
"""Reading basecalls that Guppy wrote back into FAST5.

Guppy's --fast5_out stored the basecall summary inside the FAST5 under
/Analyses. Oxford Nanopore have since removed FAST5 writing, and Dorado does
not produce FAST5 at all, so this path is for existing data.
"""

import json

import pandas as pd
import pytest

from pycoQC.common import pycoQCError
from pycoQC.Fast5_to_seq_summary import Fast5_to_seq_summary


@pytest.fixture(scope="module")
def fast5_summary(tmp_path_factory, basecalled_fast5_dir):
    out = tmp_path_factory.mktemp("fast5_out") / "sequencing_summary.txt"
    Fast5_to_seq_summary(fast5_dir=str(basecalled_fast5_dir),
                         seq_summary_fn=str(out), threads=3)
    return out


class TestBasecalledFast5:

    def test_every_read_is_extracted(self, fast5_summary):
        df = pd.read_csv(fast5_summary, sep="\t")
        assert len(df) == 30
        assert df["read_id"].nunique() == 30

    def test_basecall_summary_is_used(self, fast5_summary):
        """The length and quality come from the Analyses group Guppy wrote."""
        df = pd.read_csv(fast5_summary, sep="\t")
        assert df["sequence_length_template"].between(300, 9000).all()
        assert df["mean_qscore_template"].between(6, 16).all()

    def test_string_fields_are_decoded(self, fast5_summary):
        """Regression: HDF5 byte strings were never decoded, so values
        reached the summary as b'barcode01' and a channel of b'78' that could
        not be cast to an integer."""
        df = pd.read_csv(fast5_summary, sep="\t", dtype=str)
        for column in ("read_id", "run_id", "channel", "barcode_arrangement"):
            assert not df[column].str.startswith("b'").any(), column
        assert df["channel"].astype(int).between(1, 512).all()
        assert set(df["barcode_arrangement"]) == {
            "barcode01", "barcode02", "barcode03", "barcode04"}

    def test_start_time_is_scaled_by_sampling_rate(self, fast5_summary):
        """start_time is raw samples; it must be divided by the rate."""
        df = pd.read_csv(fast5_summary, sep="\t")
        # 400000 samples apart at 4000 Hz is 100 seconds per read.
        assert df["start_time"].max() == pytest.approx(29 * 400000 / 4000)

    def test_summary_feeds_pycoqc(self, fast5_summary, tmp_path):
        from pycoQC.pycoQC import pycoQC
        html = tmp_path / "r.html"
        js = tmp_path / "r.json"
        pycoQC(summary_file=str(fast5_summary), html_outfile=str(html),
               json_outfile=str(js), quiet=True)
        assert html.exists()
        assert json.loads(js.read_text())["All Reads"]["basecall"]["reads_number"] == 30


class TestRawFast5:

    def test_raw_only_fast5_yields_no_basecall_columns(self, tmp_path, raw_fast5_dir):
        """A FAST5 without the Analyses group has no length or quality to
        report, and should not invent any."""
        out = tmp_path / "summary.txt"
        Fast5_to_seq_summary(fast5_dir=str(raw_fast5_dir),
                             seq_summary_fn=str(out), threads=3)
        df = pd.read_csv(out, sep="\t")
        for column in ("sequence_length_template", "mean_qscore_template"):
            assert column not in df or df[column].isna().all()


class TestFileDiscovery:

    def test_single_file_directory(self, tmp_path, basecalled_fast5_dir):
        """Regression: enumerate is zero based, so a directory holding
        exactly one file reported that no files were found."""
        out = tmp_path / "one.txt"
        Fast5_to_seq_summary(fast5_dir=str(basecalled_fast5_dir),
                             seq_summary_fn=str(out), threads=3)
        assert len(pd.read_csv(out, sep="\t")) == 30

    def test_empty_directory_raises_a_clear_error(self, tmp_path):
        """Regression: with no files the loop variable was never bound and
        an UnboundLocalError escaped instead of the intended message."""
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(pycoQCError, match="[Nn]o valid fast5"):
            Fast5_to_seq_summary(fast5_dir=str(empty),
                                 seq_summary_fn=str(tmp_path / "x.txt"),
                                 threads=3)

    def test_rejects_too_few_threads(self, tmp_path, basecalled_fast5_dir):
        with pytest.raises(pycoQCError):
            Fast5_to_seq_summary(fast5_dir=str(basecalled_fast5_dir),
                                 seq_summary_fn=str(tmp_path / "x.txt"),
                                 threads=2)
