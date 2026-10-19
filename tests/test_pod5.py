# -*- coding: utf-8 -*-
"""Reading run metadata out of POD5.

POD5 stores raw signal and run information, not basecalls, so a summary
built from it describes the run rather than read quality.
"""

import pandas as pd
import pytest

from pycoQC.common import pycoQCError
from pycoQC.Pod5_to_seq_summary import Pod5_to_seq_summary


@pytest.fixture(scope="module")
def pod5_summary(tmp_path_factory, pod5_dir):
    out = tmp_path_factory.mktemp("pod5_out") / "sequencing_summary.txt"
    Pod5_to_seq_summary(pod5_dir=str(pod5_dir), seq_summary_fn=str(out),
                        threads=3)
    return out


class TestNoFabricatedData:

    def test_read_length_and_quality_are_absent(self, pod5_summary):
        """Regression: these were filled with np.random.normal draws, so every
        length and quality figure in a report built from POD5 was invented,
        down to negative read lengths."""
        df = pd.read_csv(pod5_summary, sep="\t")
        assert "read_len" not in df.columns
        assert "mean_qscore" not in df.columns
        assert "sequence_length_template" not in df.columns
        assert "mean_qscore_template" not in df.columns

    def test_output_is_reproducible(self, tmp_path, pod5_dir):
        """Random values made two runs over the same input disagree."""
        first = tmp_path / "a.txt"
        second = tmp_path / "b.txt"
        for out in (first, second):
            Pod5_to_seq_summary(pod5_dir=str(pod5_dir),
                                seq_summary_fn=str(out), threads=3)
        assert first.read_text() == second.read_text()

    def test_warns_that_it_cannot_feed_pycoqc(self, tmp_path, pod5_dir, caplog):
        with caplog.at_level("WARNING"):
            Pod5_to_seq_summary(pod5_dir=str(pod5_dir),
                                seq_summary_fn=str(tmp_path / "s.txt"),
                                threads=3)
        assert "raw signal" in caplog.text


class TestRunMetadata:

    def test_every_read_is_extracted(self, pod5_summary):
        df = pd.read_csv(pod5_summary, sep="\t")
        assert len(df) == 25
        assert df["read_id"].nunique() == 25

    def test_fields_come_from_the_run_info(self, pod5_summary):
        df = pd.read_csv(pod5_summary, sep="\t")
        assert set(df["run_id"]) == {"e2419f1214664b2b48a1f01eebf8dfa6f28fb57a"}
        assert set(df["sample_id"]) == {"test_sample"}
        assert set(df["flowcell_id"]) == {"FAW88668"}

    def test_channels_match_the_pore(self, pod5_summary):
        df = pd.read_csv(pod5_summary, sep="\t")
        expected = sorted((i * 7) % 512 + 1 for i in range(25))
        assert sorted(df["channel"]) == expected

    def test_start_time_is_seconds_not_samples(self, pod5_summary):
        """start_sample is divided by the run's sample rate of 5000 Hz, so
        reads 100000 samples apart are 20 seconds apart."""
        df = pd.read_csv(pod5_summary, sep="\t").sort_values("read_number")
        assert df["start_time"].iloc[0] == 0
        assert df["start_time"].iloc[1] == 20
        assert df["start_time"].max() == 24 * 100000 // 5000

    def test_duration_is_samples_over_rate(self, pod5_summary):
        df = pd.read_csv(pod5_summary, sep="\t")
        assert df["duration"].tolist() == pytest.approx([5000 / 5000] * 25)

    def test_read_numbers_are_preserved(self, pod5_summary):
        df = pd.read_csv(pod5_summary, sep="\t")
        assert sorted(df["read_number"]) == list(range(25))


@pytest.fixture(scope="module")
def barcoded(tmp_path_factory, pod5_barcoded_dir):
    out = tmp_path_factory.mktemp("pod5_bc") / "s.txt"
    Pod5_to_seq_summary(pod5_dir=str(pod5_barcoded_dir),
                        seq_summary_fn=str(out), threads=3)
    return pd.read_csv(out, sep="\t")


class TestFileNameHeuristics:

    def test_barcode_is_taken_from_the_file_name(self, barcoded):
        """MinKNOW puts the barcode in the file name; POD5 carries no
        barcode field of its own."""
        assert {"barcode06", "barcode11", "unclassified"} <= set(
            barcoded["barcode_arrangement"])

    def test_pass_and_fail_come_from_the_file_name(self, barcoded):
        assert set(barcoded["passes_filtering"]) == {True, False}
        failing = barcoded[barcoded["passes_filtering"] == False]  # noqa: E712
        assert len(failing) == 10

    def test_all_reads_across_files_are_kept(self, barcoded):
        assert len(barcoded) == 30
        assert barcoded["read_id"].nunique() == 30

    def test_name_without_barcode_is_marked_skip(self, pod5_summary):
        df = pd.read_csv(pod5_summary, sep="\t")
        assert set(df["barcode_arrangement"]) == {"skip"}


class TestFileDiscovery:

    def test_empty_directory_raises_a_clear_error(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(pycoQCError, match="[Nn]o valid pod5"):
            Pod5_to_seq_summary(pod5_dir=str(empty),
                                seq_summary_fn=str(tmp_path / "x.txt"),
                                threads=3)

    @pytest.mark.parametrize("max_pod5, expected_reads", [(1, 10), (2, 20), (3, 30)])
    def test_max_pod5_limits_files(self, tmp_path, pod5_barcoded_dir,
                                   max_pod5, expected_reads):
        out = tmp_path / f"s_{max_pod5}.txt"
        Pod5_to_seq_summary(pod5_dir=str(pod5_barcoded_dir),
                            seq_summary_fn=str(out), threads=3,
                            max_pod5=max_pod5)
        assert len(pd.read_csv(out, sep="\t")) == expected_reads

    def test_rejects_too_few_threads(self, tmp_path, pod5_dir):
        with pytest.raises(pycoQCError):
            Pod5_to_seq_summary(pod5_dir=str(pod5_dir),
                                seq_summary_fn=str(tmp_path / "x.txt"),
                                threads=2)
