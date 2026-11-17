# -*- coding: utf-8 -*-
"""How the run went, rather than what it produced.

Why each read ended, which of a channel's pores it came from, and how fast
the strand passed through, carried from every input that records them to
three plots: end reasons over time, active channels over time, and the
translocation speed.
"""

import json
import logging

import numpy as np
import pandas as pd
import pytest

from conftest import REPO_ROOT, _write_basecalled_fast5
from pycoQC.Bam_to_seq_summary import Bam_to_seq_summary
from pycoQC.common import run_health_fields
from pycoQC.Fast5_to_seq_summary import Fast5_to_seq_summary
from pycoQC.Fastq_to_seq_summary import FastqParser
from pycoQC.Pod5_to_seq_summary import Pod5_to_seq_summary
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC.pycoQC_plot import LEVELS, pycoQC_plot
from pycoQC.pycoQC_report import pycoQC_report

REASONS = ["signal_positive"] * 8 + ["unblock_mux_change", "data_service_unblock_mux_change"]


def write_dorado_like(path, n=1000, seed=0, **overrides):
    """A summary with the run health columns Dorado writes."""
    rng = np.random.default_rng(seed)
    length = rng.integers(200, 20000, n)
    template = length / 400.0
    data = {
        "read_id": [f"r{i:05d}" for i in range(n)],
        "run_id": "run1",
        "channel": rng.integers(1, 513, n),
        "mux": rng.integers(1, 5, n),
        "start_time": np.sort(rng.uniform(0, 48 * 3600, n)),
        "duration": template + 0.25,
        "template_duration": template,
        "sequence_length_template": length,
        "mean_qscore_template": rng.uniform(5, 20, n).round(2),
        "end_reason": [REASONS[i % len(REASONS)] for i in range(n)],
    }
    data.update(overrides)
    pd.DataFrame(data).to_csv(path, sep="\t", index=False)
    return path


def plotter_for(path, **kwargs):
    return pycoQC_plot(pycoQC_parse(summary_file=str(path), quiet=True), quiet=True, **kwargs)


@pytest.fixture(scope="module")
def dorado(tmp_path_factory):
    return plotter_for(write_dorado_like(tmp_path_factory.mktemp("dorado") / "s.txt"))


class TestParsing:

    def test_the_columns_are_kept(self, dorado):
        for col in ("end_reason", "mux", "template_duration"):
            assert col in dorado.all_df, col
        assert dorado.all_df["mux"].dtype == np.uint8

    def test_end_reasons_are_lower_case(self, tmp_path):
        path = write_dorado_like(tmp_path / "s.txt", n=10,
                                 end_reason=["Signal_Positive "] * 10)
        assert set(plotter_for(path).all_df["end_reason"]) == {"signal_positive"}

    def test_fast5_start_mux_is_the_mux(self, tmp_path):
        path = tmp_path / "s.txt"
        pd.DataFrame({"read_id": ["a", "b"], "start_mux": [1, 3],
                      "sequence_length_template": [100, 200],
                      "mean_qscore_template": [10.0, 11.0]}).to_csv(path, sep="\t", index=False)
        assert list(plotter_for(path).all_df["mux"]) == [1, 3]


class TestSpeed:

    def test_speed_uses_the_template_time(self, dorado):
        """The whole duration counts the adapter, which drags short reads
        down; the template time does not."""
        assert dorado.speed_from_template
        np.testing.assert_allclose(dorado.all_df["speed"], 400, rtol=1e-4)

    def test_whole_duration_is_the_fallback(self, tmp_path):
        path = tmp_path / "s.txt"
        pd.DataFrame({"read_id": ["a", "b", "c"], "duration": [2.0, 0.0, 4.0],
                      "sequence_length_template": [800, 500, 1600],
                      "mean_qscore_template": [10.0, 11.0, 12.0]}).to_csv(path, sep="\t", index=False)
        plotter = plotter_for(path)
        assert not plotter.speed_from_template
        assert "adapter" in plotter._speed_label
        speed = plotter.all_df["speed"].to_numpy()
        assert speed[0] == 400 and np.isnan(speed[1]) and speed[2] == 400

    def test_no_duration_no_speed_plot(self, tmp_path):
        path = tmp_path / "s.txt"
        pd.DataFrame({"read_id": ["a", "b"], "sequence_length_template": [100, 200],
                      "mean_qscore_template": [10.0, 11.0]}).to_csv(path, sep="\t", index=False)
        plotter = plotter_for(path)
        assert not plotter.has_speed
        from pycoQC.common import pycoQCError
        for method in (plotter.read_speed_1D, plotter.read_speed_over_time):
            with pytest.raises(pycoQCError, match="speed"):
                method()

    def test_speed_plots_build(self, dorado):
        assert dorado.read_speed_1D().data[3].name == "Median"
        assert len(dorado.read_speed_over_time().data) == 5

    def test_speed_is_in_the_json(self, dorado, tmp_path):
        """Through the writer too: the float32 median broke json.dump."""
        out = tmp_path / "r.json"
        pycoQC_report(dorado.parser, dorado, quiet=True).json_report(str(out))
        stats = json.loads(out.read_text())["All Reads"]["basecall"]
        assert stats["median_speed"] == pytest.approx(400, rel=1e-4)
        assert len(stats["speed_percentiles"]) == 101


class TestEndReasons:

    def test_shares_add_up_and_counts_are_exact(self, dorado):
        fig = dorado.end_reason_over_time(time_bins=20)
        shares = np.array([t.y for t in fig.data], dtype=float)
        counts = np.array([t.customdata for t in fig.data], dtype=float)
        filled = counts.sum(axis=0) > 0
        np.testing.assert_allclose(shares[:, filled].sum(axis=0), 100)
        assert counts.sum() == len(dorado.all_df)
        expected = dorado.all_df["end_reason"].value_counts()
        for trace, row in zip(fig.data, counts):
            reason = trace.name.split(" (")[0]
            assert row.sum() == expected[reason]

    def test_most_common_first_with_its_share(self, dorado):
        names = [t.name for t in dorado.end_reason_over_time().data]
        assert names[0] == "signal_positive (80.0%)"
        assert set(names[1:]) == {"unblock_mux_change (10.0%)",
                                  "data_service_unblock_mux_change (10.0%)"}

    def test_a_reason_keeps_its_colour(self, dorado):
        fig = dorado.end_reason_over_time()
        colours = {t.name.split(" (")[0]: t.line.color for t in fig.data}
        assert colours["unblock_mux_change"] == pycoQC_plot.END_REASON_COLORS["unblock_mux_change"]

    def test_every_read_set_has_a_complete_update(self, tmp_path):
        plotter = plotter_for(write_dorado_like(tmp_path / "s.txt"), min_pass_qual=0)
        assert plotter.fail_df.empty
        fig = plotter.end_reason_over_time()
        buttons = fig.layout.updatemenus[0].buttons
        assert [b.label for b in buttons] == [plotter._level_label(level) for level in LEVELS]
        for button in buttons:
            assert len(button.args[0]["y"]) == len(fig.data)
        assert all(len(y) == 0 for y in buttons[2].args[0]["y"])

    def test_missing_end_reasons_count_as_unknown(self, tmp_path):
        reasons = ["signal_positive"] * 50 + [None] * 50
        plotter = plotter_for(write_dorado_like(tmp_path / "s.txt", n=100, end_reason=reasons))
        names = [t.name for t in plotter.end_reason_over_time().data]
        assert "unknown (50.0%)" in names

    def test_end_reasons_are_in_the_json(self, dorado):
        run = dorado.summary_stats_dict()["All Reads"]["run"]
        assert run["end_reasons"]["signal_positive"] == 800


class TestActiveChannels:

    def test_counts_distinct_channels_per_interval(self, dorado):
        fig = dorado.active_channels_over_time(interval_hours=4)
        df = dorado.all_df
        t = df["start_time"].to_numpy(dtype=np.float64) / 3600
        n = int(round((t.max() - t.min()) / 4))
        width = (t.max() - t.min()) / n
        idx = np.clip(np.floor((t - t.min()) / width).astype(int), 0, n - 1)
        assert len(fig.data[0].y) == n == 12
        np.testing.assert_allclose(fig.data[0].x, t.min() + width * (np.arange(n) + 0.5))
        expected = [df["channel"][idx == i].nunique() for i in range(n)]
        assert list(fig.data[0].y) == expected
        for m in (1, 2, 3, 4):
            trace = next(tr for tr in fig.data if tr.name == "Mux {}".format(m))
            mask = df["mux"].to_numpy() == m
            assert list(trace.y) == [df["channel"][mask & (idx == i)].nunique() for i in range(n)]

    def test_the_intervals_span_the_run_evenly(self, tmp_path):
        """With fixed 1.5 hour intervals, a 23.4 hour run ended in a
        sixteenth interval holding 0.9 hours of reads, drawn as a full one,
        and its channel count read as a collapse."""
        path = write_dorado_like(tmp_path / "s.txt", n=4000,
                                 start_time=np.linspace(0, 23.4 * 3600, 4000))
        x = plotter_for(path).active_channels_over_time().data[0].x
        width = 23.4 / 16
        assert len(x) == 16
        np.testing.assert_allclose(np.diff(x), width)
        assert x[-1] + width / 2 == pytest.approx(23.4)

    def test_a_short_run_still_gets_ten_intervals(self, tmp_path):
        path = write_dorado_like(tmp_path / "s.txt", n=300,
                                 start_time=np.linspace(0, 3600, 300))
        assert len(plotter_for(path).active_channels_over_time().data[0].y) == 10

    def test_read_sets_share_the_intervals(self, dorado):
        fig = dorado.active_channels_over_time()
        xs = [tuple(b.args[0]["x"][0]) for b in fig.layout.updatemenus[0].buttons]
        assert xs[0] == xs[1] == xs[2]

    def test_the_flow_cell_size_is_drawn(self, dorado):
        fig = dorado.active_channels_over_time()
        assert fig.data[-1].name == "Flow cell (512 channels)"
        assert list(fig.data[-1].y) == [512, 512]

    def test_without_a_mux_there_is_one_line_and_the_flow_cell(self, tmp_path):
        path = write_dorado_like(tmp_path / "s.txt", n=200)
        df = pd.read_csv(path, sep="\t").drop(columns="mux")
        df.to_csv(path, sep="\t", index=False)
        fig = plotter_for(path).active_channels_over_time()
        assert [t.name for t in fig.data] == ["All channels", "Flow cell (512 channels)"]


class TestConverters:

    def test_dorado_tags(self):
        fields = run_health_fields({"du": "10.0", "ns": 40000, "ts": 4000,
                                    "er": "signal_positive", "mx": 3})
        assert fields == {"duration": 10.0, "template_duration": 9.0,
                          "end_reason": "signal_positive", "mux": 3}
        assert run_health_fields({}) == {"duration": "", "template_duration": "",
                                         "end_reason": "", "mux": ""}

    def test_fastq_header_tags(self):
        record = ["@r1 ch:i:5 du:f:10.0 ns:i:40000 ts:i:4000 er:Z:unblock_mux_change mx:i:2",
                  "ACGT", "+", "IIII"]
        entry = FastqParser.single_fastq_entry_to_dict(record, "pass")["r1"]
        assert entry["template_duration"] == 9.0
        assert entry["end_reason"] == "unblock_mux_change"
        assert entry["mux"] == 2

    def test_bam_tags(self, tmp_path):
        pysam = pytest.importorskip("pysam")
        d = tmp_path / "bam"
        d.mkdir()
        header = {"HD": {"VN": "1.6", "SO": "unknown"}}
        with pysam.AlignmentFile(str(d / "r.bam"), "wb", header=header) as bam:
            for i in range(3):
                rec = pysam.AlignedSegment()
                rec.query_name = f"r{i}"
                rec.query_sequence = "ACGT" * 100
                rec.query_qualities = pysam.qualitystring_to_array("I" * 400)
                rec.flag = 4
                rec.set_tags([("du", 2.0, "f"), ("ns", 10000, "i"), ("ts", 1000, "i"),
                              ("er", "signal_positive", "Z"), ("mx", i + 1, "i"),
                              ("ch", 7, "i"), ("st", "2025-01-06T10:06:36+00:00", "Z")])
                bam.write(rec)
        out = tmp_path / "s.txt"
        Bam_to_seq_summary(bam_dir=str(d), seq_summary_fn=str(out), threads=3)
        df = pd.read_csv(out, sep="\t")
        assert list(df["template_duration"]) == [1.8] * 3
        assert list(df["mux"]) == [1, 2, 3]
        assert set(df["end_reason"]) == {"signal_positive"}

    def test_pod5(self, tmp_path, pod5_dir):
        out = tmp_path / "s.txt"
        Pod5_to_seq_summary(pod5_dir=str(pod5_dir), seq_summary_fn=str(out), threads=3)
        df = pd.read_csv(out, sep="\t")
        assert set(df["end_reason"]) == {"signal_positive"}
        assert set(df["mux"]) == {1}

    def test_fast5_start_mux(self, tmp_path):
        _write_basecalled_fast5(tmp_path / "b.fast5", n_reads=4)
        out = tmp_path / "s.txt"
        Fast5_to_seq_summary(fast5_dir=str(tmp_path), seq_summary_fn=str(out), threads=3)
        assert set(pd.read_csv(out, sep="\t")["start_mux"]) == {1}


class TestReport:

    def test_the_default_configs_carry_the_plots(self):
        templates = REPO_ROOT / "pycoQC" / "templates"
        basecalled = json.loads((templates / "pycoQC_config.json").read_text())
        signal = json.loads((templates / "pycoQC_config_signal.json").read_text())
        for name in ("read_speed_1D", "read_speed_over_time", "end_reason_over_time",
                     "active_channels_over_time"):
            assert name in basecalled, name
        assert "end_reason_over_time" in signal and "active_channels_over_time" in signal
        assert "read_speed_1D" not in signal

    def test_the_html_report_shows_them(self, dorado, tmp_path, caplog):
        out = tmp_path / "r.html"
        reporter = pycoQC_report(dorado.parser, dorado, quiet=True)
        with caplog.at_level(logging.WARNING):
            reporter.html_report(str(out))
        html = out.read_text(encoding="utf-8")
        for title in ("Translocation speed", "Read end reason over experiment time",
                      "Active channels over experiment time"):
            assert title in html, title
        assert "failed with" not in caplog.text
