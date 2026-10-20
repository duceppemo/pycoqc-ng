# -*- coding: utf-8 -*-
"""Reports for runs that have not been basecalled.

POD5, and a FAST5 holding signal only, record which pore produced which
read and when, but no sequence. There is no read length and no quality
score until a basecaller has run, and no base count at all: an estimate
from the signal length would be a guess dressed up as a measurement. What
the run itself did is still worth reporting.
"""

import json

import numpy as np
import pandas as pd
import pytest

from pycoQC.common import pycoQCError
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC.pycoQC_plot import pycoQC_plot

BASECALL_ONLY = ["basecall_summary", "basecall_summary_per_barcode",
                 "read_len_1D", "read_qual_1D", "read_len_read_qual_2D",
                 "read_len_over_time", "read_qual_over_time",
                 "length_boxplot_per_barcode"]

SIGNAL_PLOTS = ["run_summary", "signal_summary", "read_duration_1D",
                "read_duration_over_time", "output_over_time",
                "channels_activity", "barcode_counts"]


def write_signal_summary(path, n=3000, duration=True, flags=True,
                         barcodes=True, seed=1):
    rng = np.random.default_rng(seed)
    df = pd.DataFrame({
        "read_id": [f"r{i:06d}" for i in range(n)],
        "run_id": "sig1",
        "channel": rng.integers(1, 513, n),
        "start_time": np.sort(rng.uniform(0, 36 * 3600, n)).astype(int),
        "sample_id": "sample1",
        "flowcell_id": "FAW00001",
        "read_number": np.arange(n),
    })
    if duration:
        df["duration"] = rng.lognormal(1.2, 0.6, n).round(3)
    if barcodes:
        df["barcode_arrangement"] = [f"barcode{i % 12 + 1:02d}"
                                     for i in range(n)]
    if flags:
        df["passes_filtering"] = rng.random(n) > 0.1
    df.to_csv(path, sep="\t", index=False)
    return path


@pytest.fixture
def signal_plotter(tmp_path):
    path = write_signal_summary(tmp_path / "pod5_summary.txt")
    parser = pycoQC_parse(summary_file=str(path), quiet=True)
    return pycoQC_plot(parser, quiet=True)


class TestParsing:

    def test_a_summary_with_no_sequence_is_accepted(self, tmp_path):
        path = write_signal_summary(tmp_path / "s.txt", n=3000)
        parser = pycoQC_parse(summary_file=str(path), quiet=True)
        assert len(parser.reads_df) == 3000
        assert "read_len" not in parser.reads_df
        assert "mean_qscore" not in parser.reads_df

    def test_the_plotter_knows_there_are_no_basecalls(self, signal_plotter):
        assert signal_plotter.has_basecalls is False
        assert signal_plotter.has_duration is True

    def test_a_basecalled_run_is_unaffected(self, tmp_path):
        guppy = ("docs/pycoQC/data/"
                 "Guppy-2.2.4-basecall-1D-DNA_sequencing_summary+barcode.txt.gz")
        parser = pycoQC_parse(summary_file=guppy, quiet=True)
        assert pycoQC_plot(parser, quiet=True).has_basecalls is True


class TestPassAndFail:
    """Without a quality score there is nothing to threshold on, so the
    basecaller's own verdict is used: MinKNOW sorts its files into pass and
    fail folders and the converters carry that through."""

    def test_the_recorded_flag_decides(self, signal_plotter):
        total = len(signal_plotter.all_df)
        assert 0 < len(signal_plotter.pass_df) < total
        assert (len(signal_plotter.pass_df) + len(signal_plotter.fail_df)
                == total)

    def test_without_a_flag_nothing_has_failed(self, tmp_path):
        path = write_signal_summary(tmp_path / "s.txt", flags=False)
        plotter = pycoQC_plot(pycoQC_parse(summary_file=str(path), quiet=True),
                              quiet=True)
        assert len(plotter.pass_df) == len(plotter.all_df)
        assert len(plotter.fail_df) == 0


class TestWhichPlotsAppear:

    @pytest.mark.parametrize("method", SIGNAL_PLOTS)
    def test_the_signal_plots_build(self, signal_plotter, method):
        assert getattr(signal_plotter, method)()["data"]

    @pytest.mark.parametrize("method", BASECALL_ONLY)
    def test_the_basecall_plots_are_skipped_not_crashed(self, signal_plotter,
                                                        method):
        """The report skips a plot that raises pycoQCError, the same way it
        skips the alignment plots when there is no BAM. A KeyError would
        take the whole report down."""
        with pytest.raises(pycoQCError, match="basecall"):
            getattr(signal_plotter, method)()

    def test_duration_plots_need_a_duration_column(self, tmp_path):
        path = write_signal_summary(tmp_path / "s.txt", duration=False)
        plotter = pycoQC_plot(pycoQC_parse(summary_file=str(path), quiet=True),
                              quiet=True)
        for method in ("read_duration_1D", "read_duration_over_time"):
            with pytest.raises(pycoQCError, match="duration"):
                getattr(plotter, method)()


class TestWhatIsReported:

    def test_the_signal_table_counts_reads_and_pore_time(self, signal_plotter):
        cells = signal_plotter.signal_summary()["data"][0]["cells"]["values"]
        statuses, reads = list(cells[0]), list(cells[1])
        assert statuses == ["All Reads", "Pass Reads", "Fail Reads"]
        assert reads[0] == len(signal_plotter.all_df)
        assert reads[1] + reads[2] == reads[0]
        # Pore time is the summed duration, in hours.
        expected = signal_plotter.all_df["duration"].sum() / 3600
        assert cells[3][0] == pytest.approx(expected)

    def test_reads_over_time_counts_every_read(self, signal_plotter):
        _, data, _ = signal_plotter._pycoQC_plot__output_over_time_data(
            "all", "reads")
        assert np.asarray(data["y"][0])[-1] == len(signal_plotter.all_df)

    def test_there_is_no_base_count_anywhere(self, signal_plotter):
        stats = signal_plotter.summary_stats_dict()["All Reads"]
        assert "basecall" not in stats
        assert "bases_number" not in stats.get("signal", {})

    def test_the_stats_describe_the_signal(self, signal_plotter):
        signal = signal_plotter.summary_stats_dict()["All Reads"]["signal"]
        assert signal["reads_number"] == len(signal_plotter.all_df)
        assert signal["median_duration"] == pytest.approx(
            signal_plotter.all_df["duration"].median())


class TestTheReport:

    def test_html_report_is_written(self, tmp_path):
        from pycoQC.pycoQC import pycoQC
        path = write_signal_summary(tmp_path / "s.txt")
        out = tmp_path / "r.html"
        pycoQC(summary_file=str(path), html_outfile=str(out), quiet=True)
        page = out.read_text()
        assert "Signal summary" in page
        assert "Basecalled reads length" not in page

    def test_json_report_is_written(self, tmp_path):
        from pycoQC.pycoQC import pycoQC
        path = write_signal_summary(tmp_path / "s.txt")
        out = tmp_path / "r.json"
        pycoQC(summary_file=str(path), json_outfile=str(out), quiet=True)
        stats = json.loads(out.read_text())
        assert "signal" in stats["All Reads"]

    def test_pdf_report_is_written(self, tmp_path):
        pytest.importorskip("kaleido")
        pytest.importorskip("pypdf")
        import pypdf
        from pycoQC.pycoQC import pycoQC
        path = write_signal_summary(tmp_path / "s.txt")
        out = tmp_path / "r.pdf"
        pycoQC(summary_file=str(path), pdf_outfile=str(out), quiet=True)
        # A cover plus the signal plots.
        assert len(pypdf.PdfReader(str(out)).pages) >= 5
