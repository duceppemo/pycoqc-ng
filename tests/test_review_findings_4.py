# -*- coding: utf-8 -*-
"""Defects found in the fourth review pass.

As before, almost none of these failed loudly: they produced a report whose
numbers looked reasonable and were wrong.
"""

import datetime
import random

import numpy as np
import pandas as pd
import pytest

from pycoQC.common import pycoQCError
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC.pycoQC_plot import pycoQC_plot

GUPPY_SUMMARY = ("docs/pycoQC/data/"
                 "Guppy-2.2.4-basecall-1D-DNA_sequencing_summary+barcode.txt.gz")
ALBACORE_SUMMARY = ("docs/pycoQC/data/"
                    "Albacore-2.1.10_basecall-1D-DNA_sequencing_summary.txt.gz")


def write_summary(path, n=2000, run_id="run1", barcodes=True, seed=3):
    rng = random.Random(seed)
    df = pd.DataFrame({
        "read_id": [f"r{i:06d}" for i in range(n)],
        "run_id": run_id,
        "channel": [i % 512 + 1 for i in range(n)],
        "start_time": [float(i) for i in range(n)],
        "sequence_length_template": [rng.randint(200, 5000) for _ in range(n)],
        "mean_qscore_template": [round(rng.uniform(4, 14), 2) for _ in range(n)],
    })
    if barcodes:
        df["barcode_arrangement"] = [f"barcode{i % 4 + 1:02d}" for i in range(n)]
    df.to_csv(path, sep="\t", index=False)
    return path


class TestCalibrationFiltering:
    """The column is called calibration_strand_genome_template in Albacore
    and calibration_strand_genome from Guppy 2.2, and each writes its own
    markers for "this read did not align to a control"."""

    def test_the_guppy_column_is_recognised(self, tmp_path):
        df = pd.read_csv(GUPPY_SUMMARY, sep="\t")
        df.loc[:99, "calibration_strand_genome"] = "Lambda_3.6kb"
        path = tmp_path / "s.txt"
        df.to_csv(path, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(path), filter_calibration=True,
                              quiet=True)
        assert parser.counter["Calibration reads discarded"] == 100

    def test_guppy_markers_are_not_mistaken_for_control_genomes(self, tmp_path):
        """wrong_length and * mean no control was found. Treating them as
        genome names threw away 45,541 of the file's 50,000 reads."""
        parser = pycoQC_parse(summary_file=GUPPY_SUMMARY,
                              filter_calibration=True, quiet=True)
        assert parser.counter["Calibration reads discarded"] == 0
        assert len(parser.reads_df) == 50000

    def test_albacore_markers_are_still_understood(self):
        parser = pycoQC_parse(summary_file=ALBACORE_SUMMARY,
                              filter_calibration=True, quiet=True)
        assert parser.counter["Calibration reads discarded"] == 0

    def test_an_unknown_marker_stops_rather_than_emptying_the_run(self, tmp_path):
        df = pd.read_csv(GUPPY_SUMMARY, sep="\t", nrows=3000)
        df["calibration_strand_genome"] = "some_future_marker"
        path = tmp_path / "s.txt"
        df.to_csv(path, sep="\t", index=False)

        with pytest.raises(pycoQCError, match="not recognised"):
            pycoQC_parse(summary_file=str(path), filter_calibration=True,
                         quiet=True)


class TestBarcodesAreNames:
    """Deepbinner writes barcodes as bare numbers. Typed as numbers they
    split in two across chunks, became 6.0 once a merge introduced a missing
    value, and could not be sorted against "unclassified"."""

    @staticmethod
    def _deepbinner(path, read_ids):
        pd.DataFrame({
            "read_ID": read_ids,
            "barcode_call": [i % 6 + 1 for i in range(len(read_ids))],
        }).to_csv(path, sep="\t", index=False)
        return path

    def test_numeric_barcodes_are_read_as_labels(self, tmp_path):
        summary = write_summary(tmp_path / "s.txt", n=2000, barcodes=False)
        ids = pd.read_csv(summary, sep="\t", usecols=["read_id"])["read_id"]
        # Short of the summary, so the merge introduces missing values.
        barcodes = self._deepbinner(tmp_path / "bc.txt", ids[:1500])

        parser = pycoQC_parse(summary_file=str(summary),
                              barcode_file=str(barcodes), quiet=True)
        labels = sorted(parser.reads_df["barcode"].unique())
        assert labels == ["1", "2", "3", "4", "5", "6", "unclassified"]

    def test_split_names_its_files_after_the_barcode(self, tmp_path):
        from pycoQC.Barcode_split import barcode_split
        summary = write_summary(tmp_path / "s.txt", n=2000, barcodes=False)
        ids = pd.read_csv(summary, sep="\t", usecols=["read_id"])["read_id"]
        barcodes = self._deepbinner(tmp_path / "bc.txt", ids)
        out = tmp_path / "split"

        barcode_split(summary_file=str(summary), barcode_file=str(barcodes),
                      output_dir=str(out), quiet=True)

        assert sorted(p.name for p in out.iterdir()) == [
            f"sequencing_summary_{i}.txt" for i in range(1, 7)]


class TestElapsedTimeIsPerRun:
    """A basecaller's summary counts from the start of its own run. The
    converters counted from the first read in the directory, which charged
    the gap between two runs to the second one."""

    @staticmethod
    def _two_runs(path):
        def record(i, run_id, when):
            return ("@read_{}_{} runid={} ch={} start_time={}\n{}\n+\n{}\n"
                    .format(run_id, i, run_id, i % 512 + 1,
                            when.strftime("%Y-%m-%dT%H:%M:%SZ"),
                            "ACGT" * 50, "I" * 200))
        with open(path, "w") as fh:
            for i in range(200):
                fh.write(record(i, "AAAA",
                                datetime.datetime(2024, 1, 1)
                                + datetime.timedelta(seconds=18 * i)))
            for i in range(200):
                # Two days later: a reloaded flow cell, or one parent folder.
                fh.write(record(i, "BBBB",
                                datetime.datetime(2024, 1, 3)
                                + datetime.timedelta(seconds=18 * i)))

    def test_each_run_starts_at_zero(self, tmp_path):
        from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary
        d = tmp_path / "fastq"
        d.mkdir()
        self._two_runs(d / "reads.fastq")
        out = tmp_path / "s.txt"

        Fastq_to_seq_summary(fastq_dir=str(d), seq_summary_fn=str(out),
                             threads=3)

        df = pd.read_csv(out, sep="\t")
        spans = df.groupby("run_id")["start_time"].agg(["min", "max"])
        assert (spans["min"] == 0).all()
        assert (spans["max"] == 3582).all()

    def test_each_run_keeps_its_own_start(self, tmp_path):
        from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary
        d = tmp_path / "fastq"
        d.mkdir()
        self._two_runs(d / "reads.fastq")
        out = tmp_path / "s.txt"

        Fastq_to_seq_summary(fastq_dir=str(d), seq_summary_fn=str(out),
                             threads=3)

        df = pd.read_csv(out, sep="\t")
        starts = df.groupby("run_id")["run_start_time"].first()
        assert starts["AAAA"].startswith("2024-01-01")
        assert starts["BBBB"].startswith("2024-01-03")


class TestHistogramBins:

    def test_the_longest_and_shortest_reads_are_inside_the_bins(self):
        """Going through log10 and back moved the bottom edge up, so the
        shortest read in the run fell outside the histogram."""
        bins = pycoQC_plot._hist_bins(5, 183021, 200, "log")
        counted, _ = np.histogram([5, 183021], bins=bins)
        assert counted.sum() == 2

    def test_the_axis_stops_at_the_data(self):
        """The top edge was extended by 0.1 in log10, a quarter of the way
        to the next decade, so the density ran out well past the longest
        read in the run."""
        bins = pycoQC_plot._hist_bins(5, 100946, 200, "log")
        assert bins[-1] == pytest.approx(100946, rel=1e-6)

    def test_identical_values_still_have_a_bin_to_fall_in(self):
        for scale in ("log", "linear"):
            bins = pycoQC_plot._hist_bins(500, 500, 50, scale)
            counted, _ = np.histogram([500] * 20, bins=bins)
            assert counted.sum() == 20

    def test_an_all_missing_column_gives_an_empty_histogram(self):
        """A BAM aligned against the wrong reference leaves every align_len
        missing; nanmin of that is nan, the edges became nan, and numpy
        counted the missing values into the top bins."""
        x, y = pycoQC_plot._compute_hist(pd.Series([np.nan] * 3), "log")
        assert (x, y) == ([], [])


class TestDensitiesAreCentredOnTheirBins:

    @pytest.fixture(scope="class")
    @classmethod
    def plotter(cls, tmp_path_factory):
        path = write_summary(tmp_path_factory.mktemp("d") / "s.txt", n=3000)
        return pycoQC_plot(pycoQC_parse(summary_file=str(path), quiet=True),
                           quiet=True)

    def test_the_1d_density_stops_at_the_longest_read(self, plotter):
        x = np.asarray(plotter.read_len_1D()["data"][0]["x"], dtype=float)
        assert x[-1] < plotter.all_df["read_len"].max()

    def test_the_2d_density_grid_matches_its_cells(self, plotter):
        """A contour trace reads the coordinates it is given as cell
        centres, so handing it the edges drew the density half a bin away
        from the median marker on the same axes."""
        trace = plotter.read_len_read_qual_2D()["data"][0]
        z = np.asarray(trace["z"])
        assert len(trace["x"]) == z.shape[1]
        assert len(trace["y"]) == z.shape[0]


class TestReportedTotalsAreTheRunsOwn:

    @pytest.fixture(scope="class")
    @classmethod
    def plotter(cls, tmp_path_factory):
        path = write_summary(tmp_path_factory.mktemp("t") / "s.txt", n=4000)
        parser = pycoQC_parse(summary_file=str(path), quiet=True)
        # Sample hard, so anything taken from the sample is obviously wrong.
        return pycoQC_plot(parser, sample=400, quiet=True)

    @pytest.mark.parametrize("count_level", ["reads", "bases"])
    def test_output_over_time_ends_at_the_real_total(self, plotter, count_level):
        """Bases were scaled by the ratio of read counts, which the sampled
        reads do not obey: the yield on the plot contradicted the summary
        table on the same page."""
        _, data, _ = plotter._pycoQC_plot__output_over_time_data(
            "all", count_level)
        exact = (len(plotter.all_df) if count_level == "reads"
                 else int(plotter.all_df["read_len"].sum()))
        assert np.asarray(data["y"][0])[-1] == pytest.approx(exact)

    def test_the_barcode_table_counts_the_whole_run(self, plotter):
        cells = plotter.basecall_summary_per_barcode()["data"][0]["cells"]["values"]
        assert sum(cells[1]) == len(plotter.all_df)
        assert sum(cells[2]) == int(plotter.all_df["read_len"].sum())

    def test_pass_and_fail_barcode_counts_add_up(self, plotter):
        """Pass and fail are sampled separately, so their scaled counts did
        not add up to the all reads row for any barcode."""
        def counts(level):
            _, data = plotter._pycoQC_plot__summary_barcode_data(level)
            return {row[0]: row[1] for row in data["y"][0]}

        all_c, pass_c, fail_c = counts("all"), counts("pass"), counts("fail")
        for bc in all_c:
            assert pass_c.get(bc, 0) + fail_c.get(bc, 0) == all_c[bc]

    def test_the_barcode_pie_counts_the_whole_run(self, plotter):
        _, data = plotter._pycoQC_plot__barcode_counts_data("all")
        assert sum(data["values"][0]) == len(plotter.all_df)


class TestOptionalColumnsDoNotCostReads:
    """run_id, channel and start_time decide which plots a read can appear
    in. Requiring them discarded reads that had a length and a quality, and
    those are the numbers people quote."""

    @pytest.fixture
    def gappy(self, tmp_path):
        path = write_summary(tmp_path / "s.txt", n=3000)
        df = pd.read_csv(path, sep="\t")
        df.loc[::10, "channel"] = None
        df.loc[1::25, "start_time"] = None
        df.to_csv(path, sep="\t", index=False)
        return path

    def test_reads_missing_a_channel_are_kept(self, gappy):
        parser = pycoQC_parse(summary_file=str(gappy), quiet=True)
        assert len(parser.reads_df) == 3000

    @pytest.mark.parametrize("method", ["output_over_time", "channels_activity",
                                        "read_len_over_time", "read_len_1D"])
    def test_the_plots_still_build(self, gappy, method):
        parser = pycoQC_parse(summary_file=str(gappy), quiet=True)
        plotter = pycoQC_plot(parser, quiet=True)
        assert getattr(plotter, method)()["data"]


class TestOneSpellingPerColumn:
    """A file carrying two spellings of one field left the frame with two
    columns of the same name, and pandas then answered an ordinary question
    about it with "The truth value of a Series is ambiguous"."""

    @pytest.mark.parametrize("extra,value", [
        ("barcode", "legacy"),
        ("sequence_length", 1),
        ("mean_qscore", 0.0),
        ("sequence_length_2d", 1),
    ])
    def test_the_preferred_spelling_wins(self, tmp_path, extra, value):
        df = pd.read_csv(GUPPY_SUMMARY, sep="\t", nrows=3000)
        expected_len = df["sequence_length_template"].median()
        df[extra] = value
        path = tmp_path / "s.txt"
        df.to_csv(path, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(path), quiet=True)
        assert parser.reads_df["read_len"].median() == expected_len
        assert "legacy" not in set(parser.reads_df["barcode"])


class TestInputProblemsAreExplained:

    def test_an_empty_file_says_so(self, tmp_path):
        path = tmp_path / "empty.txt"
        path.touch()
        with pytest.raises(pycoQCError, match="No data found"):
            pycoQC_parse(summary_file=str(path), quiet=True)

    def test_columns_missing_from_one_file_are_reported(self, tmp_path, caplog):
        """Only the columns every file has survive the stack, which used to
        happen in silence."""
        full = write_summary(tmp_path / "a.txt", n=2000)
        partial = tmp_path / "b.txt"
        pd.read_csv(full, sep="\t").drop(columns="start_time").to_csv(
            partial, sep="\t", index=False)

        with caplog.at_level("WARNING"):
            pycoQC_parse(summary_file=[str(full), str(partial)])

        assert "start_time" in caplog.text

    def test_a_single_read_file_says_what_is_wrong(self, tmp_path):
        path = write_summary(tmp_path / "s.txt", n=1)
        with pytest.raises(pycoQCError, match="at least two reads"):
            pycoQC_parse(summary_file=str(path), quiet=True)


class TestOneMeanCoverage:

    def test_the_table_and_the_plot_agree(self, tmp_path):
        """The table divided aligned query bases by the reference length
        while the plot integrates reference span, so one report gave two
        different mean coverages."""
        pysam = pytest.importorskip("pysam")
        ref_len = 60000
        rng = random.Random(7)
        unsorted = tmp_path / "unsorted.bam"
        header = {"HD": {"VN": "1.6", "SO": "unknown"},
                  "SQ": [{"SN": "chr1", "LN": ref_len}]}
        n = 300
        with pysam.AlignmentFile(str(unsorted), "wb", header=header) as bam:
            for i in range(n):
                rec = pysam.AlignedSegment()
                rec.query_name = f"r{i:06d}"
                rec.query_sequence = "".join(rng.choice("ACGT")
                                             for _ in range(1000))
                rec.query_qualities = pysam.qualitystring_to_array("I" * 1000)
                rec.flag = 0
                rec.reference_id = 0
                rec.reference_start = rng.randint(0, ref_len - 1600)
                rec.mapping_quality = 60
                # Insertions and deletions, so the two definitions differ.
                rec.cigartuples = [(0, 400), (1, 100), (0, 400), (2, 300),
                                   (0, 100)]
                rec.set_tag("NM", 400, "i")
                bam.write(rec)
        final = tmp_path / "aligned.bam"
        pysam.sort("-o", str(final), str(unsorted))
        pysam.index(str(final))

        summary = write_summary(tmp_path / "s.txt", n=n, barcodes=False)

        parser = pycoQC_parse(summary_file=str(summary), bam_file=str(final),
                              quiet=True)
        plotter = pycoQC_plot(parser, quiet=True)
        table = plotter.alignment_summary()["data"][0]["cells"]["values"][3][0]
        plot = plotter.alignment_coverage()["data"][1]["y"][0]
        assert table == pytest.approx(plot, rel=0.01)


class TestBarcodeFigureControls:
    """The three read-set buttons sat over the pie and the legend, only
    ever restyled the table, and because each carried both args and args2
    plotly treated them as toggles, so a switch took two clicks."""

    @pytest.fixture
    def plotter(self, tmp_path):
        path = write_summary(tmp_path / "s.txt", n=4000)
        return pycoQC_plot(pycoQC_parse(summary_file=str(path), quiet=True),
                           quiet=True)

    def test_one_click_switches_read_set(self, plotter):
        menu = plotter.barcode_counts()["layout"]["updatemenus"][0]
        for button in menu["buttons"]:
            assert getattr(button, "args2", None) is None

    def test_both_the_table_and_the_pie_follow_the_button(self, plotter):
        fig = plotter.barcode_counts()
        kinds = [trace.type for trace in fig["data"]]
        assert kinds == ["table", "pie"] * 3
        for i, button in enumerate(fig["layout"]["updatemenus"][0]["buttons"]):
            visible = button["args"][0]["visible"]
            assert [n for n, v in enumerate(visible) if v] == [2 * i, 2 * i + 1]

    def test_each_read_set_shows_its_own_counts(self, plotter):
        fig = plotter.barcode_counts()
        totals = [sum(fig["data"][2 * i]["cells"]["values"][1])
                  for i in range(3)]
        assert totals == [len(plotter.all_df), len(plotter.pass_df),
                          len(plotter.fail_df)]

    def test_the_buttons_sit_below_the_pie(self, plotter):
        menu = plotter.barcode_counts()["layout"]["updatemenus"][0]
        assert menu["direction"] == "right"
        assert menu["y"] <= 0          # under the plotting area
        assert 0 <= menu["x"] <= 0.4   # across the pie, not over the legend

    def test_the_figure_grows_with_the_number_of_barcodes(self, tmp_path):
        """At a fixed height a 96 barcode run had its legend and its table
        cut off partway down."""
        path = tmp_path / "many.txt"
        n = 4800
        pd.DataFrame({
            "read_id": [f"r{i:06d}" for i in range(n)],
            "run_id": "run1",
            "channel": [i % 512 + 1 for i in range(n)],
            "start_time": [float(i) for i in range(n)],
            "sequence_length_template": 1000,
            "mean_qscore_template": 10.0,
            "barcode_arrangement": [f"barcode{i % 96 + 1:02d}"
                                    for i in range(n)],
        }).to_csv(path, sep="\t", index=False)

        plotter = pycoQC_plot(
            pycoQC_parse(summary_file=str(path), min_barcode_percent=0,
                         quiet=True), quiet=True)

        for fig in (plotter.barcode_counts(),
                    plotter.basecall_summary_per_barcode()):
            # Room for a header and 96 rows.
            assert fig["layout"]["height"] >= 96 * 30

    def test_an_explicit_height_is_respected(self, plotter):
        assert plotter.barcode_counts(height=640)["layout"]["height"] == 640


class TestTheLargestBarcodeKit:
    """96 is as many as a standard Oxford Nanopore kit goes (SQK-NBD114-96,
    SQK-RBK114-96, EXP-PBC096), so the report has to hold that many."""

    @pytest.fixture(scope="class")
    @classmethod
    def plotter(cls, tmp_path_factory):
        path = tmp_path_factory.mktemp("bc96") / "s.txt"
        n = 9600
        pd.DataFrame({
            "read_id": [f"r{i:06d}" for i in range(n)],
            "run_id": "run1",
            "channel": [i % 512 + 1 for i in range(n)],
            "start_time": [float(i) for i in range(n)],
            "sequence_length_template": [500 + (i % 4000) for i in range(n)],
            "mean_qscore_template": [5.0 + (i % 10) for i in range(n)],
            "barcode_arrangement": [f"barcode{i % 96 + 1:02d}"
                                    for i in range(n)],
        }).to_csv(path, sep="\t", index=False)
        parser = pycoQC_parse(summary_file=str(path), min_barcode_percent=0,
                              quiet=True)
        return pycoQC_plot(parser, quiet=True)

    def test_every_barcode_reaches_the_table(self, plotter):
        cells = plotter.basecall_summary_per_barcode()["data"][0]["cells"]["values"]
        assert len(cells[0]) == 96

    def test_the_table_is_tall_enough_to_show_them(self, plotter):
        fig = plotter.basecall_summary_per_barcode()
        # Header plus 96 rows of 30 pixels.
        assert fig["layout"]["height"] >= 40 + 96 * 30

    def test_counts_are_drawn_as_bars_not_a_pie(self, plotter):
        """A pie divided 96 ways is a ring of slivers under a thicket of
        leader lines, with a legend longer than the plot is tall."""
        kinds = [trace.type for trace in plotter.barcode_counts()["data"]]
        assert kinds == ["table", "bar"] * 3

    def test_every_barcode_gets_a_bar(self, plotter):
        bar = plotter.barcode_counts()["data"][1]
        assert len(bar["y"]) == 96
        assert sum(bar["x"]) == len(plotter.all_df)

    def test_a_small_run_keeps_its_pie(self, tmp_path):
        path = write_summary(tmp_path / "s.txt", n=2000)   # four barcodes
        plotter = pycoQC_plot(pycoQC_parse(summary_file=str(path), quiet=True),
                              quiet=True)
        kinds = [trace.type for trace in plotter.barcode_counts()["data"]]
        assert kinds == ["table", "pie"] * 3

    def test_the_boxplots_cover_every_barcode(self, plotter):
        for method in ("length_boxplot_per_barcode", "gc_boxplot_per_barcode"):
            try:
                fig = getattr(plotter, method)()
            except pycoQCError:
                continue          # %GC is not in this summary
            # One box trace, grouped by the barcode on the x axis.
            assert len(set(fig["data"][0]["x"])) == 96


class TestPdfPagesFitTheirPlot:

    def test_a_tall_plot_gets_a_tall_page(self, tmp_path):
        """Every page was rendered at one size, so a 96 barcode table was
        cut off after its twelfth row."""
        pytest.importorskip("pypdf")
        pytest.importorskip("kaleido")
        from pycoQC.pycoQC_report import pycoQC_report
        import pypdf

        path = tmp_path / "s.txt"
        n = 4800
        pd.DataFrame({
            "read_id": [f"r{i:06d}" for i in range(n)],
            "run_id": "run1",
            "channel": [i % 512 + 1 for i in range(n)],
            "start_time": [float(i) for i in range(n)],
            "sequence_length_template": [500 + (i % 4000) for i in range(n)],
            "mean_qscore_template": [5.0 + (i % 10) for i in range(n)],
            "barcode_arrangement": [f"barcode{i % 96 + 1:02d}"
                                    for i in range(n)],
        }).to_csv(path, sep="\t", index=False)

        parser = pycoQC_parse(summary_file=str(path), min_barcode_percent=0,
                              quiet=True)
        report = pycoQC_report(parser, pycoQC_plot(parser, quiet=True),
                               quiet=True)
        out = tmp_path / "r.pdf"
        report.pdf_report(str(out), height=600)

        heights = [float(page.mediabox.height)
                   for page in pypdf.PdfReader(str(out)).pages]
        assert max(heights) > 2 * min(heights)
