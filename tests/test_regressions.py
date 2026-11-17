# -*- coding: utf-8 -*-
"""Regressions found running every entry point against real runs.

Albacore, Guppy and MinKNOW summaries, Dorado BAM and FASTQ, a PromethION
POD5 and MinKNOW multi-read FAST5 files were each taken through conversion
and reporting. Each class below covers one thing that broke, or one thing a
restructuring of the code must keep doing.
"""

import gzip
import json
import logging
import os
import subprocess
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd
import pytest

from conftest import REPO_ROOT, _write_basecalled_fast5, _write_pod5
from pycoQC.__main__ import main_pycoQC
from pycoQC.Bam_to_seq_summary import Bam_to_seq_summary
from pycoQC.Barcode_split import barcode_split
from pycoQC.common import (
    check_arg, expand_file_names, format_worker_error, head, pycoQCError)
from pycoQC.Fast5_to_seq_summary import Fast5_to_seq_summary
from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary
from pycoQC.Pod5_to_seq_summary import Pod5_to_seq_summary
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC.pycoQC_plot import COVERAGE_MAX_LABELLED_REFS, pycoQC_plot
from pycoQC.pycoQC_report import pycoQC_report


def write_summary(path, n=600, **columns):
    """A basecalled summary, with any column overridden by keyword."""
    data = {
        "read_id": [f"r{i:05d}" for i in range(n)],
        "run_id": "run1",
        "channel": [i % 512 + 1 for i in range(n)],
        "start_time": [float(i * 10) for i in range(n)],
        "sequence_length_template": [200 + (i * 37) % 5000 for i in range(n)],
        "mean_qscore_template": [4.0 + (i % 12) for i in range(n)],
    }
    data.update(columns)
    pd.DataFrame(data).to_csv(path, sep="\t", index=False)
    return path


def plotter_for(path, **kwargs):
    return pycoQC_plot(pycoQC_parse(summary_file=str(path), quiet=True),
                       quiet=True, **kwargs)


@pytest.fixture(scope="module")
def all_pass(tmp_path_factory):
    """A run where every read passes, leaving the fail set empty."""
    path = write_summary(tmp_path_factory.mktemp("allpass") / "s.txt")
    plotter = plotter_for(path, min_pass_qual=0)
    assert len(plotter.fail_df) == 0
    return plotter


@pytest.fixture(scope="module")
def box_plotter(tmp_path_factory):
    names = ["barcode2", "barcode10", "barcode1", "unclassified"]
    path = write_summary(tmp_path_factory.mktemp("box") / "s.txt",
                         barcode_arrangement=[names[i % 4] for i in range(600)])
    return plotter_for(path)


class TestMultiReadFast5WithoutFileType:
    """MinKNOW wrote multi-read FAST5 before it added the file_type
    attribute. Those files were taken for single read ones and the
    conversion died looking for /Raw/Reads."""

    @pytest.fixture
    def fast5_dir(self, tmp_path):
        h5py = pytest.importorskip("h5py")
        path = _write_basecalled_fast5(tmp_path / "batch_0.fast5", n_reads=12)
        with h5py.File(path, "a") as fh:
            del fh.attrs["file_type"]
        return tmp_path

    def test_every_read_is_converted(self, fast5_dir, tmp_path):
        out = tmp_path / "summary.txt"
        Fast5_to_seq_summary(fast5_dir=str(fast5_dir), seq_summary_fn=str(out),
                             threads=3)
        df = pd.read_csv(out, sep="\t")
        assert len(df) == 12
        assert df["read_id"].str.endswith("-b178-468d-ab62-af921129c079").all()

    def test_a_file_without_reads_says_so(self, tmp_path):
        h5py = pytest.importorskip("h5py")
        d = tmp_path / "empty"
        d.mkdir()
        with h5py.File(d / "empty.fast5", "w") as fh:
            fh.attrs["file_type"] = "multi-read"
        with pytest.raises(pycoQCError, match="No reads found in the fast5"):
            Fast5_to_seq_summary(fast5_dir=str(d),
                                 seq_summary_fn=str(tmp_path / "s.txt"), threads=3)

    def test_the_reason_is_not_buried_in_a_traceback(self, tmp_path):
        """A pycoQCError raised in a worker reaches the caller as its own
        message, not wrapped in thirty lines of the worker's stack."""
        d = tmp_path / "nothing"
        d.mkdir()
        with pytest.raises(pycoQCError) as excinfo:
            Fast5_to_seq_summary(fast5_dir=str(d),
                                 seq_summary_fn=str(tmp_path / "s.txt"), threads=3)
        assert "Traceback" not in str(excinfo.value)
        assert "No valid fast5" in str(excinfo.value)


class TestEmptyReadSets:
    """A run with no failing reads sent an empty update for the Fail Reads
    button, and plotly then kept drawing whatever was on screen before:
    the passing reads, labelled as the failing ones."""

    @pytest.mark.parametrize("method", [
        "read_len_1D", "read_qual_1D", "read_len_read_qual_2D",
        "read_len_over_time", "read_qual_over_time", "output_over_time",
        "channels_activity"])
    def test_every_button_carries_a_full_update(self, all_pass, method):
        fig = getattr(all_pass, method)()
        for button in fig.layout.updatemenus[0].buttons:
            assert button.args[0], (method, button.label)
        fail = [b for b in fig.layout.updatemenus[0].buttons if b.label.startswith("Fail")]
        assert fail
        for button in fail:
            for key in ("x", "y"):
                assert all(len(v) == 0 for v in button.args[0][key] if v is not None), \
                    (method, button.label, key)

    def test_a_field_with_no_values_is_skipped_not_drawn_empty(self, all_pass):
        """With nothing to plot at all the plot says so, and the report
        leaves it out, instead of raising KeyError on the empty data."""
        all_pass.all_sample_df = all_pass.all_sample_df.assign(mean_qscore=np.nan)
        try:
            with pytest.raises(pycoQCError):
                all_pass.read_qual_1D()
        finally:
            all_pass.all_sample_df = all_pass.all_df

    def test_a_flat_2d_density_still_has_contour_levels(self, tmp_path):
        path = write_summary(tmp_path / "s.txt", n=50,
                             sequence_length_template=[1000] * 50,
                             mean_qscore_template=[10.0] * 50)
        fig = plotter_for(path).read_len_read_qual_2D()
        assert fig.data[0].contours.size > 0


class TestUnbasecalledRunDescription:
    """`pycoQC -v` printed the plotter, which looked up read_len, so a POD5
    or MinKNOW summary dumped a traceback in the middle of the run."""

    def test_str_describes_a_signal_only_run(self, tmp_path):
        path = tmp_path / "s.txt"
        pd.DataFrame({"read_id": ["a", "b", "c"], "channel": [1, 2, 3],
                      "start_time": [0, 5, 9], "duration": [1.0, 2.0, 3.0]}
                     ).to_csv(path, sep="\t", index=False)
        text = str(plotter_for(path))
        assert "All reads: 3" in text
        assert "bases" not in text


class TestChannelActivity:

    def test_a_channel_past_the_platform_count_widens_the_map(self, tmp_path):
        """IndexError used to take down the whole report."""
        channels = [i % 400 + 1 for i in range(600)]
        channels[10] = 3100
        fig = plotter_for(write_summary(tmp_path / "s.txt", channel=channels)
                          ).channels_activity()
        assert len(fig.data[0].x) == 3100

    def test_channel_zero_is_not_drawn_on_the_last_channel(self, tmp_path):
        channels = [i % 400 + 1 for i in range(600)]
        channels[:50] = [0] * 50
        fig = plotter_for(write_summary(tmp_path / "s.txt", channel=channels)
                          ).channels_activity(smooth_sigma=0)
        z = np.asarray(fig.data[0].z)
        assert z.shape[1] == 512
        assert z[:, -1].sum() == 0
        assert z.sum() == 550

    def test_counts_match_a_plain_loop(self, tmp_path):
        rng = np.random.default_rng(3)
        path = write_summary(tmp_path / "s.txt",
                             channel=list(rng.integers(1, 513, 600)),
                             start_time=list(rng.uniform(0, 7200, 600)))
        plotter = plotter_for(path)
        _, data = plotter._channels_activity_data(
            "all", count_level="bases", n_channels=512, smooth_sigma=0, time_bins=20)
        df = plotter.all_sample_df
        t = df["start_time"].to_numpy(dtype=np.float64) / 3600
        edges = np.linspace(t.min(), t.max(), 21)
        idx = np.clip(np.digitize(t, edges[1:-1]), 0, 19)
        expected = np.zeros((20, 512))
        for i, ch, n in zip(idx, df["channel"], df["read_len"]):
            expected[i, ch - 1] += n
        np.testing.assert_allclose(data["z"][0], expected, rtol=1e-6)


class TestPerBarcodeBoxplot:

    def test_barcodes_sort_as_numbers(self, box_plotter):
        """category ascending sorted as text: barcode10 before barcode2."""
        fig = box_plotter.length_boxplot_per_barcode()
        assert list(fig.layout.xaxis.categoryarray) == [
            "barcode1", "barcode2", "barcode10", "unclassified"]

    def test_boxes_are_the_quartiles_of_the_reads(self, box_plotter):
        plotter = box_plotter
        fig = plotter.length_boxplot_per_barcode()
        box = fig.data[0]
        df = plotter.all_sample_df
        for i, bc in enumerate(box.x):
            values = df.loc[df["barcode"] == bc, "read_len"].to_numpy(dtype=float)
            q1, med, q3 = np.percentile(values, (25, 50, 75))
            assert box.q1[i] == pytest.approx(q1)
            assert box.median[i] == pytest.approx(med)
            assert box.q3[i] == pytest.approx(q3)
            assert box.lowerfence[i] >= values.min()
            assert box.upperfence[i] <= values.max()

    def test_the_reads_are_not_shipped_in_the_report(self, box_plotter):
        """The sample went into the page four times, 6 MB per boxplot on a
        real run."""
        fig = box_plotter.length_boxplot_per_barcode()
        assert fig.data[0].y is None
        assert len(fig.to_json()) < 20000

    def test_boxplots_have_help(self):
        for method in (pycoQC_plot.length_boxplot_per_barcode,
                       pycoQC_plot.gc_boxplot_per_barcode):
            assert method.__doc__ and "* plot_title" in method.__doc__


class TestAlignmentColours:

    def test_a_category_keeps_its_colour_when_another_is_absent(self, tmp_path, aligned_bam):
        """Colours went by position, so Unmapped took Supplementary's colour
        on a run without supplementary alignments."""
        summary = write_summary(tmp_path / "s.txt")
        parser = pycoQC_parse(summary_file=str(summary), bam_file=str(aligned_bam), quiet=True)
        plotter = pycoQC_plot(parser, quiet=True)
        colours = ["p", "s", "x", "d", "u"]
        full = dict(zip(plotter.alignments_df["Alignments"],
                        plotter.alignment_reads_status(colors=colours).data[1].marker.colors))
        assert full == {"Primary": "p", "Secondary": "s", "Supplementary": "x", "Unmapped": "u"}
        plotter.alignments_df = plotter.alignments_df[
            plotter.alignments_df["Alignments"] != "Supplementary"]
        pie = plotter.alignment_reads_status(colors=colours).data[1]
        assert dict(zip(pie.labels, pie.marker.colors))["Unmapped"] == "u"


class TestCoverageReferences:

    @pytest.fixture
    def plotter(self, tmp_path, aligned_bam):
        summary = write_summary(tmp_path / "s.txt")
        parser = pycoQC_parse(summary_file=str(summary), bam_file=str(aligned_bam),
                              quiet=True)
        return pycoQC_plot(parser, quiet=True)

    def test_the_last_of_an_even_number_of_references_is_shaded(self, plotter):
        plotter.ref_len_dict = OrderedDict(
            [("chr1", 50000), ("chr2", 50000), ("chr3", 50000), ("chr4", 50000)])
        fig = plotter.alignment_coverage(nbins=400)
        bounds = sorted((s.x0, s.x1) for s in fig.layout.shapes)
        assert bounds == [(100, 200), (300, 400)]

    def test_many_references_are_not_labelled_one_by_one(self, plotter):
        """A fragmented assembly with 233,000 contigs made a 31 MB plot."""
        n = COVERAGE_MAX_LABELLED_REFS * 10
        refs = OrderedDict([("chr1", 50000)])
        refs.update((f"contig{i}", 1000) for i in range(n))
        plotter.ref_len_dict = refs
        fig = plotter.alignment_coverage()
        assert not fig.layout.shapes
        assert not fig.layout.xaxis.ticktext
        assert "{:,} references".format(n + 1) in fig.layout.xaxis.title.text


class TestParsing:

    def test_runid_list_applies_without_start_times(self, tmp_path):
        """The selection sat behind the check for start_time, so a summary
        without one was reported in full whatever runid_list said."""
        path = tmp_path / "s.txt"
        pd.DataFrame({"read_id": [f"r{i}" for i in range(40)],
                      "run_id": ["a", "b"] * 20,
                      "sequence_length_template": [1000] * 40,
                      "mean_qscore_template": [10.0] * 40}
                     ).to_csv(path, sep="\t", index=False)
        parser = pycoQC_parse(summary_file=str(path), runid_list=["a"], quiet=True)
        assert set(parser.reads_df["run_id"]) == {"a"}
        assert len(parser.reads_df) == 20

    def test_runid_list_needs_a_run_id_column(self, tmp_path):
        path = tmp_path / "s.txt"
        pd.DataFrame({"read_id": ["x", "y", "z"],
                      "sequence_length_template": [10, 20, 30],
                      "mean_qscore_template": [9.0, 9.0, 9.0]}
                     ).to_csv(path, sep="\t", index=False)
        with pytest.raises(pycoQCError, match="no run_id"):
            pycoQC_parse(summary_file=str(path), runid_list=["a"], quiet=True)

    def test_an_unrecognised_barcode_file_is_named(self, tmp_path):
        summary = write_summary(tmp_path / "s.txt")
        with pytest.raises(pycoQCError, match="not in a recognised format"):
            pycoQC_parse(summary_file=str(summary), barcode_file=str(summary), quiet=True)

    def test_a_plain_read_id_and_barcode_file_is_accepted(self, tmp_path):
        """The check for unknown layouts must not reject the simplest one."""
        summary = write_summary(tmp_path / "s.txt")
        bc = tmp_path / "bc.tsv"
        pd.DataFrame({"read_id": [f"r{i:05d}" for i in range(600)],
                      "barcode": ["barcode01", "barcode02"] * 300}
                     ).to_csv(bc, sep="\t", index=False)
        parser = pycoQC_parse(summary_file=str(summary), barcode_file=str(bc), quiet=True)
        assert parser.reads_df["barcode"].value_counts().to_dict() == {
            "barcode01": 300, "barcode02": 300}

    def test_alignment_categories_come_in_a_fixed_order(self, tmp_path, aligned_bam):
        summary = write_summary(tmp_path / "s.txt")
        parser = pycoQC_parse(summary_file=str(summary), bam_file=str(aligned_bam),
                              quiet=True)
        assert list(parser.alignments_df["Alignments"]) == [
            "Primary", "Secondary", "Supplementary", "Unmapped"]

    def test_an_n_in_md_is_a_mismatch(self):
        pysam = pytest.importorskip("pysam")
        header = pysam.AlignmentHeader.from_dict(
            {"SQ": [{"SN": "chr1", "LN": 1000}]})
        rec = pysam.AlignedSegment(header)
        rec.query_name = "r"
        rec.query_sequence = "A" * 10
        rec.reference_id = 0
        rec.reference_start = 0
        rec.cigartuples = [(0, 10)]
        rec.set_tag("MD", "5N4", "Z")
        stats = pycoQC_parse._get_read_stats(None, rec)
        assert stats["mismatch"] == 1
        assert stats["identity_freq"] == pytest.approx(0.9)

    def test_quiet_means_no_progress(self, tmp_path, caplog):
        """Progress was logged as warnings, so quiet could not hide it."""
        path = write_summary(tmp_path / "s.txt")
        with caplog.at_level(logging.INFO):
            pycoQC_plot(pycoQC_parse(summary_file=str(path), quiet=True), quiet=True)
        assert "Check input data files" not in caplog.text
        assert "Loading plotting interface" not in caplog.text


class TestFileNames:

    def test_a_name_with_brackets_is_taken_literally(self, tmp_path):
        """glob read [1] as a character class, so this file could not be
        named at all."""
        path = write_summary(tmp_path / "run[1].txt")
        assert expand_file_names(str(path)) == [str(path)]

    def test_a_file_named_twice_is_read_once(self, tmp_path):
        """Reading it twice doubled every count."""
        path = write_summary(tmp_path / "s.txt")
        assert expand_file_names([str(path), str(tmp_path / "*.txt")]) == [str(path)]
        parser = pycoQC_parse(summary_file=[str(path), str(path)], quiet=True)
        assert len(parser.reads_df) == 600

    def test_glob_matches_are_sorted(self, tmp_path):
        for name in ("c.txt", "a.txt", "b.txt"):
            write_summary(tmp_path / name, n=3)
        found = expand_file_names(str(tmp_path / "*.txt"))
        assert [os.path.basename(f) for f in found] == ["a.txt", "b.txt", "c.txt"]

    def test_head_reads_gzipped_files(self, tmp_path, capsys):
        path = tmp_path / "s.txt.gz"
        with gzip.open(path, "wt") as fh:
            fh.write("a\tb\n1\t2\n")
        head(str(path))
        assert "a" in capsys.readouterr().out


class TestErrors:

    def test_check_arg_raises_the_package_error(self):
        with pytest.raises(pycoQCError):
            check_arg("min_pass_len", "many", int)
        with pytest.raises(pycoQCError):
            check_arg("min_pass_len", -1, int, min=0)

    def test_worker_errors_keep_only_the_message_for_input_problems(self):
        try:
            raise pycoQCError("bad input")
        except pycoQCError as E:
            assert format_worker_error(E) == "bad input"
        try:
            raise KeyError("bug")
        except KeyError as E:
            assert "Traceback" in format_worker_error(E)

    def test_cli_prints_the_reason_and_exits_1(self, tmp_path, capsys):
        with pytest.raises(SystemExit) as excinfo:
            main_pycoQC(["-f", str(tmp_path / "missing.txt"), "-j", str(tmp_path / "o.json")])
        assert excinfo.value.code == 1
        err = capsys.readouterr().err
        assert err.startswith("ERROR: No files found")
        assert "Traceback" not in err

    def test_cli_selects_runs(self, tmp_path):
        """--runid_list was in the Python API and missing from the command."""
        path = write_summary(tmp_path / "s.txt",
                             run_id=["a", "b", "c"] * 200)
        out = tmp_path / "o.json"
        main_pycoQC(["-f", str(path), "-j", str(out), "--runid_list", "a", "c", "-q"])
        stats = json.loads(out.read_text())
        assert stats["All Reads"]["basecall"]["reads_number"] == 400
        assert stats["All Reads"]["run"]["runid_number"] == 2


class TestReport:

    @pytest.fixture
    def reporter(self, tmp_path):
        path = write_summary(tmp_path / "s.txt", sample_id=["<b>A&B</b>"] * 600)
        parser = pycoQC_parse(summary_file=str(path), quiet=True)
        return pycoQC_report(parser, pycoQC_plot(parser, quiet=True), quiet=True)

    def test_run_details_are_escaped(self, reporter, tmp_path):
        """The template's name ends in .j2, which select_autoescape(["html"])
        does not match, so input values went into the page as markup."""
        out = tmp_path / "r.html"
        reporter.html_report(str(out))
        html = out.read_text(encoding="utf-8")
        assert "&lt;b&gt;A&amp;B&lt;/b&gt;" in html
        assert "<b>A&B</b>" not in html
        # The generated parts are still raw markup.
        assert "<script>" in html and "Plotly" in html
        assert "<li>{}</li>".format(os.path.abspath(reporter.parser.summary_files_list[0])) in html

    def test_non_plot_methods_in_a_config_are_skipped(self, reporter, tmp_path):
        config = tmp_path / "c.json"
        config.write_text(json.dumps({"summary_stats_dict": {}, "_compute_stats": {},
                                      "run_summary": {}}))
        out = tmp_path / "r.html"
        reporter.html_report(str(out), config_file=str(config))
        assert out.read_text(encoding="utf-8").count('class="plotly-graph-div"') == 1

    def test_a_failing_plot_does_not_lose_the_others(self, reporter, tmp_path,
                                                     caplog, monkeypatch):
        def broken(**kwargs):
            raise ValueError("unexpected")
        monkeypatch.setattr(reporter.plotter, "read_len_1D", broken)
        config = tmp_path / "c.json"
        config.write_text(json.dumps({"read_len_1D": {}, "run_summary": {}}))
        out = tmp_path / "r.html"
        with caplog.at_level(logging.WARNING):
            reporter.html_report(str(out), config_file=str(config))
        assert "read_len_1D failed with ValueError: unexpected" in caplog.text
        assert out.read_text(encoding="utf-8").count('class="plotly-graph-div"') == 1


class TestBarcodeSplit:

    def test_numeric_barcodes_with_blanks_keep_their_names(self, tmp_path):
        """A blank turned the column to float: sequencing_summary_6.0.txt."""
        barcodes = [6.0, 12.0, np.nan] * 40
        path = write_summary(tmp_path / "s.txt", n=120, barcode_arrangement=barcodes)
        out = tmp_path / "out"
        barcode_split(str(path), output_dir=str(out), output_unclassified=True, quiet=True)
        assert sorted(os.listdir(out)) == [
            "sequencing_summary_12.txt", "sequencing_summary_6.txt",
            "sequencing_summary_unclassified.txt"]

    def test_a_barcode_name_cannot_leave_the_output_directory(self, tmp_path):
        path = write_summary(tmp_path / "s.txt", n=100,
                             barcode_arrangement=["../escape"] * 100)
        out = tmp_path / "out"
        barcode_split(str(path), output_dir=str(out), quiet=True)
        assert os.listdir(out) == ["sequencing_summary_.._escape.txt"]
        assert not (tmp_path / "escape.txt").exists()


class TestPod5Duplicates:

    def test_a_read_in_two_files_is_reported(self, tmp_path, caplog):
        pytest.importorskip("pod5")
        d = tmp_path / "pod5"
        d.mkdir()
        _write_pod5(d / "a.pod5", n_reads=5)
        _write_pod5(d / "b.pod5", n_reads=5)
        with caplog.at_level(logging.WARNING):
            Pod5_to_seq_summary(pod5_dir=str(d), seq_summary_fn=str(tmp_path / "s.txt"),
                                threads=3)
        assert "more than one pod5 file" in caplog.text
        assert len(pd.read_csv(tmp_path / "s.txt", sep="\t")) == 5


class TestSpawnedWorkers:
    """Spawned workers, the default on macOS and from Python 3.14 on Linux,
    start with fresh loggers. Once the converters stopped configuring
    logging on import, their field counts went missing."""

    def test_field_counts_are_reported(self, fastq_dir, tmp_path):
        code = (
            "import multiprocessing as mp, sys\n"
            "from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary\n"
            "if __name__ == '__main__':\n"
            "    mp.set_start_method('spawn')\n"
            "    Fastq_to_seq_summary(fastq_dir=sys.argv[1], seq_summary_fn=sys.argv[2],\n"
            "                         threads=3, verbose_level=1)\n")
        script = tmp_path / "run.py"
        script.write_text(code)
        env = dict(os.environ, PYTHONPATH=str(REPO_ROOT))
        result = subprocess.run(
            [sys.executable, str(script), str(fastq_dir), str(tmp_path / "s.txt")],
            env=env, capture_output=True, text=True, timeout=300)
        assert result.returncode == 0, result.stderr
        assert "fields found" in result.stderr
        assert "fields not found" in result.stderr


class TestWorkersAreNotForked:
    """Forking a process that has a thread running can deadlock the child
    on a lock that thread held, and Python 3.12 warns about it. pytest, a
    notebook or an earlier progress bar is enough to leave such a thread."""

    def test_converters_use_a_fork_server_or_spawn(self):
        from pycoQC._pipeline import worker_context
        assert worker_context().get_start_method() in ("forkserver", "spawn")

    @pytest.mark.parametrize("converter", ["fastq", "fast5"])
    def test_no_fork_warning_with_a_thread_running(self, converter, fastq_dir,
                                                   basecalled_fast5_dir, tmp_path):
        call = {
            "fastq": "Fastq_to_seq_summary(fastq_dir=sys.argv[1], ",
            "fast5": "Fast5_to_seq_summary(fast5_dir=sys.argv[1], ",
        }[converter]
        code = (
            "import sys, threading\n"
            "from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary\n"
            "from pycoQC.Fast5_to_seq_summary import Fast5_to_seq_summary\n"
            "if __name__ == '__main__':\n"
            "    stop = threading.Event()\n"
            "    threading.Thread(target=stop.wait, daemon=True).start()\n"
            "    " + call + "seq_summary_fn=sys.argv[2], threads=3)\n"
            "    stop.set()\n")
        script = tmp_path / "run.py"
        script.write_text(code)
        in_dir = fastq_dir if converter == "fastq" else basecalled_fast5_dir
        out = tmp_path / "s.txt"
        env = dict(os.environ, PYTHONPATH=str(REPO_ROOT))
        result = subprocess.run(
            # Not -W error: os.fork swallows the exception that would make.
            [sys.executable, "-W", "always::DeprecationWarning", str(script),
             str(in_dir), str(out)],
            env=env, capture_output=True, text=True, timeout=300)
        assert result.returncode == 0, result.stderr
        assert "use of fork()" not in result.stderr
        assert len(pd.read_csv(out, sep="\t")) > 0


def _scripted_worker(conn, steps):
    """A stand-in worker: ("sleep", s), ("send", message), ("set", event),
    ("wait", event) or ("exit", code)."""
    import os
    import time
    for action, value in steps:
        if action == "sleep":
            time.sleep(value)
        elif action == "send":
            conn.send(value)
        elif action == "set":
            value.set()
        elif action == "wait":
            value.wait(60)
        else:
            conn.close()
            os._exit(value)


def _event():
    """An event the stand-in workers can share, to put their steps in order."""
    from pycoQC._pipeline import worker_context
    return worker_context().Event()


def _run_scripted(*scripts, grace=60, started=None):
    """Run stand-in workers, each reporting on a pipe of its own, as the
    pipeline does; returns wait_for_workers' result, or raises its error.
    The processes are added to `started`, when given, for the caller to
    check once they have been joined."""
    from pycoQC._pipeline import wait_for_workers, worker_context
    ctx = worker_context()
    workers = []
    try:
        for steps in scripts:
            reader, writer = ctx.Pipe(duplex=False)
            proc = ctx.Process(target=_scripted_worker, args=(writer, steps))
            proc.start()
            writer.close()
            workers.append((proc, reader))
        return wait_for_workers(workers, grace=grace)
    finally:
        for proc, reader in workers:
            proc.join(10)
            reader.close()
        if started is not None:
            started.extend(proc for proc, _ in workers)


class TestWhichErrorIsReported:
    """The first message to arrive decided the outcome, and messages from
    different processes arrive in no fixed order."""

    def test_the_cause_wins_over_its_consequence(self):
        """A lister finding no files leaves the writer with no reads, and
        whichever said so first was the error the user saw. The grace is
        long so that the order the two arrive in cannot matter."""
        from pycoQC._pipeline import LISTER, WRITER
        with pytest.raises(pycoQCError, match="no files found"):
            _run_scripted([("send", ("error", WRITER, "no reads to write"))],
                          [("sleep", 0.3), ("send", ("error", LISTER, "no files found"))])

    def test_an_error_after_the_report_still_fails_and_is_read(self):
        """The parent returned on the writer's report and then joined the
        workers, and one with a large message still unread could not exit."""
        from pycoQC._pipeline import PARSER
        started = []
        with pytest.raises(pycoQCError):
            _run_scripted([("send", ("done", {}))],
                          [("sleep", 0.3), ("send", ("error", PARSER, "x" * 300000))],
                          started=started)
        assert [p.exitcode for p in started] == [0, 0]

    def test_a_clean_run_returns_the_report(self):
        assert _run_scripted([("send", ("done", {"reads": 3}))]) == {"reads": 3}

    def test_a_crash_after_the_report_does_not_fail_the_run(self):
        """Every file had been read and the summary written, but a worker
        crashing on its way out failed the run and lost the summary. The
        crash waits for the report to be sent: with a fixed head start, a
        slow start of the first worker let the crash come first."""
        sent = _event()
        report = _run_scripted(
            [("send", ("done", {"reads": 3})), ("set", sent), ("exit", 0)],
            [("wait", sent), ("exit", 3)])
        assert report == {"reads": 3}

    def test_a_report_sent_just_before_a_crash_is_read(self):
        """An exit was judged before the pipe of the worker that made it had
        been read to its end."""
        report = _run_scripted([("send", ("done", {"reads": 3})), ("exit", 3)])
        assert report == {"reads": 3}

    def test_processes_never_started_can_be_stopped(self):
        """stop_processes joined every process in the list, and joining one
        that never started raises."""
        from pycoQC._pipeline import stop_processes, worker_context
        stop_processes([worker_context().Process(target=_scripted_worker, args=(None, []))])


class TestWindowsHandleLimit:
    """The parent waits on a pipe and an exit per worker, and before Python
    3.13 Windows waits on at most 63: --threads 32 failed at once. The cap
    first stopped at 3.12, which has the same limit."""

    @pytest.mark.parametrize("platform, version, expected", [
        ("win32", (3, 9), 31), ("win32", (3, 11), 31), ("win32", (3, 12), 31),
        ("win32", (3, 13), 40), ("linux", (3, 9), 40), ("darwin", (3, 12), 40)])
    def test_threads_are_capped_only_where_needed(self, platform, version, expected):
        from pycoQC._pipeline import usable_threads
        assert usable_threads(40, platform, version) == expected

    def test_fewer_threads_are_left_alone(self):
        from pycoQC._pipeline import usable_threads
        assert usable_threads(8, "win32", (3, 9)) == 8

    def test_too_few_threads_are_refused(self):
        from pycoQC._pipeline import usable_threads
        with pytest.raises(pycoQCError, match="At least 3 threads"):
            usable_threads(2)

    def test_the_cap_is_applied_and_said(self, fastq_dir, tmp_path, monkeypatch, caplog):
        import pycoQC._pipeline as pipeline
        monkeypatch.setattr(pipeline, "usable_threads", lambda threads: 3)
        conv = Fastq_to_seq_summary.__new__(Fastq_to_seq_summary)
        with caplog.at_level(logging.WARNING):
            conv._configure(str(fastq_dir), str(tmp_path / "s.txt"), 0, 40, 0, False, 0)
        assert conv.threads == 1
        assert "--threads lowered from 40 to 3" in caplog.text


@pytest.mark.skipif(sys.getfilesystemencoding().lower() != "utf-8",
                    reason="the byte counts are those of UTF-8")
class TestNameStart:
    """The staging directory kept 150 characters of the summary's name, and
    the limit is 255 bytes: 80 Chinese characters is 240."""

    def test_a_multibyte_name_is_cut_in_bytes(self):
        from pycoQC._pipeline import name_start
        start = name_start("测" * 80)
        assert start == "测" * 16
        assert len(os.fsencode(start)) <= 50

    def test_a_character_is_never_cut_in_half(self):
        from pycoQC._pipeline import name_start
        assert name_start("a" + "测" * 80) == "a" + "测" * 16

    def test_a_short_name_is_kept_whole(self):
        from pycoQC._pipeline import name_start
        assert name_start("summary.txt") == "summary.txt"

    @pytest.mark.skipif(os.name != "posix", reason="undecodable bytes are a POSIX matter")
    def test_an_undecodable_byte_is_kept(self):
        """Decoding with errors="ignore" dropped such a byte from the middle
        of the name, so a left over directory no longer matched its summary."""
        from pycoQC._pipeline import name_start
        name = os.fsdecode(b"caf\xe9.txt")
        assert name_start(name) == name


class TestSendError:

    def test_a_closed_pipe_does_not_hide_the_interrupt(self):
        """On Ctrl-C the parent closes its ends, and the BrokenPipeError of
        reporting replaced the KeyboardInterrupt that stopped the worker."""
        from multiprocessing import Pipe
        from pycoQC._pipeline import LISTER, _send_error
        reader, writer = Pipe(duplex=False)
        reader.close()
        try:
            with pytest.raises(KeyboardInterrupt):
                _send_error(writer, LISTER, KeyboardInterrupt())
        finally:
            writer.close()

    def test_an_exception_is_reported_and_not_raised(self):
        from multiprocessing import Pipe
        from pycoQC._pipeline import PARSER, _send_error
        reader, writer = Pipe(duplex=False)
        _send_error(writer, PARSER, pycoQCError("no reads"))
        assert reader.recv() == ("error", PARSER, "no reads")
        reader.close()
        writer.close()


    def test_an_exception_is_not_lost_with_the_pipe(self):
        """A failed send was swallowed, so a lister whose error could not be
        reported returned, exited cleanly, and its short summary passed for
        a complete one."""
        from multiprocessing import Pipe
        from pycoQC._pipeline import LISTER, _send_error
        reader, writer = Pipe(duplex=False)
        reader.close()
        try:
            with pytest.raises(pycoQCError, match="unreadable"):
                _send_error(writer, LISTER, pycoQCError("unreadable"))
        finally:
            writer.close()


class _FakeConn:
    """A pipe holding `messages`, which reports its end only once closed and
    fails where a real one would block."""

    def __init__(self, *messages):
        self.messages = list(messages)
        self.closed = False

    def recv(self):
        if self.messages:
            return self.messages.pop(0)
        if self.closed:
            raise EOFError
        raise AssertionError("recv on a drained, open pipe would block")


class _FakeProc:

    def __init__(self, conn):
        self.conn = conn
        self.sentinel = object()
        self.exitcode = 0

    def join(self, timeout=None):
        self.conn.closed = True


class TestDuplicateReadyHandles:

    def test_a_handle_listed_twice_is_read_once(self, monkeypatch):
        """On Windows, CPython 3.13's batched wait can list a ready handle
        twice, and the second recv on the drained pipe blocked."""
        import pycoQC._pipeline as pipeline
        conn = _FakeConn(("done", {"reads": 3}))
        proc = _FakeProc(conn)
        calls = iter([[conn, conn], [proc.sentinel, proc.sentinel]])
        monkeypatch.setattr(pipeline, "wait_ready", lambda objects, timeout: next(calls))
        assert pipeline.wait_for_workers([(proc, conn)]) == {"reads": 3}


class TestWorkersThatDieSilently:
    """A worker that died without reporting, failing as it started or
    killed for memory, left the parent waiting on the error queue for ever."""

    def test_the_parent_raises_instead_of_waiting(self):
        with pytest.raises(pycoQCError, match="exit code 3"):
            _run_scripted([("exit", 3)], grace=0.5)

    def test_a_script_without_a_main_guard_is_told_so(self, fastq_dir, tmp_path):
        script = tmp_path / "unguarded.py"
        script.write_text(
            "import sys\n"
            "from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary\n"
            "Fastq_to_seq_summary(fastq_dir=sys.argv[1], seq_summary_fn=sys.argv[2],\n"
            "                     threads=3)\n")
        env = dict(os.environ, PYTHONPATH=str(REPO_ROOT))
        result = subprocess.run(
            [sys.executable, str(script), str(fastq_dir), str(tmp_path / "s.txt")],
            env=env, capture_output=True, text=True, timeout=120)
        assert result.returncode != 0
        assert '__name__ == "__main__"' in result.stderr


@pytest.mark.skipif(not os.path.isdir("/proc/self/task"),
                    reason="counts threads through /proc")
class TestForkServerStaysSingleThreaded:
    """The fork server preloads pycoQC so workers share it. Importing numpy
    starts OpenBLAS threads and pandas 3 starts pyarrow's jemalloc thread,
    and a server holding either warns on every fork."""

    def test_preloading_starts_no_thread(self):
        code = ("import os, pycoQC._forkserver_preload\n"
                "print(len(os.listdir('/proc/self/task')))\n")
        env = dict(os.environ, PYTHONPATH=str(REPO_ROOT))
        result = subprocess.run([sys.executable, "-c", code], env=env,
                                capture_output=True, text=True, timeout=120)
        assert result.returncode == 0, result.stderr
        assert result.stdout.strip() == "1"



class _LateListerFailure(Fastq_to_seq_summary):
    """Lists every file, then fails: the writer finishes all the same."""

    def _work_items(self):
        yield from super()._work_items()
        raise pycoQCError("listing failed after the last file")


class _ForwardingFastq(Fastq_to_seq_summary):
    """A subclass whose __init__ names no fields parameter of its own."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


class _TwoFieldFastq(Fastq_to_seq_summary):
    """A subclass with a default of its own."""
    DEFAULT_FIELDS = ("read_id", "channel")


class _ExitingLister(Fastq_to_seq_summary):
    """Lists one file, then is stopped by something that is not an
    Exception. A class of its own per exception: a worker re-imports this
    module, so a class attribute set by a test would not reach it."""

    def _work_items(self):
        items = iter(super()._work_items())
        yield next(items)
        raise SystemExit(0)


class _InterruptedLister(_ExitingLister):

    def _work_items(self):
        items = iter(Fastq_to_seq_summary._work_items(self))
        yield next(items)
        raise KeyboardInterrupt


class _UniqueFastq(Fastq_to_seq_summary):
    """The FASTQ converter dropping repeated read ids, in blocks of 3."""
    unique_read_ids = True
    READS_PER_CHUNK = 3


def _report_env(q):
    import os
    q.put(os.environ.get("OPENBLAS_NUM_THREADS"))


def _write_bam(path, records):
    """An unaligned BAM holding (read_id, seq, qual) records, in order."""
    import pysam
    header = {"HD": {"VN": "1.6", "SO": "unknown"}}
    with pysam.AlignmentFile(str(path), "wb", header=header) as bam:
        for read_id, seq, qual in records:
            rec = pysam.AlignedSegment()
            rec.query_name = read_id
            rec.query_sequence = seq
            rec.query_qualities = pysam.qualitystring_to_array(qual)
            rec.flag = 4
            bam.write(rec)


class TestFailedConversionsWriteNothing:
    """The writer saved the summary before every worker had finished, so a
    run where one file failed left a summary missing that file's reads,
    which a later report could take for a complete one."""

    @pytest.fixture
    def broken_dir(self, tmp_path, bam_dir):
        import shutil
        d = tmp_path / "in"
        d.mkdir()
        shutil.copy(bam_dir / "reads.bam", d / "a.bam")
        (d / "b.bam").write_bytes(b"not a bam file at all")
        return d

    def test_no_summary_is_left_behind(self, broken_dir, tmp_path):
        out = tmp_path / "s.txt"
        with pytest.raises(pycoQCError):
            Bam_to_seq_summary(bam_dir=str(broken_dir), seq_summary_fn=str(out),
                               threads=4)
        assert not out.exists()
        assert sorted(p.name for p in tmp_path.iterdir()) == ["in"]

    def test_an_earlier_summary_is_kept(self, broken_dir, tmp_path):
        out = tmp_path / "s.txt"
        out.write_text("earlier summary\n")
        with pytest.raises(pycoQCError):
            Bam_to_seq_summary(bam_dir=str(broken_dir), seq_summary_fn=str(out),
                               threads=4)
        assert out.read_text() == "earlier summary\n"

    def test_a_file_of_the_users_named_like_the_temporary_one_is_kept(
            self, broken_dir, tmp_path):
        """The temporary file was <summary>.partial, and a failed run deleted
        whatever had that name, whether or not the run had written it."""
        mine = tmp_path / "s.txt.partial"
        mine.write_text("mine\n")
        with pytest.raises(pycoQCError):
            Bam_to_seq_summary(bam_dir=str(broken_dir),
                               seq_summary_fn=str(tmp_path / "s.txt"), threads=4)
        assert mine.read_text() == "mine\n"

    def test_an_error_after_the_writer_finished_still_fails(self, fastq_dir, tmp_path):
        """The writer's report could reach the parent before an error the
        lister sent first, and the summary was then published."""
        out = tmp_path / "s.txt"
        with pytest.raises(pycoQCError, match="listing failed after the last file"):
            _LateListerFailure(fastq_dir=str(fastq_dir), seq_summary_fn=str(out),
                               threads=3)
        assert not out.exists()

    @pytest.mark.parametrize("suffix", [".gz", ".bz2", ".xz", ".zip"])
    def test_a_compressed_summary_is_compressed(self, suffix, fastq_dir, tmp_path):
        """The summary is first written under a temporary name, which gave
        pandas no extension to infer compression from: a summary.txt.gz held
        plain text, and pycoQC could not read it."""
        out = tmp_path / ("s.txt" + suffix)
        Fastq_to_seq_summary(fastq_dir=str(fastq_dir), seq_summary_fn=str(out), threads=3)
        if suffix == ".gz":
            assert out.read_bytes()[:2] == b"\x1f\x8b"
        assert len(pd.read_csv(out, sep="\t")) == 200
        assert sorted(p.name for p in tmp_path.iterdir()) == [out.name]

    def test_a_long_summary_name_still_fits(self, fastq_dir, tmp_path):
        """The temporary directory was named after the whole summary name,
        and with its suffix passed the 255 byte limit of a file name."""
        out = tmp_path / ("s" * 240 + ".txt")
        Fastq_to_seq_summary(fastq_dir=str(fastq_dir), seq_summary_fn=str(out), threads=3)
        assert sorted(p.name for p in tmp_path.iterdir()) == [out.name]

    @pytest.mark.skipif(sys.getfilesystemencoding().lower() != "utf-8",
                        reason="the name needs a UTF-8 file system encoding")
    def test_a_long_multibyte_summary_name_still_fits(self, fastq_dir, tmp_path):
        out = tmp_path / ("测" * 80 + ".txt")
        Fastq_to_seq_summary(fastq_dir=str(fastq_dir), seq_summary_fn=str(out), threads=3)
        assert sorted(p.name for p in tmp_path.iterdir()) == [out.name]

    @pytest.mark.parametrize("lister, stop", [(_ExitingLister, "SystemExit"),
                                              (_InterruptedLister, "KeyboardInterrupt")])
    def test_a_lister_stopped_part_way_writes_nothing(self, lister, stop, fastq_dir,
                                                       tmp_path):
        """The lister caught only an Exception. Stopped by anything else, it
        still sent the parsers their stop signals, so the writer wrote and
        reported a summary missing the files not yet listed. Stopped by
        SystemExit(0), it even exited cleanly."""
        import shutil
        in_dir = tmp_path / "in"
        in_dir.mkdir()
        for name in ("a", "b"):
            shutil.copy(fastq_dir / "reads.fastq.gz", in_dir / (name + ".fastq.gz"))
        out = tmp_path / "s.txt"
        with pytest.raises(pycoQCError, match=stop):
            lister(fastq_dir=str(in_dir), seq_summary_fn=str(out), threads=3)
        assert not out.exists()

    def test_a_clean_run_leaves_only_the_summary(self, bam_dir, tmp_path):
        out = tmp_path / "s.txt"
        Bam_to_seq_summary(bam_dir=str(bam_dir), seq_summary_fn=str(out), threads=3)
        assert out.exists()
        assert sorted(p.name for p in tmp_path.iterdir()) == ["s.txt"]


class TestRepeatedReadIds:
    """Read ids repeated within one file were counted with those repeated
    across files, under a warning about the directory holding a run twice."""

    def test_repeats_across_files(self, reads, tmp_path, caplog):
        d = tmp_path / "in"
        d.mkdir()
        _write_bam(d / "a.bam", reads[:20])
        _write_bam(d / "b.bam", reads[:20])
        with caplog.at_level(logging.WARNING):
            Bam_to_seq_summary(bam_dir=str(d), seq_summary_fn=str(tmp_path / "s.txt"),
                               threads=3)
        assert "20 read id(s) appear in more than one bam file" in caplog.text
        assert "within one" not in caplog.text
        assert len(pd.read_csv(tmp_path / "s.txt", sep="\t")) == 20

    def test_repeats_within_a_split_file_are_within_one_file(self, reads, tmp_path,
                                                              caplog):
        """Repeats were told apart by work item, and the blocks of a split
        file are several work items: a read twice in one file read as a
        directory holding the run twice."""
        d = tmp_path / "in"
        d.mkdir()
        with gzip.open(d / "reads.fastq.gz", "wt") as fh:
            for read_id, seq, qual in reads[:10] + reads[:2]:
                fh.write(f"@{read_id} ch=1\n{seq}\n+\n{qual}\n")
        with caplog.at_level(logging.WARNING):
            # Six parsers and one file: it is split into blocks of 3 reads.
            _UniqueFastq(fastq_dir=str(d), seq_summary_fn=str(tmp_path / "s.txt"),
                         threads=8)
        assert "2 read id(s) appear more than once within one fastq file" in caplog.text
        assert "more than one fastq file" not in caplog.text

    def test_repeats_within_a_file(self, reads, tmp_path, caplog):
        d = tmp_path / "in"
        d.mkdir()
        _write_bam(d / "a.bam", reads[:20] + reads[:5])
        with caplog.at_level(logging.WARNING):
            Bam_to_seq_summary(bam_dir=str(d), seq_summary_fn=str(tmp_path / "s.txt"),
                               threads=3)
        assert "5 read id(s) appear more than once within one bam file" in caplog.text
        assert "more than one bam file" not in caplog.text
        assert len(pd.read_csv(tmp_path / "s.txt", sep="\t")) == 20


class TestMetricsFollowTheReadSets:
    """The figures were cached on id(frame), and an id is reused once the
    frame is freed, so a replaced read set could be given the old figures."""

    def test_a_replaced_read_set_is_not_kept_in_memory(self, guppy_summary):
        """The cache held each frame it measured, so a replaced read set
        stayed in memory for as long as the plotter."""
        import gc
        import weakref
        plotter = pycoQC_plot(pycoQC_parse(summary_file=str(guppy_summary), quiet=True),
                              quiet=True)
        old = plotter.pass_df
        ref = weakref.ref(old)
        plotter._metrics("pass")
        # A run this small is not sampled, so the sample is the frame too.
        plotter.pass_df = plotter.pass_sample_df = old.iloc[:10].copy()
        del old
        gc.collect()
        assert ref() is None

    def test_a_plotter_can_still_be_pickled_and_copied(self, guppy_summary):
        """The cache's weak references made a plotter that had drawn a table
        impossible to pickle or deep copy."""
        import copy
        import pickle
        plotter = pycoQC_plot(pycoQC_parse(summary_file=str(guppy_summary), quiet=True),
                              quiet=True)
        plotter.basecall_summary()
        for clone in (pickle.loads(pickle.dumps(plotter)), copy.deepcopy(plotter)):
            assert clone._metrics("all")["reads"] == len(plotter.all_df)

    def test_a_replaced_read_set_is_measured_again(self, guppy_summary):
        plotter = pycoQC_plot(pycoQC_parse(summary_file=str(guppy_summary), quiet=True),
                              quiet=True)
        assert plotter._metrics("pass")["reads"] == len(plotter.pass_df)
        plotter.pass_df = plotter.pass_df.iloc[:10]
        assert plotter._metrics("pass")["reads"] == 10


class TestRecordsCarryOnlyWhatIsAsked:

    @pytest.mark.parametrize("converter", ["fastq", "fast5"])
    def test_path_in_fields_is_taken_as_include_path(self, converter, fastq_dir,
                                                     basecalled_fast5_dir, tmp_path, caplog):
        """Only the BAM converter knew that path is not an input field: the
        FASTQ one counted it missing from every read and the FAST5 one
        rejected it as invalid."""
        out = tmp_path / "s.txt"
        with caplog.at_level(logging.WARNING):
            if converter == "fastq":
                Fastq_to_seq_summary(fastq_dir=str(fastq_dir), seq_summary_fn=str(out),
                                     threads=3, fields=["read_id", "channel", "path"])
            else:
                Fast5_to_seq_summary(fast5_dir=str(basecalled_fast5_dir),
                                     seq_summary_fn=str(out), threads=3,
                                     fields=["read_id", "channel", "path"])
        assert "only written with --include_path" in caplog.text
        assert "path" not in pd.read_csv(out, sep="\t").columns

    def test_asking_for_path_without_include_path_says_so(self, bam_dir, tmp_path, caplog):
        """The path is not in the BAM, so 'a BAM carries no field called path'
        blamed the file for an option that was not given."""
        with caplog.at_level(logging.WARNING):
            Bam_to_seq_summary(bam_dir=str(bam_dir), seq_summary_fn=str(tmp_path / "s.txt"),
                               threads=3, fields=["channel", "path", "path"])
        assert "only written with --include_path" in caplog.text
        assert "carries no field" not in caplog.text

    def test_a_field_no_bam_carries_is_named(self, bam_dir, tmp_path, caplog):
        with caplog.at_level(logging.WARNING):
            Bam_to_seq_summary(bam_dir=str(bam_dir), seq_summary_fn=str(tmp_path / "s.txt"),
                               threads=3, fields=["channel", "nonesuch"])
        assert ("a BAM carries no field called nonesuch, so it is left out of "
                "the summary") in caplog.text

    def test_a_field_named_twice_is_one_column(self, bam_dir, tmp_path):
        out = tmp_path / "s.txt"
        Bam_to_seq_summary(bam_dir=str(bam_dir), seq_summary_fn=str(out),
                           threads=3, fields=["channel", "channel"])
        assert list(pd.read_csv(out, sep="\t").columns) == ["read_id", "channel"]

    def test_a_field_named_twice_is_counted_once(self, fastq_dir, reads, tmp_path,
                                                 caplog):
        with caplog.at_level(logging.INFO):
            Fastq_to_seq_summary(fastq_dir=str(fastq_dir),
                                 seq_summary_fn=str(tmp_path / "s.txt"), threads=3,
                                 fields=["read_id", "channel", "channel"],
                                 verbose_level=2)
        text = caplog.text
        found = text[text.index("fields found"):text.index("fields not found")]
        assert "channel: {:,}".format(len(reads)) in found

    def test_path_stays_where_it_was_listed(self, bam_dir, tmp_path):
        out = tmp_path / "s.txt"
        Bam_to_seq_summary(bam_dir=str(bam_dir), seq_summary_fn=str(out), threads=3,
                           include_path=True, fields=["channel", "path", "mux"])
        columns = list(pd.read_csv(out, sep="\t").columns)
        assert columns[:3] == ["read_id", "channel", "path"]

    def test_an_empty_field_list_is_every_field(self, bam_dir, tmp_path):
        """fields=[] meant every field until the check became `is not None`."""
        every, empty = tmp_path / "every.txt", tmp_path / "empty.txt"
        Bam_to_seq_summary(bam_dir=str(bam_dir), seq_summary_fn=str(every), threads=3)
        Bam_to_seq_summary(bam_dir=str(bam_dir), seq_summary_fn=str(empty), threads=3,
                           fields=[])
        assert (list(pd.read_csv(empty, sep="\t").columns)
                == list(pd.read_csv(every, sep="\t").columns))

    def test_an_empty_field_list_is_the_fastq_default(self, fastq_dir, tmp_path):
        """BAM took an empty list for its default and FASTQ for no fields at
        all: every read came out empty, and the run failed on finding none."""
        every, empty = tmp_path / "every.txt", tmp_path / "empty.txt"
        Fastq_to_seq_summary(fastq_dir=str(fastq_dir), seq_summary_fn=str(every), threads=3)
        Fastq_to_seq_summary(fastq_dir=str(fastq_dir), seq_summary_fn=str(empty), threads=3,
                             fields=[])
        assert (list(pd.read_csv(empty, sep="\t").columns)
                == list(pd.read_csv(every, sep="\t").columns))

    def test_no_fields_is_the_default_for_every_converter(
            self, fastq_dir, basecalled_fast5_dir, tmp_path):
        """None was every field for BAM; FAST5 failed on it at once and
        FASTQ in every worker."""
        for name, run in [
                ("fastq", lambda fields, out: Fastq_to_seq_summary(
                    fastq_dir=str(fastq_dir), seq_summary_fn=str(out), threads=3,
                    fields=fields)),
                ("fast5", lambda fields, out: Fast5_to_seq_summary(
                    fast5_dir=str(basecalled_fast5_dir), seq_summary_fn=str(out),
                    threads=3, fields=fields))]:
            default, none = tmp_path / (name + "_default.txt"), tmp_path / (name + "_none.txt")
            run(Fastq_to_seq_summary.DEFAULT_FIELDS if name == "fastq"
                else Fast5_to_seq_summary.DEFAULT_FIELDS, default)
            run(None, none)
            assert (list(pd.read_csv(none, sep="\t").columns)
                    == list(pd.read_csv(default, sep="\t").columns))

    def test_an_empty_list_works_for_a_subclass(self, fastq_dir, tmp_path):
        """The default was looked up on the __init__ signature, and a
        subclass taking *args and **kwargs has no fields parameter there."""
        every, empty = tmp_path / "every.txt", tmp_path / "empty.txt"
        Fastq_to_seq_summary(fastq_dir=str(fastq_dir), seq_summary_fn=str(every), threads=3)
        _ForwardingFastq(fastq_dir=str(fastq_dir), seq_summary_fn=str(empty), threads=3,
                         fields=[])
        assert (list(pd.read_csv(empty, sep="\t").columns)
                == list(pd.read_csv(every, sep="\t").columns))

    def test_an_array_of_fields_is_a_field_list(self, fastq_dir, tmp_path):
        """`not fields` asked for an array's truth value, which raises."""
        out = tmp_path / "s.txt"
        Fastq_to_seq_summary(fastq_dir=str(fastq_dir), seq_summary_fn=str(out), threads=3,
                             fields=np.array(["read_id", "channel"]))
        assert list(pd.read_csv(out, sep="\t").columns) == ["read_id", "channel"]

    @pytest.mark.parametrize("fields", [None, []])
    def test_a_subclass_default_is_its_own(self, fields, fastq_dir, tmp_path):
        """The default was both the signature's and DEFAULT_FIELDS, and a
        subclass changing one got the other."""
        out = tmp_path / "s.txt"
        _TwoFieldFastq(fastq_dir=str(fastq_dir), seq_summary_fn=str(out), threads=3,
                       fields=fields)
        assert list(pd.read_csv(out, sep="\t").columns) == ["read_id", "channel"]

    @pytest.mark.parametrize("converter", [Fastq_to_seq_summary, Fast5_to_seq_summary,
                                           Bam_to_seq_summary])
    def test_default_fields_have_one_source(self, converter):
        import inspect
        assert inspect.signature(converter).parameters["fields"].default is None

    @pytest.mark.parametrize("converter", ["fastq", "bam"])
    def test_path_alone_without_include_path_is_the_read_ids(
            self, converter, fastq_dir, bam_dir, tmp_path):
        """Path is left out, which left BAM with an empty list it took for
        every field, and FASTQ with no fields at all."""
        out = tmp_path / "s.txt"
        if converter == "fastq":
            Fastq_to_seq_summary(fastq_dir=str(fastq_dir), seq_summary_fn=str(out),
                                 threads=3, fields=["path"])
        else:
            Bam_to_seq_summary(bam_dir=str(bam_dir), seq_summary_fn=str(out),
                               threads=3, fields=["path"])
        assert list(pd.read_csv(out, sep="\t").columns) == ["read_id"]

    def test_a_bam_record_has_no_path_unless_asked(self, bam_dir):
        from pycoQC.Bam_to_seq_summary import Bam_to_seq_summary as Bam
        from pycoQC._pipeline import _new_counts
        conv = Bam.__new__(Bam)
        conv.in_dir = str(bam_dir)
        for include_path in (False, True):
            conv.include_path = include_path
            record = next(conv._records(str(bam_dir / "reads.bam"), _new_counts()))
            assert ("path" in record) is include_path

    def test_blank_run_health_fields_count_as_not_found(self, fastq_dir):
        """A Guppy style header carries no Dorado tags, and the four run
        health fields were reported found, blank, for every read."""
        from pycoQC.Fastq_to_seq_summary import Fastq_to_seq_summary as Fq
        from pycoQC._pipeline import _new_counts
        conv = Fq.__new__(Fq)
        conv.include_path = False
        conv.reads_per_chunk = 1000
        conv.fields = ("read_id", "duration", "end_reason")
        counts = _new_counts()
        n = sum(1 for _ in conv._records((str(fastq_dir / "reads.fastq.gz"), "pass", None),
                                         counts))
        assert counts["fields_found"]["read_id"] == n
        assert counts["fields_not_found"]["duration"] == n
        assert counts["fields_not_found"]["end_reason"] == n
        assert counts["fields_found"]["duration"] == 0


class TestWorkerThreadCaps:
    """Spawned workers, on macOS and Windows, started with numpy's thread
    pools uncapped: idle threads per core in every worker."""

    @pytest.mark.parametrize("method", ["spawn", "worker_context"])
    def test_workers_start_capped_and_the_parent_is_restored(self, method, monkeypatch):
        import multiprocessing as mp
        from pycoQC._pipeline import _worker_environment, worker_context
        monkeypatch.delenv("OPENBLAS_NUM_THREADS", raising=False)
        ctx = mp.get_context("spawn") if method == "spawn" else worker_context()
        q = ctx.Queue()
        with _worker_environment():
            proc = ctx.Process(target=_report_env, args=(q,))
            proc.start()
        try:
            assert q.get(timeout=120) == "1"
        finally:
            proc.join(30)
        assert "OPENBLAS_NUM_THREADS" not in os.environ


class TestImportHasNoSideEffects:
    """Importing a converter configured the root logger at INFO, switching
    on every other library's chatter, set BLAS thread counts for the whole
    process, and pycoQC silenced FutureWarning everywhere."""

    def test_importing_the_package_leaves_the_process_alone(self):
        code = (
            "import logging, os, warnings\n"
            "before = list(warnings.filters)\n"
            "import pycoQC.__main__, pycoQC.pycoQC_report\n"
            "assert not logging.getLogger().handlers, logging.getLogger().handlers\n"
            "assert logging.getLogger().level == logging.WARNING\n"
            "assert 'OMP_NUM_THREADS' not in os.environ\n"
            "assert not any(f[2] is FutureWarning and f[0] == 'ignore' "
            "for f in warnings.filters if f not in before)\n")
        env = {k: v for k, v in os.environ.items() if not k.endswith("_NUM_THREADS")}
        env["PYTHONPATH"] = str(REPO_ROOT)
        result = subprocess.run([sys.executable, "-c", code], env=env,
                                capture_output=True, text=True)
        assert result.returncode == 0, result.stderr


class TestPdfRenderRetry:
    """kaleido now and then evaluates its script before the page holding it
    has loaded, failing with "kaleido_scopes is not defined". One failed
    page used to abort the whole PDF report."""

    @staticmethod
    def flaky(failures):
        calls = []

        def to_image(fig, **kwargs):
            calls.append(kwargs)
            if len(calls) <= failures:
                raise RuntimeError("ReferenceError: kaleido_scopes is not defined")
            return b"%PDF"
        return to_image, calls

    def test_a_transient_failure_is_retried(self, monkeypatch):
        import plotly.graph_objs as go
        from pycoQC.pycoQC_report import PDF_RENDER_ATTEMPTS, _render_pdf_page
        to_image, calls = self.flaky(PDF_RENDER_ATTEMPTS - 1)
        monkeypatch.setattr(go.Figure, "to_image", to_image)
        assert _render_pdf_page(go.Figure(), 300, 200) == b"%PDF"
        assert len(calls) == PDF_RENDER_ATTEMPTS
        assert calls[-1] == {"format": "pdf", "width": 300, "height": 200}

    def test_a_broken_engine_still_fails_the_report(self, monkeypatch, tmp_path,
                                                    guppy_summary):
        pytest.importorskip("pypdf")
        import plotly.graph_objs as go
        from pycoQC.pycoQC import pycoQC
        from pycoQC.pycoQC_report import PDF_RENDER_ATTEMPTS
        to_image, calls = self.flaky(10 ** 6)
        monkeypatch.setattr(go.Figure, "to_image", to_image)
        cfg = tmp_path / "config.json"
        cfg.write_text(json.dumps({"read_len_1D": {"plot_title": "Read length"}}))
        with pytest.raises(pycoQCError, match="kaleido_scopes"):
            pycoQC(summary_file=str(guppy_summary), config_file=str(cfg),
                   pdf_outfile=str(tmp_path / "r.pdf"), quiet=True)
        # The cover's failure is logged and skipped, then the first plot
        # exhausts its attempts and ends the report.
        assert len(calls) == 2 * PDF_RENDER_ATTEMPTS
