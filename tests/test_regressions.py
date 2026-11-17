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
from pycoQC.Barcode_split import barcode_split
from pycoQC.common import (
    check_arg, expand_file_names, format_worker_error, head, pycoQCError)
from pycoQC.Fast5_to_seq_summary import Fast5_to_seq_summary
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
