# -*- coding: utf-8 -*-
"""The console entry points, driven in process.

Each main_* function takes an argument list. It used to accept one and then
call parse_args() with no argument, so the list was ignored and the real
command line was parsed instead.
"""

import json

import pytest

from pycoQC.__main__ import (
    main_Bam_to_seq_summary,
    main_Barcode_split,
    main_Fast5_to_seq_summary,
    main_Fastq_to_seq_summary,
    main_Pod5_to_seq_summary,
    main_pycoQC,
)

ENTRY_POINTS = [
    main_pycoQC,
    main_Fast5_to_seq_summary,
    main_Fastq_to_seq_summary,
    main_Pod5_to_seq_summary,
    main_Bam_to_seq_summary,
    main_Barcode_split,
]


class TestArgumentHandling:

    @pytest.mark.parametrize("entry", ENTRY_POINTS, ids=lambda f: f.__name__)
    def test_help_exits_cleanly(self, entry, capsys):
        with pytest.raises(SystemExit) as excinfo:
            entry(["--help"])
        assert excinfo.value.code == 0
        assert "usage:" in capsys.readouterr().out

    @pytest.mark.parametrize("entry", ENTRY_POINTS, ids=lambda f: f.__name__)
    def test_version_reports_the_package(self, entry, capsys):
        with pytest.raises(SystemExit) as excinfo:
            entry(["--version"])
        assert excinfo.value.code == 0
        assert "pycoqc-ng" in capsys.readouterr().out

    @pytest.mark.parametrize("entry", ENTRY_POINTS, ids=lambda f: f.__name__)
    def test_passed_args_are_used_not_sys_argv(self, entry, monkeypatch, capsys):
        """Regression: the args parameter was accepted and then ignored, so a
        caller's list had no effect and sys.argv was parsed instead."""
        monkeypatch.setattr("sys.argv", ["unused", "--definitely-not-a-flag"])
        with pytest.raises(SystemExit) as excinfo:
            entry(["--help"])
        # Reaching a clean exit means --help was honoured; parsing sys.argv
        # would have failed on the unknown flag with code 2.
        assert excinfo.value.code == 0
        assert "usage:" in capsys.readouterr().out

    def test_unknown_flag_is_rejected(self):
        with pytest.raises(SystemExit) as excinfo:
            main_pycoQC(["--no-such-option"])
        assert excinfo.value.code == 2

    def test_default_config_is_valid_json(self, capsys):
        with pytest.raises(SystemExit) as excinfo:
            main_pycoQC(["--default_config"])
        assert excinfo.value.code == 0
        assert isinstance(json.loads(capsys.readouterr().out), dict)

    def test_missing_summary_file_exits_nonzero(self):
        """Regression: usage errors called sys.exit() with no status, which
        the interpreter reports as 0, so a workflow manager treated a missing
        required argument as a successful run."""
        with pytest.raises(SystemExit) as excinfo:
            main_pycoQC(["-o", "unused.html"])
        assert excinfo.value.code == 2

    def test_missing_output_file_exits_nonzero(self, guppy_summary):
        with pytest.raises(SystemExit) as excinfo:
            main_pycoQC(["-f", str(guppy_summary)])
        assert excinfo.value.code == 2


class TestReportThroughCli:

    def test_generates_html_and_json(self, tmp_path, guppy_summary):
        html = tmp_path / "report.html"
        js = tmp_path / "report.json"
        main_pycoQC(["-f", str(guppy_summary), "-o", str(html),
                     "-j", str(js), "--quiet"])
        assert html.exists() and html.stat().st_size > 100_000
        data = json.loads(js.read_text())
        assert data["All Reads"]["basecall"]["reads_number"] > 0

    def test_report_title_is_applied(self, tmp_path, guppy_summary):
        html = tmp_path / "titled.html"
        main_pycoQC(["-f", str(guppy_summary), "-o", str(html),
                     "--report_title", "A Distinctive Title", "--quiet"])
        assert "A Distinctive Title" in html.read_text(errors="replace")


class TestConverterCli:

    def test_fastq_converter_via_cli(self, tmp_path, fastq_dir):
        out = tmp_path / "fq.txt"
        main_Fastq_to_seq_summary(["-f", str(fastq_dir), "-s", str(out),
                                   "-t", "3"])
        assert out.exists()
        assert len(out.read_text().splitlines()) > 1

    def test_bam_converter_via_cli(self, tmp_path, bam_dir):
        out = tmp_path / "bam.txt"
        main_Bam_to_seq_summary(["-f", str(bam_dir), "-s", str(out), "-t", "3"])
        assert out.exists()
        assert len(out.read_text().splitlines()) > 1
