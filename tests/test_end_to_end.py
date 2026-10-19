# -*- coding: utf-8 -*-
"""End-to-end report generation, through the API and the command line."""

import json
import subprocess
import sys

import pytest

from pycoQC.pycoQC import pycoQC


@pytest.fixture(scope="module")
def reports(tmp_path_factory, guppy_summary):
    """Build the HTML and JSON reports once for the whole module."""
    out = tmp_path_factory.mktemp("report")
    html = out / "report.html"
    js = out / "report.json"
    pycoQC(summary_file=str(guppy_summary),
           html_outfile=str(html), json_outfile=str(js), quiet=True)
    return html, js


class TestHtmlAndJsonReport:

    def test_html_report_is_written(self, reports):
        html, _ = reports
        assert html.exists()
        # A report carries every plot inline, so it is never a stub.
        assert html.stat().st_size > 100_000
        assert "<html" in html.read_text(errors="replace").lower()

    def test_json_report_has_expected_structure(self, reports):
        _, js = reports
        data = json.loads(js.read_text())
        assert "pycoqc" in data
        assert "All Reads" in data
        basecall = data["All Reads"]["basecall"]
        assert basecall["reads_number"] > 0
        assert basecall["N50"] > 0
        assert basecall["bases_number"] > 0

    def test_pass_reads_are_a_subset_of_all_reads(self, reports):
        _, js = reports
        data = json.loads(js.read_text())
        assert (data["Pass Reads"]["basecall"]["reads_number"]
                <= data["All Reads"]["basecall"]["reads_number"])


class TestCommandLine:

    def test_version_flag(self):
        out = subprocess.run(
            [sys.executable, "-m", "pycoQC", "--version"],
            capture_output=True, text=True)
        assert out.returncode == 0
        assert "pycoqc-ng" in (out.stdout + out.stderr)

    def test_default_config_is_valid_json(self):
        """--default_config reads a packaged resource, which breaks if the
        package data is not installed alongside the code."""
        out = subprocess.run(
            [sys.executable, "-m", "pycoQC", "--default_config"],
            capture_output=True, text=True)
        assert out.returncode == 0
        assert isinstance(json.loads(out.stdout), dict)

    def test_generates_a_report_from_the_cli(self, tmp_path, guppy_summary):
        html = tmp_path / "cli_report.html"
        out = subprocess.run(
            [sys.executable, "-m", "pycoQC", "-f", str(guppy_summary),
             "-o", str(html), "--quiet"],
            capture_output=True, text=True)
        assert out.returncode == 0, out.stderr
        assert html.exists()


class TestPackageResources:

    def test_config_template_is_packaged(self):
        from importlib.resources import files
        config = files("pycoQC").joinpath("templates/pycoQC_config.json")
        assert isinstance(json.loads(config.read_text()), dict)

    def test_jinja_template_is_packaged(self):
        import jinja2
        env = jinja2.Environment(
            loader=jinja2.PackageLoader("pycoQC", "templates"),
            autoescape=jinja2.select_autoescape(["html"]))
        assert env.get_template("spectre.html.j2") is not None
