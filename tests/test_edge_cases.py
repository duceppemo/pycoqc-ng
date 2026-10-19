# -*- coding: utf-8 -*-
"""Degenerate inputs that still have to produce a report.

Real runs reach these: a clean run has no failing reads, a bad flowcell or
the wrong basecalling model leaves none passing, and a small test run has
only a handful of reads.
"""

import json
import random

import numpy as np
import pandas as pd
import pytest

from pycoQC.common import pycoQCError
from pycoQC.pycoQC import pycoQC
from pycoQC.pycoQC_plot import pycoQC_plot


def write_summary(path, n=800, qual=(2.0, 20.0), length=(200, 3000), seed=2,
                  **columns):
    rng = random.Random(seed)
    rows = []
    for i in range(n):
        row = {
            "read_id": f"r{i:05d}",
            "run_id": "run1",
            "channel": rng.randint(1, 512),
            "start_time": i * 2.0,
            "sequence_length_template": rng.randint(*length),
            "mean_qscore_template": round(rng.uniform(*qual), 2),
        }
        row.update({k: v(i, rng) if callable(v) else v
                    for k, v in columns.items()})
        rows.append(row)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
    return path


def build_report(tmp_path, summary, **kwargs):
    html = tmp_path / "report.html"
    js = tmp_path / "report.json"
    pycoQC(summary_file=str(summary), html_outfile=str(html),
           json_outfile=str(js), quiet=True, **kwargs)
    return html, json.loads(js.read_text())


class TestEmptyReadSets:

    def test_run_where_every_read_fails(self, tmp_path):
        """Regression: np.quantile raised IndexError on the empty pass set,
        so a bad run produced no report at all."""
        summary = write_summary(tmp_path / "fail.txt", qual=(2.0, 6.5))
        html, data = build_report(tmp_path, summary)
        assert html.exists()
        assert data["All Reads"]["basecall"]["reads_number"] == 800
        assert data["Pass Reads"]["basecall"]["reads_number"] == 0

    def test_run_where_every_read_passes(self, tmp_path):
        summary = write_summary(tmp_path / "pass.txt", qual=(12.0, 20.0))
        html, data = build_report(tmp_path, summary)
        assert html.exists()
        assert data["Pass Reads"]["basecall"]["reads_number"] == 800
        assert (data["All Reads"]["basecall"]["reads_number"]
                == data["Pass Reads"]["basecall"]["reads_number"])

    def test_empty_pass_set_reports_na_not_zero(self, tmp_path):
        """An absent median must not be reported as a real measurement."""
        summary = write_summary(tmp_path / "fail.txt", qual=(2.0, 6.5))
        _, data = build_report(tmp_path, summary)
        n50 = data["Pass Reads"]["basecall"]["N50"]
        assert n50 is None or (isinstance(n50, float) and np.isnan(n50))

    def test_threshold_exactly_on_the_boundary(self, tmp_path):
        """min_pass_qual is inclusive, so a run at exactly the threshold
        passes rather than falling in a gap."""
        summary = write_summary(tmp_path / "edge.txt", qual=(7.0, 7.0))
        _, data = build_report(tmp_path, summary, min_pass_qual=7)
        assert data["Pass Reads"]["basecall"]["reads_number"] == 800


class TestStatisticHelpers:

    @pytest.mark.parametrize("method", [
        "_compute_percentiles", "_compute_N50"])
    def test_helpers_accept_empty_input(self, method):
        result = getattr(pycoQC_plot, method)(pd.Series([], dtype=float))
        if isinstance(result, list):
            assert len(result) == 101 and all(np.isnan(v) for v in result)
        else:
            assert np.isnan(result)

    def test_percentiles_of_all_nan(self):
        values = pd.Series([np.nan, np.nan, np.nan])
        assert all(np.isnan(v) for v in pycoQC_plot._compute_percentiles(values))

    def test_histogram_of_empty_input(self):
        x, y = pycoQC_plot._compute_hist(pd.Series([], dtype=float))
        assert list(x) == [] and list(y) == []

    def test_median_or_nan_ignores_missing(self):
        assert pycoQC_plot._median_or_nan(pd.Series([1.0, np.nan, 3.0])) == 2.0
        assert np.isnan(pycoQC_plot._median_or_nan(pd.Series([], dtype=float)))
        assert np.isnan(
            pycoQC_plot._median_or_nan(pd.Series([np.nan, np.nan])))

    def test_single_value_percentiles_are_constant(self):
        result = pycoQC_plot._compute_percentiles(pd.Series([42.0]))
        assert len(result) == 101
        assert all(v == 42.0 for v in result)


class TestDegenerateDistributions:

    def test_every_read_the_same_length(self, tmp_path):
        """Fixed-length reads, as an untrimmed Illumina run gives."""
        summary = write_summary(tmp_path / "const.txt", length=(151, 151))
        _, data = build_report(tmp_path, summary)
        assert data["All Reads"]["basecall"]["N50"] == 151

    def test_zero_length_reads_are_filtered(self, tmp_path):
        summary = write_summary(tmp_path / "zero.txt", n=800)
        df = pd.read_csv(summary, sep="\t")
        df.loc[:99, "sequence_length_template"] = 0
        df.to_csv(summary, sep="\t", index=False)
        _, data = build_report(tmp_path, summary)
        assert data["All Reads"]["basecall"]["reads_number"] == 700

    def test_all_reads_in_one_channel(self, tmp_path):
        summary = write_summary(tmp_path / "chan.txt", channel=1)
        _, data = build_report(tmp_path, summary)
        assert data["All Reads"]["run"]["active_channels"] == 1

    def test_all_reads_at_the_same_time(self, tmp_path):
        """A zero-length run gives a zero time span, which must not divide."""
        summary = write_summary(tmp_path / "time.txt", start_time=0.0)
        html, data = build_report(tmp_path, summary)
        assert html.exists()
        assert data["All Reads"]["run"]["run_duration"] == 0


class TestFilteringRejectsEverything:

    def test_impossible_quality_threshold_is_reported(self, tmp_path):
        """Nothing can pass q60, and the failure should say so rather than
        raising something opaque."""
        summary = write_summary(tmp_path / "s.txt")
        html = tmp_path / "r.html"
        pycoQC(summary_file=str(summary), html_outfile=str(html),
               min_pass_qual=60, quiet=True)
        assert html.exists()

    def test_impossible_length_threshold_is_reported(self, tmp_path):
        summary = write_summary(tmp_path / "s.txt")
        html = tmp_path / "r.html"
        pycoQC(summary_file=str(summary), html_outfile=str(html),
               min_pass_len=10_000_000, quiet=True)
        assert html.exists()

    def test_unknown_runid_selection_is_rejected(self, tmp_path):
        summary = write_summary(tmp_path / "s.txt")
        with pytest.raises(pycoQCError):
            pycoQC(summary_file=str(summary),
                   html_outfile=str(tmp_path / "r.html"),
                   runid_list=["no_such_run"], quiet=True)
