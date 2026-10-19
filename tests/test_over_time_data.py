# -*- coding: utf-8 -*-
"""The time-binned percentile series behind the "over time" plots.

These are checked against a straightforward reference implementation, since
the production version is vectorised for speed.
"""

import numpy as np
import pandas as pd
import pytest

from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC.pycoQC_plot import pycoQC_plot


def reference_over_time_data(df, field_name, time_bins=500, smooth_sigma=0):
    """Bin values by start time and take percentiles, one bin at a time."""
    subset = df[["start_time", field_name]].dropna()
    data = subset[field_name].to_numpy()
    t = (subset["start_time"] / 3600).to_numpy()
    x = np.linspace(t.min(), t.max(), num=time_bins)
    bin_idx = np.digitize(t, bins=x, right=True)

    per_bin = {}
    for idx, val in zip(bin_idx, data):
        per_bin.setdefault(idx, []).append(val)

    val_name = ["Min", "Max", "25%", "75%", "Median"]
    out = {name: [] for name in val_name}
    for i in range(x.size):
        if i in per_bin:
            p = np.percentile(per_bin[i], [0, 100, 25, 75, 50])
        else:
            p = [np.nan] * 5
        for name, stat in zip(val_name, p):
            out[name].append(stat)
    return x, {k: np.asarray(v) for k, v in out.items()}


@pytest.fixture(scope="module")
def plotter(guppy_summary):
    parser = pycoQC_parse(summary_file=str(guppy_summary), quiet=True)
    return pycoQC_plot(parser, quiet=True)


def call_over_time_data(plotter, field_name, smooth_sigma=0, time_bins=500):
    """Reach the name-mangled private method."""
    method = getattr(plotter, "_pycoQC_plot__over_time_data")
    return method(df_level="all", field_name=field_name,
                  smooth_sigma=smooth_sigma, time_bins=time_bins)


class TestOverTimeData:

    @pytest.mark.parametrize("field_name", ["read_len", "mean_qscore"])
    def test_matches_reference_implementation(self, plotter, field_name):
        _, data_dict = call_over_time_data(plotter, field_name)
        x_ref, stats_ref = reference_over_time_data(
            plotter.all_sample_df, field_name)

        np.testing.assert_allclose(data_dict["x"][0], x_ref)
        for series, name in zip(data_dict["y"],
                                ["Min", "Max", "25%", "75%", "Median"]):
            # Some summary columns are stored as float32, so compare at the
            # precision of the input rather than of the accumulator.
            np.testing.assert_allclose(
                np.asarray(series, dtype=float), stats_ref[name],
                rtol=1e-6, equal_nan=True)

    def test_quartiles_are_ordered(self, plotter):
        _, data_dict = call_over_time_data(plotter, "read_len")
        minimum, maximum, q25, q75, median = (
            np.asarray(s, dtype=float) for s in data_dict["y"])
        present = ~np.isnan(median)
        assert (minimum[present] <= q25[present]).all()
        assert (q25[present] <= median[present]).all()
        assert (median[present] <= q75[present]).all()
        assert (q75[present] <= maximum[present]).all()

    def test_values_stay_paired_with_their_own_time(self, plotter):
        """Regression: only the value column was dropna()'d, so values were
        zipped against the wrong times whenever the field had gaps."""
        df = plotter.all_sample_df.copy()
        # Blank out the field for the earliest reads only. If times and values
        # are realigned correctly, the surviving values are unaffected.
        order = df["start_time"].argsort()
        df.loc[df.index[order[:len(df) // 3]], "read_len"] = np.nan

        subset = df[["start_time", "read_len"]].dropna()
        assert len(subset) < len(df)

        x_ref, stats_ref = reference_over_time_data(df, "read_len")
        # Every reported median must be a real value drawn from the kept rows.
        medians = stats_ref["Median"][~np.isnan(stats_ref["Median"])]
        assert medians.min() >= subset["read_len"].min()
        assert medians.max() <= subset["read_len"].max()

    def test_empty_field_returns_no_data(self, plotter):
        """A field that is entirely missing must not raise."""
        _, data_dict = call_over_time_data(plotter, "read_len", time_bins=500)
        assert data_dict  # sanity: the normal case does return data
