# -*- coding: utf-8 -*-
"""Tests for the statistics helpers used to build the plots."""

import numpy as np
import pandas as pd
import pytest

from pycoQC.pycoQC_plot import pycoQC_plot


def reference_n50(values):
    """N50 by its usual definition: the largest length for which the reads
    at least that long hold half the bases."""
    values = sorted(values, reverse=True)
    half = sum(values) / 2
    total = 0
    for v in values:
        total += v
        if total >= half:
            return int(v)
    return int(values[-1]) if values else None


class TestComputeN50:

    def test_known_value(self):
        # 1+2+3+4+5 = 15, half = 7.5, cumulative reaches it at 4.
        assert pycoQC_plot._compute_N50(pd.Series([1, 2, 3, 4, 5])) == 4

    def test_single_read(self):
        assert pycoQC_plot._compute_N50(pd.Series([1234])) == 1234

    def test_ignores_missing_values(self):
        with_nan = pd.Series([1, 2, np.nan, 3, 4, 5])
        assert (pycoQC_plot._compute_N50(with_nan)
                == pycoQC_plot._compute_N50(pd.Series([1, 2, 3, 4, 5])))

    def test_empty_input_is_nan(self):
        assert np.isnan(pycoQC_plot._compute_N50(pd.Series([], dtype=float)))

    def test_all_zero_lengths(self):
        assert pycoQC_plot._compute_N50(pd.Series([0.0, 0.0, 0.0])) == 0

    def test_matches_reference_on_random_input(self):
        rng = np.random.default_rng(0)
        for _ in range(200):
            n = int(rng.integers(1, 500))
            values = rng.integers(1, 100_000, size=n).astype(float)
            assert (pycoQC_plot._compute_N50(pd.Series(values))
                    == reference_n50(values))

    def test_accepts_a_read_only_series(self):
        """pandas may hand out immutable arrays; sorting must not be in place."""
        values = np.arange(1, 100, dtype=float)
        values.flags.writeable = False
        series = pd.Series(values, copy=False)
        assert pycoQC_plot._compute_N50(series) == reference_n50(values)


class TestComputePercentiles:

    def test_returns_101_points_from_min_to_max(self):
        result = pycoQC_plot._compute_percentiles(pd.Series(np.arange(101.0)))
        assert len(result) == 101
        assert result[0] == pytest.approx(0.0)
        assert result[-1] == pytest.approx(100.0)
        assert result[50] == pytest.approx(50.0)


class TestN50Ties:
    """Where a prefix sums to exactly half the bases, the read above the
    boundary is the answer."""

    @pytest.mark.parametrize("values, expected", [
        ([1, 2, 3, 4, 10], 10),
        ([2, 2, 2, 2, 4, 4], 4),
        ([5, 5], 5),
        ([0, 0, 0], 0),
        ([7], 7),
    ])
    def test_known_values(self, values, expected):
        assert pycoQC_plot._compute_N50(pd.Series(values, dtype=float)) == expected
