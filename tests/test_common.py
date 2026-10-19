# -*- coding: utf-8 -*-
"""Tests for the shared helpers in pycoQC.common."""

import math

import numpy as np
import pandas as pd
import pytest

from pycoQC.common import (mean_qscore, mkdir, parse_timestamp, pycoQCError,
                           start_time_to_elapsed_seconds)
from conftest import reference_mean_qscore


class TestMeanQscore:

    def test_matches_reference_for_every_read(self, reads):
        """The table lookup must agree with a direct error-probability average."""
        for _, _, qual in reads:
            assert mean_qscore(qual, phred33=True) == pytest.approx(
                reference_mean_qscore(qual), abs=1e-9)

    def test_accepts_bytes_and_str_alike(self):
        qual = "IIIII$$###"
        assert mean_qscore(qual, phred33=True) == pytest.approx(
            mean_qscore(qual.encode("ascii"), phred33=True))

    def test_raw_phred_values_match_phred33_string(self):
        """Raw scores and the Phred+33 encoding of them must agree."""
        scores = [2, 7, 14, 30, 40]
        encoded = "".join(chr(q + 33) for q in scores)
        assert mean_qscore(bytes(scores)) == pytest.approx(
            mean_qscore(encoded, phred33=True))

    def test_uniform_quality_returns_that_quality(self):
        """Averaging equal scores cannot drift away from the score itself."""
        for q in (2, 10, 20, 40):
            assert mean_qscore(bytes([q] * 50)) == pytest.approx(q)

    def test_dominated_by_low_quality_bases(self):
        """A log-scale average must sit below the arithmetic mean."""
        scores = [40] * 99 + [2]
        got = mean_qscore(bytes(scores))
        assert got < np.mean(scores)

    @pytest.mark.parametrize("empty", ["", b"", None, []])
    def test_missing_qualities_are_zero(self, empty):
        """pysam hands out None when a record carries no qualities."""
        assert mean_qscore(empty) == 0.0


class TestStartTimeToElapsedSeconds:

    def test_counts_seconds_from_the_first_read(self):
        times = pd.Series([
            "2022-06-27T23:00:00Z",
            "2022-06-27T23:00:30Z",
            "2022-06-27T23:02:00Z"])
        assert list(start_time_to_elapsed_seconds(times)) == [0, 30, 120]

    def test_unordered_input_is_relative_to_the_minimum(self):
        times = pd.Series([
            "2022-06-27T23:02:00Z",
            "2022-06-27T23:00:00Z",
            "2022-06-27T23:00:30Z"])
        assert list(start_time_to_elapsed_seconds(times)) == [120, 0, 30]

    def test_integer_dtype_when_nothing_is_missing(self):
        times = pd.Series(["2022-06-27T23:00:00Z", "2022-06-27T23:00:30Z"])
        assert start_time_to_elapsed_seconds(times).dtype == np.dtype("int64")

    def test_missing_times_fall_back_to_float(self):
        """A NaT cannot be cast to int, so the column stays floating point."""
        times = pd.Series(["2022-06-27T23:00:00Z", None])
        result = start_time_to_elapsed_seconds(times)
        assert result.dtype == np.dtype("float64")
        assert result.iloc[0] == 0.0
        assert math.isnan(result.iloc[1])


class TestParseTimestamp:

    @pytest.mark.parametrize("value", [
        "2022-06-27T23:14:28Z",
        "2021-03-29T15:18:10.066+00:00",
        "2022-06-27T23:14:28",
        "2022-06-27T23:14:28.123456Z",
        "2021-03-29T15:18:10.066-04:00",
    ])
    def test_agrees_with_dateutil(self, value):
        """The fast path must not change how timestamps are interpreted."""
        from dateutil.parser import parse as dateutil_parse
        assert parse_timestamp(value) == dateutil_parse(value)

    def test_utc_designator_is_an_offset_not_a_naive_time(self):
        parsed = parse_timestamp("2022-06-27T23:14:28Z")
        assert parsed.utcoffset().total_seconds() == 0

    def test_falls_back_for_non_iso_formats(self):
        """Unusual basecaller output still has to parse."""
        from dateutil.parser import parse as dateutil_parse
        value = "27 June 2022 23:14:28"
        assert parse_timestamp(value) == dateutil_parse(value)

    def test_ordering_is_preserved(self):
        earlier = parse_timestamp("2022-06-27T23:14:28Z")
        later = parse_timestamp("2022-06-27T23:15:28Z")
        assert (later - earlier).total_seconds() == 60


class TestMkdir:

    def test_unwritable_path_raises_pycoqc_error(self):
        """Regression: this raised NameError, because the handler referred to
        pycoMethError, a class from a different project."""
        with pytest.raises(pycoQCError) as excinfo:
            mkdir("/proc/this/cannot/be/created", exist_ok=True)
        assert "Error creating output folder" in str(excinfo.value)

    def test_creates_nested_directories(self, tmp_path):
        target = tmp_path / "a" / "b" / "c"
        mkdir(str(target), exist_ok=True)
        assert target.is_dir()
