# -*- coding: utf-8 -*-
"""The read sets every plot draws from, the run's properties, and the helpers the plot families share.

Part of pycoQC_plot, which composes the plot families in pycoQC.plots into
one class. Methods here reach the read sets, the helpers and the other
families through self.
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library imports
import re

# Third party imports
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d

# Local lib import
from pycoQC.common import get_logger, pycoQCError
from pycoQC.pycoQC_parse import pycoQC_parse

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~GLOBAL SETTINGS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Above this many barcodes the read counts are drawn as bars rather than
# a pie, which stops being legible well before a 96 barcode kit.
PIE_MAX_BARCODES = 20

# Seed for the deterministic sampling of large runs.
SEED = 42

# The read sets every interactive plot can switch between, in button order.
LEVELS = ("all", "pass", "fail")

# Past this many references the coverage plot stops labelling and shading
# each one. Nobody can read a hundred labels along one axis, and a
# transcriptome or a fragmented assembly has hundreds of thousands, each of
# which cost a tick label and a rectangle in the report.
COVERAGE_MAX_LABELLED_REFS = 100


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class PlotBase:
    def __init__(self,
                 parser: pycoQC_parse,
                 min_pass_qual: int = 7,
                 min_pass_len: int = 0,
                 sample: int = 100000,
                 verbose: bool = False,
                 quiet: bool = False):
        """
        * parser
            A pycoQC_parse object
        * min_pass_qual
            Minimum quality to consider a read as 'pass'
        * min_pass_len
            Minimum read length to consider a read as 'pass'
        * sample
            If not None a n number of reads will be randomly selected instead of the entire dataset for plotting function(deterministic sampling)
        """

        # Set logging level
        self.logger = get_logger(name=__name__, verbose=verbose, quiet=quiet)
        self.logger.info("Loading plotting interface")

        # Save args to self values
        self.min_pass_qual = min_pass_qual
        self.sample = sample

        # Check that parser is a valid instance of pycoQC_parse
        if not isinstance(parser, pycoQC_parse):
            raise pycoQCError("{} is not a valid pycoQC_parse object".format(parser))
        self.parser = parser

        # Extract values from parser object
        self.all_df = self._with_speed(parser.reads_df)
        if self.has_alignment:
            self.ref_len_dict = parser.ref_len_dict
            self.alignments_df = parser.alignments_df
        self.logger.info("\tFound {:,} total reads".format(len(self.all_df)))

        # Pass and fail. Selected by mask rather than by dropping the
        # passing index labels: read ids are not guaranteed unique, and
        # dropping a repeated label takes every row that carries it, which
        # emptied the fail set on a run whose files repeat ids. Without
        # basecalls there is no quality score to threshold on, so the
        # basecaller's own verdict is used where it recorded one: MinKNOW
        # sorts its files into pass and fail folders, and the POD5 and FAST5
        # converters carry that through as passes_filtering.
        passing = self._passing_mask(min_pass_qual, min_pass_len)
        self.pass_df = self.all_df[passing]
        self.fail_df = self.all_df[~passing]

        # Each read set is sampled on its own, so a small fail set is kept
        # whole rather than thinned along with the run.
        self.all_sample_df, self.all_scaling_factor = self._sample(self.all_df, sample)
        self.pass_sample_df, self.pass_scaling_factor = self._sample(self.pass_df, sample)
        self.fail_sample_df, self.fail_scaling_factor = self._sample(self.fail_df, sample)

        if self.has_basecalls:
            self.logger.info("\tFound {:,} pass reads (qual >= {} and length >= {})".format(
                len(self.pass_df), min_pass_qual, min_pass_len))
        else:
            self.logger.info("\tFound {:,} pass reads".format(len(self.pass_df)))
        self.logger.info("\tFound {:,} fail reads".format(len(self.fail_df)))

    def _with_speed(self, df):
        """Add each read's translocation speed, in bases per second.

        The time is the strand's own where the summary records it, as Guppy
        and Dorado do in template_duration and the converters derive from
        Dorado's ns and ts tags. The whole read's duration also counts the
        adapter, which a short read spends a large share of its time on:
        on a real Dorado run, reads under 500 bases came out at 298 bases a
        second against 392 from the template time, while every length
        agreed at 392 to 394 from the template time. The whole duration is
        used only when nothing better is recorded.
        """
        self.speed_from_template = "template_duration" in df
        if "read_len" not in df:
            return df
        if self.speed_from_template:
            time = df["template_duration"]
            if "duration" in df:
                time = time.fillna(df["duration"])
        elif "duration" in df:
            time = df["duration"]
        else:
            return df
        time = pd.to_numeric(time, errors="coerce").astype(np.float64)
        speed = df["read_len"].astype(np.float64) / time.where(time > 0)
        return df.assign(speed=speed.astype(np.float32))

    @staticmethod
    def _sample(df, sample):
        """A deterministic sample of at most `sample` reads, and how many
        reads each sampled one stands for."""
        if sample and len(df) > sample:
            return df.sample(n=sample, random_state=SEED), len(df) / sample
        return df, 1

    def _passing_mask(self, min_pass_qual, min_pass_len):
        """Which reads count as passing, however this run records that."""
        if self.has_basecalls:
            return ((self.all_df["mean_qscore"] >= min_pass_qual)
                    & (self.all_df["read_len"] >= min_pass_len))

        if "passes_filtering" in self.all_df:
            flag = self.all_df["passes_filtering"]
            if flag.dtype == bool:
                return flag
            # A blank verdict means nothing judged this read, so it has not
            # failed; counting it as a failure invented a fail count.
            text = flag.astype(str).str.strip().str.lower()
            passed = text.isin(("true", "1", "yes", "pass"))
            return passed | flag.isna() | text.isin(("", "nan", "none"))

        # Nothing has judged these reads, so none of them has failed.
        return pd.Series(True, index=self.all_df.index)

    def __str__(self):
        m = ""
        m += "\tBarcode: {}\n".format(self.has_barcodes)
        m += "\tAlignment: {}\n".format(self.has_alignment)
        m += "\tPromethion: {}\n".format(self.is_promethion)
        m += "\tFlongle: {}\n".format(self.is_flongle)

        # A run that has not been basecalled has no read_len, and this
        # raised KeyError from inside the logging call that printed it, so
        # `pycoQC -v` on a POD5 summary dumped a traceback mid run.
        for label, df in (("All", self.all_df), ("Pass", self.pass_df), ("Fail", self.fail_df)):
            m += "\t{} reads: {:,}\n".format(label, len(df))
            if self.has_basecalls:
                m += "\t{} bases: {:,}\n".format(label, int(df["read_len"].sum()))
                m += "\t{} median read length: {}\n".format(
                    label, self._median_or_nan(df["read_len"]))
            if self.has_gc_percent:
                m += "\t{} median %GC: {}\n".format(
                    label, self._median_or_nan(df["%GC"]))

        return m

    def __repr__(self):
        return "[{}]\n".format(self.__class__.__name__)

    def __getstate__(self):
        # The figures cache holds weak references, which cannot be pickled
        # or deep copied, and is rebuilt on demand anyway.
        state = self.__dict__.copy()
        state.pop("_metrics_cache", None)
        return state

    @property
    def has_barcodes(self):
        return "barcode" in self.all_df

    @property
    def has_basecalls(self):
        """Whether anything has been basecalled.

        A POD5 file, or a FAST5 holding signal only, records which pore
        produced which read and when, but no sequence: there is no read
        length and no quality score until a basecaller has run. The run is
        still worth describing, so the plots that need a sequence say so and
        are skipped rather than the whole report failing.
        """
        return "read_len" in self.all_df and "mean_qscore" in self.all_df

    @property
    def has_duration(self):
        return "duration" in self.all_df

    def _require_duration(self):
        if not self.has_duration:
            raise pycoQCError("No read duration information available")

    def _require_basecalls(self):
        if not self.has_basecalls:
            raise pycoQCError(
                "No basecall information available. This summary holds raw "
                "signal only: basecall the run to get read lengths and "
                "quality scores.")

    @property
    def has_speed(self):
        return "speed" in self.all_df and self.all_df["speed"].notna().any()

    @property
    def has_end_reason(self):
        return "end_reason" in self.all_df

    @property
    def has_mux(self):
        return "mux" in self.all_df

    @property
    def _speed_label(self):
        return ("Bases per second" if self.speed_from_template
                else "Bases per second (whole read, with adapter)")

    @property
    def _flow_cell_channels(self):
        """How many channels the flow cell has, judged from the highest one
        seen. A channel numbered past the platform's count, from a device or
        a summary this does not know, widens it."""
        if self.is_promethion:
            n_channels = 3000
        elif self.is_flongle:
            n_channels = 126
        else:
            n_channels = 512  # MinION / GridION
        return max(n_channels, int(self.all_df["channel"].max()))

    @property
    def has_gc_percent(self):
        return "%GC" in self.all_df

    @property
    def has_alignment(self):
        return "ref_id" in self.all_df

    @property
    def has_identity_freq(self):
        return "identity_freq" in self.all_df

    @property
    def is_promethion(self):
        if "channel" not in self.all_df:
            return False
        return self.all_df["channel"].max() > 512

    @property
    def is_flongle(self):
        if "channel" not in self.all_df:
            return False
        return self.all_df["channel"].max() < 127

    @property
    def total_ref_len(self):
        if self.has_alignment:
            return np.sum(list(self.ref_len_dict.values()))

    @staticmethod
    def _barcode_order(labels):
        """Barcode names in the order a reader expects them.

        Barcodes are names, not numbers, so plain sorting put a Deepbinner
        run's barcode 10 between 1 and 2. Names are compared with their
        digits as numbers, so barcode2 comes before barcode10 too, which
        plain text sorting got wrong for every kit not padded to two digits.
        Bare numbers come first and unclassified last.
        """
        def key(label):
            text = str(label)
            if text == "unclassified":
                return (2, 0.0, [])
            try:
                return (0, float(text), [])
            except ValueError:
                # Alternating text and digit runs, starting with text, so
                # each position always compares like with like.
                chunks = re.split(r"(\d+)", text)
                return (1, 0.0, [int(c) if i % 2 else c for i, c in enumerate(chunks)])

        return sorted(labels, key=key)

    def _level_dfs(self, df_level):
        """The whole frame and the sample drawn from it, for one read set.

        Counting reads and summing bases is a single pass over a column, so
        those come from the whole run and are exact. Medians, N50s and the
        densities stay on the sample, which is what sampling is for.
        """
        if df_level == "pass":
            return self.pass_df, self.pass_sample_df
        if df_level == "fail":
            return self.fail_df, self.fail_sample_df
        return self.all_df, self.all_sample_df

    def _level_sample(self, df_level):
        """The sample for one read set, and how many reads each sampled
        read stands for."""
        if df_level == "pass":
            return self.pass_sample_df, self.pass_scaling_factor
        if df_level == "fail":
            return self.fail_sample_df, self.fail_scaling_factor
        return self.all_sample_df, self.all_scaling_factor

    @staticmethod
    def _level_label(df_level, unit="Reads"):
        return "{} {}".format(df_level.capitalize(), unit)

    @staticmethod
    def _level_menu(prepared, method, **position):
        """The row of buttons switching a plot between read sets.

        * prepared
            (label, *args) per button, args being what plotly's method takes:
            a trace update, and for "update" a layout update after it.
        * method
            "restyle" or "update".
        * position
            Placement of the buttons, passed to the updatemenu.
        """
        return [dict(type="buttons", active=0,
                     buttons=[dict(label=label, method=method, args=list(args))
                              for label, *args in prepared],
                     **position)]

    @staticmethod
    def _median_or_nan(values):
        """Median of a column, or NA when there is nothing to take it of

        A read set can legitimately be empty: a run where every read passes
        has no failing reads. np.median warns and returns nan for an empty
        input, so the emptiness is checked rather than relied upon.
        """
        values = values.dropna()
        if values.empty:
            return np.nan
        return np.median(values)

    @staticmethod
    def _time_bins(t, time_bins):
        """Divide the span of `t` into time_bins equal intervals.

        Returns the centre of each interval and the interval each value
        falls in. Every over time plot bins this way, so their x axes line
        up with one another and with the channel heatmap's y axis.
        """
        edges = np.linspace(t.min(), t.max(), num=time_bins + 1)
        centres = (edges[:-1] + edges[1:]) / 2
        idx = np.clip(np.digitize(t, bins=edges[1:-1]), 0, time_bins - 1)
        return centres, idx

    def _require_speed(self):
        if not self.has_speed:
            raise pycoQCError(
                "No translocation speed available: it needs a read length "
                "and a duration for each read")

    @staticmethod
    def _compute_percentiles(data):
        # A read set can legitimately be empty, in either direction: a clean
        # run has no failing reads, and a bad flowcell or the wrong
        # basecalling model can leave none passing. np.quantile raises on an
        # empty input, so report NA for each percentile instead.
        values = data.dropna()
        if values.empty:
            return [np.nan] * 101
        return list(np.quantile(values, q=np.linspace(0, 1, 101)))

    @staticmethod
    def _compute_N50(data):
        # np.sort returns a new array: pandas hands out read-only views, so an
        # in-place sort would fail. The cumulative sum replaces a Python loop.
        values = np.sort(data.dropna().to_numpy())
        if values.size == 0:
            return np.nan
        cum_sum = np.cumsum(values)
        # N50 is the largest length for which the reads at least that long
        # hold half the bases. When a prefix sums to exactly half, "left"
        # picked the read below the boundary: [1,2,3,4,10] gave 4 rather
        # than 10.
        # All-zero lengths put the target at the end of the cumulative sum,
        # so the index has to be held inside the array.
        index = np.searchsorted(cum_sum, cum_sum[-1] / 2, side="right")
        return int(values[min(index, values.size - 1)])

    @staticmethod
    def _hist_bins(min_value, max_value, nbins, scale="linear"):
        """Bin edges spanning the data, with the largest value inside them.

        Extending the top edge by 0.1 in log10 made room for it, but that is
        26% of the way to the next decade: the length density ran out to
        127 kb on a run whose longest read was 101 kb. The edge only has to
        clear the rounding error from the log10 round trip, which a relative
        nudge does. Equal values still need a width to be drawn at all.
        """
        if scale == "log":
            # A log axis cannot place zero, and log10(0) is -inf, which made
            # every edge NaN and emptied the histogram. Read lengths are
            # filtered for this; durations are not, and a zero sample read
            # is a real thing to find in a pod5.
            if min_value <= 0:
                min_value = max_value * 1e-6 if max_value > 0 else 1.0
                max_value = max(max_value, min_value)
            if max_value <= min_value:
                return np.logspace(np.log10(min_value * 0.95),
                                   np.log10(max_value * 1.05), nbins)
            # Both ends need the nudge. Going through log10 and back moves
            # the edge either way, and the bottom one moved up, which put
            # the shortest read in the run outside the histogram entirely.
            return np.logspace(np.log10(min_value * (1 - 1e-9)),
                               np.log10(max_value * (1 + 1e-9)), nbins)

        if max_value <= min_value:
            return np.linspace(min_value - 0.5, max_value + 0.5, nbins)
        return np.linspace(min_value, max_value, nbins)

    @staticmethod
    def _compute_hist(data, x_scale="linear", smooth_sigma=2, nbins=200):

        # Nothing to bin when the read set is empty, or when none of it is a
        # number. A BAM aligned against the wrong reference leaves every
        # align_len missing, and nanmin of that is itself nan: the bin edges
        # all became nan and numpy then counted the missing values into the
        # top bins, so the JSON reported a histogram of reads that did not
        # align to anything.
        data = np.asarray(data, dtype=np.float64).ravel()
        data = data[np.isfinite(data)]
        if data.size == 0:
            return [], []

        # Count each categories in log or linear space
        min_value = data.min()
        max_value = data.max()

        if x_scale == "log":
            count_y, bins = np.histogram(
                a=data, bins=PlotBase._hist_bins(min_value, max_value, nbins, "log"))
        elif x_scale == "linear":
            count_y, bins = np.histogram(
                a=data, bins=PlotBase._hist_bins(min_value, max_value, nbins))
        else:  # Should not get here
            count_y = ''
            bins = ''

        # Label each bin by its centre. Using the right edge put every
        # count half a bin too far along the axis, which is visible on a
        # log scale where the first bins are wide: a 200 base read was
        # drawn at the position of a 260 base one.
        count_x = (bins[:-1] + bins[1:]) / 2

        # Smooth results with a gaussian filter. np.histogram returns
        # integer counts and gaussian_filter1d keeps the input dtype, so
        # smoothing truncated every bin towards zero: a six read set came
        # out as all zeros, and the real run lost 1,217 of 975,715 reads
        # from its length histogram.
        if smooth_sigma:
            count_y = gaussian_filter1d(count_y.astype(np.float64),
                                        sigma=smooth_sigma)

        # Convert to python list
        count_x = [float(i) for i in count_x]
        count_y = [float(i) for i in count_y]

        return count_x, count_y
