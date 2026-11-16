# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library imports
from collections import OrderedDict, namedtuple
import datetime
import re

# Third party imports
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, gaussian_filter1d
import plotly.graph_objs as go
from plotly.subplots import make_subplots

# Local lib import
from pycoQC.common import get_logger, pycoQCError
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC import __version__ as package_version

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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class pycoQC_plot():
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

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PROPERTY METHODS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
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

    # Run metadata is optional. A summary built from a plain fastq, or from an
    # Illumina run, carries no timing, channel or run id, so these report NA
    # rather than failing and taking the whole summary table with them.
    def _run_duration(self, df):
        if "start_time" not in df:
            return np.nan
        times = df["start_time"].dropna()
        if times.empty:
            return np.nan
        return float(times.max() - times.min()) / 3600

    def _active_channels(self, df):
        if "channel" not in df:
            return np.nan
        return int(df["channel"].nunique())

    def _runid_number(self, df):
        if "run_id" not in df:
            return np.nan
        return int(df["run_id"].nunique())

    def _barcodes_number(self, df):
        n_bc = 0
        if self.has_barcodes:
            n_bc = int(df["barcode"].nunique())
            if 'unclassified' in df["barcode"].values:
                n_bc -= 1
        return n_bc

    def _barcodes_count(self, df):
        counts = df["barcode"].value_counts()
        counts = counts.reindex(self._barcode_order(counts.index))
        return counts.to_frame()

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

    def _basecalled_reads(self, df):
        return len(df)

    def _basecalled_bases(self, df):
        return int(df["read_len"].sum())

    def _basecall_N50(self, df):
        return self._compute_N50(df["read_len"])

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

    def _basecall_median_read_len(self, df):
        return self._median_or_nan(df["read_len"])

    def _basecall_median_read_qscore(self, df):
        return self._median_or_nan(df["mean_qscore"])

    def _basecall_median_read_gc(self, df):
        return self._median_or_nan(df["%GC"])

    def _alignment_mean_coverage(self, df):
        """Mean depth over the reference.

        Measured as the reference each alignment spans, which is what the
        coverage plot integrates. Summing align_len instead counts aligned
        query bases, which insertions inflate and deletions deflate, and the
        report then gave two different mean coverages for one run: 13.77X in
        this table against 13.37X on the plot.
        """
        if not self.has_alignment:
            return np.nan
        if "ref_start" not in df or "ref_end" not in df:
            return np.nan
        spans = (df["ref_end"] - df["ref_start"]).dropna()
        return spans.sum() / self.total_ref_len

    def _aligned_reads(self, df):
        return len(df["align_len"].dropna()) if self.has_alignment else np.nan

    def _aligned_bases(self, df):
        return int(df["align_len"].dropna().sum()) if self.has_alignment else np.nan

    def _alignment_N50(self, df):
        return self._compute_N50(df["align_len"]) if self.has_alignment else np.nan

    def _alignment_median_read_len(self, df):
        return self._median_or_nan(df["align_len"]) if self.has_alignment else np.nan

    def _alignment_median_identity(self, df):
        return self._median_or_nan(df["identity_freq"]) if self.has_identity_freq else np.nan

    def _alignment_insertion_rate(self, df):
        return df["insertion"].dropna().sum()/self._aligned_bases(df) if self.has_identity_freq else np.nan

    def _alignment_deletion_rate(self, df):
        return df["deletion"].dropna().sum() / self._aligned_bases(df) if self.has_identity_freq else np.nan

    def _alignment_mismatch_rate(self, df):
        return df["mismatch"].dropna().sum() / self._aligned_bases(df) if self.has_identity_freq else np.nan

    #~~~~~~~SUMMARY_STATS_DICT METHOD AND HELPER~~~~~~~#

    def summary_stats_dict(self):
        """
        Return a dictionnary containing exhaustive information about the run.
        """
        self.logger.info("\tCompute overall summary statistics")
        d = OrderedDict()

        d["pycoqc"] = OrderedDict()
        d["pycoqc"]["version"] = package_version
        d["pycoqc"]["date"] = datetime.datetime.now().strftime("%d/%m/%y")

        for df, lab in ((self.all_df, "All Reads"), (self.pass_df, "Pass Reads")):
            d[lab] = self._compute_stats(df)
        return d

    def _compute_stats(self, df):
        d = OrderedDict()
        # run information
        d["run"] = OrderedDict()
        d["run"]["run_duration"] = self._run_duration(df)
        d["run"]["active_channels"] = self._active_channels(df)
        d["run"]["runid_number"] = self._runid_number(df)
        d["run"]["barcodes_number"] = self._barcodes_number(df)
        if self.has_end_reason:
            counts = df["end_reason"].fillna("unknown").value_counts()
            d["run"]["end_reasons"] = OrderedDict(
                (str(k), int(v)) for k, v in counts.items())
        # Everything under "basecall" describes a sequence, which a raw
        # signal run does not have yet. What it does have is how long each
        # read occupied its pore, which is reported in its place.
        if self.has_basecalls:
            d["basecall"] = OrderedDict()
            d["basecall"]["reads_number"] = self._basecalled_reads(df)
            d["basecall"]["bases_number"] = self._basecalled_bases(df)
            d["basecall"]["N50"] = self._basecall_N50(df)
            d["basecall"]["len_percentiles"] = self._compute_percentiles(df["read_len"])
            d["basecall"]["qual_score_percentiles"] = self._compute_percentiles(df["mean_qscore"])

            x, y = self._compute_hist(data=df["read_len"], x_scale="log", smooth_sigma=2, nbins=100)
            d["basecall"]["len_hist"] = OrderedDict()
            d["basecall"]["len_hist"]["x"] = x
            d["basecall"]["len_hist"]["y"] = y
            x, y = self._compute_hist(data=df["mean_qscore"], x_scale="linear", smooth_sigma=2, nbins=100)
            d["basecall"]["qual_score_hist"] = OrderedDict()
            d["basecall"]["qual_score_hist"]["x"] = x
            d["basecall"]["qual_score_hist"]["y"] = y
            if self.has_speed:
                d["basecall"]["median_speed"] = self._median_or_nan(df["speed"])
                d["basecall"]["speed_percentiles"] = self._compute_percentiles(df["speed"])
        else:
            d["signal"] = OrderedDict()
            d["signal"]["reads_number"] = len(df)
            if self.has_duration:
                d["signal"]["median_duration"] = self._median_or_nan(df["duration"])
                d["signal"]["total_duration"] = float(df["duration"].sum())
                d["signal"]["duration_percentiles"] = self._compute_percentiles(df["duration"])
                x, y = self._compute_hist(data=df["duration"], x_scale="log",
                                          smooth_sigma=2, nbins=100)
                d["signal"]["duration_hist"] = OrderedDict()
                d["signal"]["duration_hist"]["x"] = x
                d["signal"]["duration_hist"]["y"] = y

        if self.has_alignment:
            d["alignment"] = OrderedDict()
            d["alignment"]["reads_number"] = self._aligned_reads(df)
            d["alignment"]["bases_number"] = self._aligned_bases(df)
            d["alignment"]["mean_coverage"] = self._alignment_mean_coverage(df)
            d["alignment"]["N50"] = self._alignment_N50(df)
            d["alignment"]["len_percentiles"] = self._compute_percentiles(df["align_len"])
            x, y = self._compute_hist(data=df["align_len"], x_scale="log", smooth_sigma=2, nbins=100)
            d["alignment"]["len_hist"] = OrderedDict()
            d["alignment"]["len_hist"]["x"] = x
            d["alignment"]["len_hist"]["y"] = y

            if self.has_identity_freq:
                d["alignment"]["identity_freq_percentiles"] = self._compute_percentiles(df["identity_freq"])
                d["alignment"]["insertion_rate"] = self._alignment_insertion_rate(df)
                d["alignment"]["deletion_rate"] = self._alignment_deletion_rate(df)
                d["alignment"]["mismatch_rate"] = self._alignment_mismatch_rate(df)
                x, y = self._compute_hist(data=df["identity_freq"], x_scale="linear", smooth_sigma=2, nbins=100)
                d["alignment"]["identity_freq_hist"] = OrderedDict()
                d["alignment"]["identity_freq_hist"]["x"] = x
                d["alignment"]["identity_freq_hist"]["y"] = y

        return d

    #~~~~~~~SUMMARY METHODS AND HELPER~~~~~~~#

    def run_summary(self,
                    width: int = None,
                    height: int = 350,
                    plot_title: str = "General run summary"):
        """
        Plot an interactive overall summary table
        * groupby
            Value of field to group the data in the table
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Extract data
        data = []
        for status, df in (("All Reads", self.all_df), ("Pass Reads", self.pass_df), ("Fail Reads", self.fail_df)):
            data.append([
                status,
                self._run_duration(df),
                self._active_channels(df),
                self._runid_number(df),
                self._barcodes_number(df)])

        fig = self._summary_plot(
            width=width,
            height=height,
            plot_title=plot_title,
            header=["Status", "Run Duration(h)", "Active Channels", "Number of RunIDs", "Number of Barcodes"],
            data_format=["", ".2f", "", "", ""],
            data=[*zip(*data)])

        return fig

    def basecall_summary(self,
                         width: int = None,
                         height: int = 350,
                         plot_title: str = "Basecall summary"):
        """
        Plot an interactive basecall summary table
        * groupby
            Value of field to group the data in the table
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self._require_basecalls()
        # Extract data
        data = []
        for status, df in (("All Reads", self.all_df), ("Pass Reads", self.pass_df), ("Fail Reads", self.fail_df)):
            row = [
                status,
                self._basecalled_reads(df),
                self._basecalled_bases(df),
                self._basecall_N50(df),
                self._basecall_median_read_len(df),
                self._basecall_median_read_qscore(df)]
            if self.has_gc_percent:
                row.append(self._basecall_median_read_gc(df))
            data.append(row)

        fig_header_list = ["Status", "Reads", "Bases", "N50", "Median Length", "Median Q-score"]
        data_format_list = ["", ",", ",", ",", ",", ".1f"]
        if '%GC' in self.all_df:
            fig_header_list.append("Median %GC")
            data_format_list.append(".1f")

        fig = self._summary_plot(
            width=width,
            height=height,
            plot_title=plot_title,
            header=fig_header_list,
            data_format=data_format_list,
            data=[*zip(*data)])

        return fig

    def alignment_summary(self,
                          width: int = None,
                          height: int = 350,
                          plot_title: str = "Alignment summary"):
        """
        Plot an interactive alignment summary table
        * groupby
            Value of field to group the data in the table
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignment information are available
        if not self.has_alignment:
            raise pycoQCError("No alignment information available")

        data = []
        for status, df in (("All Reads", self.all_df), ("Pass Reads", self.pass_df), ("Fail Reads", self.fail_df)):
            data.append([
                status,
                self._aligned_reads(df),
                self._aligned_bases(df),
                self._alignment_mean_coverage(df),
                self._alignment_N50(df),
                self._alignment_median_read_len(df),
                self._alignment_median_identity(df)])

        fig = self._summary_plot(
            width=width,
            height=height,
            plot_title=plot_title,
            header=["Status", "Reads", "Bases", "Mean Coverage", "N50", "Median Read Length", "Median Identity Freq"],
            # d3's "r" rounds to significant digits, so ".1r" turned a median
            # read length of 4,732 into 5,000. Identity is a fraction between
            # 0 and 1, which ".1f" rendered as a flat 1.0.
            data_format=["", ",", ",", ".1f", ",", ",", ".4f"],
            data=[*zip(*data)])

        return fig

    def _summary_plot(self, width, height, plot_title, header, data_format, data):
        """Private function generating summary table plots"""
        self.logger.info("\t\tComputing plot")

        # Plot data
        data = [go.Table(
            header={
                "values": header,
                "align": "center", "fill": {"color": "grey"},
                "font": {"size": 14, "color": "white"},
                "height": 40},
            cells={
                "values": data,
                "format": data_format,
                "align": "center",
                "fill": {"color": "whitesmoke"},
                "font": {"size": 12}, "height": 30})]

        # tweak plot layout
        layout = go.Layout(
            width=width,
            height=height,
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"})

        return go.Figure(data=data, layout=layout)

    def basecall_summary_per_barcode(self,
                                     width: int = None,
                                     height: int = None,
                                     plot_title: str = "Basecall summary per barcode"):
        """
        Plot an interactive basecall summary table
        * groupby
            Value of field to group the data in the table
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel. Grows with the number of
            barcodes unless set.
        * plot_title
            Title to display on top of the plot
        """
        self._require_basecalls()

        # Verify that barcode information are available
        if not self.has_barcodes:
            raise pycoQCError("No barcode information available")

        self.logger.info("\t\tComputing plot")

        fig_header_list = ["Barcode", "Reads", "Bases", "N50", "Median Length", "Median Q-score"]
        data_format_list = ["", ",", ",", ",", ",", ".1f"]
        if '%GC' in self.all_df:
            fig_header_list.append("Median %GC")
            data_format_list.append(".1f")

        fig = self._summary_barcode_plot(
            width=width,
            height=height,
            plot_title=plot_title,
            header=fig_header_list,
            data_format=data_format_list)

        return fig

    def _summary_barcode_plot(self, width, height, plot_title, header, data_format):
        """Private function generating summary table plot per barcode"""
        self.logger.info("\t\tComputing plot")

        # Prepare all data
        prepared = [self._summary_barcode_data(level) for level in LEVELS]
        _, dd1 = prepared[0]

        # One row per barcode, so a 96 barcode run needs a taller table than
        # a 12 barcode one rather than being cut off at the same place.
        if height is None:
            height = max(500, 290 + 30 * len(dd1["x"][0]))

        def cells(dd):
            return {
                "values": [*zip(*dd['y'][0])],
                "format": data_format,
                "align": "center",
                "fill": {"color": "whitesmoke"},
                "font": {"size": 12}, "height": 30}

        # Plot initial data
        data = [go.Table(
            header={
                "values": header,
                "align": "center", "fill": {"color": "grey"},
                "font": {"size": 14, "color": "white"},
                "height": 40},
            cells=cells(dd1))]

        # tweak plot layout
        layout = go.Layout(
            width=width,
            height=height,
            updatemenus=self._level_menu(
                [(label, dict(cells=cells(dd))) for label, dd in prepared], "restyle",
                x=-0.2, y=0.5, xanchor='left', yanchor='bottom'),
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"})

        return go.Figure(data=data, layout=layout)

    def _summary_barcode_data(self, df_level):

        # Reads and bases are counted over the whole run, so this table and
        # the basecall summary on the same page agree to the read. Scaling a
        # sample up was unbiased but noisy, and because pass and fail are
        # sampled separately their barcode counts did not add up to the all
        # reads row. The distribution statistics still come from the sample.
        df, sample_df = self._level_dfs(df_level)

        # Extract data
        data = []
        bc_list = self._barcode_order(df['barcode'].unique())

        for bc in bc_list:
            # Extract data for barcode
            df1 = df[df['barcode'] == bc]
            sample1 = sample_df[sample_df['barcode'] == bc]
            # A barcode with only a handful of reads can miss the sample
            # entirely, and its row then showed a real read count beside a
            # blank N50 and median. Such a barcode is small enough to
            # measure in full.
            if sample1.empty:
                sample1 = df1
            row = [
                bc,
                self._basecalled_reads(df1),
                self._basecalled_bases(df1),
                # Medians and N50 are distribution statistics, so the sample
                # estimates them directly.
                self._basecall_N50(sample1),
                self._basecall_median_read_len(sample1),
                self._basecall_median_read_qscore(sample1)]
            if self.has_gc_percent:
                row.append(self._basecall_median_read_gc(sample1))
            data.append(row)

        data_dict = dict(
            x=[bc_list],
            y=[data]
        )

        return self._level_label(df_level), data_dict

    def signal_summary(self,
                       width: int = None,
                       height: int = 350,
                       plot_title: str = "Signal summary"):
        """Summarise a run that has not been basecalled.

        What a POD5 or a raw FAST5 holds is the signal each pore produced
        and when, so this counts reads and the time they spent in a pore.
        There is deliberately no base count: bases do not exist until a
        basecaller has read the signal, and an estimate from the signal
        length would be a guess dressed up as a measurement.

        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info("\t\tComputing plot")

        header = ["Status", "Reads", "Median Duration (s)", "Pore Time (h)",
                  "Active Channels"]
        data_format = ["", ",", ".2f", ",.1f", ","]

        rows = []
        for df, status in ((self.all_df, "All Reads"),
                           (self.pass_df, "Pass Reads"),
                           (self.fail_df, "Fail Reads")):
            rows.append([
                status,
                len(df),
                self._median_or_nan(df["duration"]) if self.has_duration else np.nan,
                float(df["duration"].sum()) / 3600 if self.has_duration else np.nan,
                self._active_channels(df)])

        # The table takes its values column by column.
        return self._summary_plot(
            data=[*zip(*rows)], width=width, height=height,
            plot_title=plot_title, header=header, data_format=data_format)

    def read_duration_1D(self,
                         color: str = "mediumaquamarine",
                         nbins: int = 200,
                         smooth_sigma: int = 2,
                         width: int = None,
                         height: int = 500,
                         plot_title: str = "Read duration"):
        """
        Plot a distribution of how long each read occupied its pore
        * color
            Color of the area(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x-axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self._require_duration()
        return self._1D_density_plot(
            field_name="duration",
            plot_title=plot_title,
            x_lab="Time in pore (s)",
            color=color,
            x_scale="log",
            nbins=nbins,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height)

    def read_duration_over_time(
            self,
            median_color: str = "rgb(102,168,255)",
            quartile_color: str = "rgb(153,197,255)",
            extreme_color: str = "rgba(153,197,255,0.5)",
            smooth_sigma: int = 1,
            time_bins: int = 500,
            width: int = None,
            height: int = 550,
            plot_title: str = "Read duration over experiment time"):
        """
        Plot read duration over experiment time
        * median_color
            Color of median line color(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in(x-axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self._require_duration()
        return self._over_time_plot(
            field_name="duration",
            plot_title=plot_title,
            y_lab="Time in pore (s)",
            y_scale="log",
            median_color=median_color,
            quartile_color=quartile_color,
            extreme_color=extreme_color,
            smooth_sigma=smooth_sigma,
            time_bins=time_bins,
            width=width,
            height=height)

    #~~~~~~~1D DISTRIBUTION METHODS AND HELPER~~~~~~~#
    def read_len_1D(self,
                    color: str = "lightsteelblue",
                    nbins: int = 200,
                    smooth_sigma: int = 2,
                    width: int = None,
                    height: int = 500,
                    plot_title: str = "Basecalled reads length"):
        """
        Plot a distribution of read length(log scale)
        * color
            Color of the area(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x-axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self._require_basecalls()
        fig = self._1D_density_plot(
            field_name="read_len",
            plot_title=plot_title,
            x_lab="Basecalled length",
            color=color,
            x_scale="log",
            nbins=nbins,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height)
        return fig

    def read_qual_1D(self,
                     color: str = "salmon",
                     nbins: int = 200,
                     smooth_sigma: int = 2,
                     width: int = None,
                     height: int = 500,
                     plot_title: str = "Basecalled reads PHRED quality"):
        """
        Plot a distribution of quality scores
        * color
            Color of the area(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x-axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self._require_basecalls()
        fig = self._1D_density_plot(
            field_name="mean_qscore",
            plot_title=plot_title,
            x_lab="Read quality scores",
            color=color,
            x_scale="linear",
            nbins=nbins,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height)
        return fig

    def align_len_1D(
            self,
            color: str = "mediumseagreen",
            nbins: int = 200,
            smooth_sigma: int = 2,
            width: int = None,
            height: int = 500,
            plot_title: str = "Aligned reads length"):
        """
        Plot a distribution of read length(log scale)
        * color
            Color of the area(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x-axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignment information are available
        if not self.has_alignment:
            raise pycoQCError("No Alignment information available")

        fig = self._1D_density_plot(
            field_name="align_len",
            plot_title=plot_title,
            x_lab="Alignment length",
            color=color,
            x_scale="log",
            nbins=nbins,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height)
        return fig

    def identity_freq_1D(
            self,
            color: str = "sandybrown",
            nbins: int = 200,
            smooth_sigma: int = 2,
            width: int = None,
            height: int = 500,
            plot_title: str = "Aligned reads identity"):
        """
        Plot a distribution of alignments identity
        * color
            Color of the area(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x-axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignment information are available
        if not self.has_identity_freq:
            raise pycoQCError("No identity frequency information available")

        fig = self._1D_density_plot(
            field_name="identity_freq",
            plot_title=plot_title,
            x_lab="Identity frequency",
            color=color,
            x_scale="linear",
            nbins=nbins,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height)
        return fig

    def _1D_density_plot(self, field_name, plot_title, x_lab, color, x_scale, nbins, smooth_sigma, width, height):
        """Private function generating density plots for all 1D distribution functions"""
        self.logger.info("\t\tComputing plot")

        # Prepare all data
        prepared = [self._1D_density_data(level, field_name, x_scale, nbins, smooth_sigma)
                    for level in LEVELS]
        _, dd1, ld1 = prepared[0]
        if not len(dd1["x"][0]):
            raise pycoQCError("No {} values to plot".format(field_name))

        # Plot initial data
        common = {
            "mode": "lines+text",
            "hoverinfo": "skip",
            "textposition": 'top center',
            "line":  {'color': 'gray', 'width': 1, 'dash': 'dot'}}
        data = [go.Scatter(x=dd1["x"][0], y=dd1["y"][0], name=dd1["name"][0], fill='tozeroy',
                           fillcolor=color, mode='none', showlegend=True)]
        data += [go.Scatter(x=dd1["x"][i], y=dd1["y"][i], name=dd1["name"][i], text=dd1["text"][i], **common)
                 for i in range(1, len(dd1["x"]))]

        # tweak plot layout
        layout = go.Layout(
            hovermode="closest",
            plot_bgcolor="whitesmoke",
            legend={"x": -0.2, "y": 1, "xanchor": 'left', "yanchor": 'top'},
            updatemenus=self._level_menu(prepared, "update", x=-0.2, y=0,
                                         xanchor='left', yanchor='bottom'),
            width=width,
            height=height,
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"},
            xaxis={"title": x_lab, "type": x_scale, "zeroline": False, "showline": True},
            yaxis={"title": "Read density", "zeroline": False, "showline": True, "fixedrange": True,
                   "range": ld1["yaxis.range"]})

        return go.Figure(data=data, layout=layout)

    def _1D_density_data(self, df_level, field_name, x_scale, nbins, smooth_sigma):
        """Private function preparing data for reads_1D"""

        self.logger.debug("\t\tPreparing data for {} reads and {}".format(df_level, field_name))

        df, _ = self._level_sample(df_level)
        data = df[field_name].to_numpy(dtype=np.float64, na_value=np.nan)
        data = data[np.isfinite(data)]
        names = ["Density", "10%", "25%", "Median", "75%", "90%"]
        label = self._level_label(df_level)

        # An empty read set still sends a complete update. The empty dict it
        # used to send left the previous set's curve on screen, so a run with
        # no failing reads showed its passing reads under "Fail Reads".
        if data.size == 0:
            data_dict = dict(x=[[] for _ in names], y=[[] for _ in names], name=names,
                             text=[""] + [["", ""] for _ in names[1:]])
            return label, data_dict, {"yaxis.range": [0, 1]}

        # Count each categories in log or linear space
        count_y, bins = np.histogram(
            a=data, bins=self._hist_bins(data.min(), data.max(), nbins, x_scale))

        # Label each bin by its centre. The right edge drew every point
        # half a bin too far along: on the log length axis that is a
        # constant +2.6% on a 200 bin plot, and it put the density's
        # right hand tail 26% past the longest read in the run.
        count_x = (bins[:-1] + bins[1:]) / 2

        # Smooth results with a gaussian filter, in floating point: the
        # counts are integers and smoothing in place truncated them, so
        # a small read set lost most of its density and whole bins were
        # flattened to zero.
        count_y = count_y.astype(np.float64)
        if smooth_sigma:
            count_y = gaussian_filter1d(count_y, sigma=smooth_sigma)

        # Percentile markers, each a vertical line labelled at its top
        stat = np.percentile(data, [10, 25, 50, 75, 90])
        y_max = count_y.max()

        data_dict = dict(
            x=[count_x] + [[v, v] for v in stat],
            y=[count_y] + [[0, y_max] for _ in stat],
            name=names,
            text=[""] + [["", "{}<br>{:,.2f}".format(n, v)] for n, v in zip(names[1:], stat)])

        # Make layout dict = Off set for labels on top
        layout_dict = {"yaxis.range": [0, (y_max + y_max / 6) or 1]}
        return label, data_dict, layout_dict

    #~~~~~~~2D DISTRIBUTION METHOD AND HELPER~~~~~~~#
    def read_len_read_qual_2D(self,
                              colorscale=(
                                  [0.0, 'rgba(255,255,255,0)'],
                                  [0.1, 'rgba(255,150,0,0)'],
                                  [0.25, 'rgb(255,100,0)'],
                                  [0.5, 'rgb(200,0,0)'],
                                  [0.75, 'rgb(120,0,0)'],
                                  [1.0, 'rgb(70,0,0)']),
                              x_nbins: int = 200,
                              y_nbins: int = 100,
                              smooth_sigma: int = 2,
                              width: int = None,
                              height: int = 600,
                              plot_title: str = "Basecalled reads length vs reads PHRED quality"):
        """
        Plot a 2D distribution of quality scores vs length of the reads
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/(Not recommanded to change)
        * x_nbins
            Number of bins to divide the read length values in(x-axis)
        * y_nbins
            Number of bins to divide the read quality values in(y-axis)
        * smooth_sigma
            standard deviation for 2D Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self._require_basecalls()
        fig = self._2D_density_plot(
            x_field_name="read_len",
            y_field_name="mean_qscore",
            x_lab="Basecalled length",
            y_lab="PHRED quality scores",
            x_scale="log",
            y_scale="linear",
            x_nbins=x_nbins,
            y_nbins=y_nbins,
            colorscale=colorscale,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height,
            plot_title=plot_title)
        return fig

    def read_len_read_gc_2D(
            self,
            colorscale=(
                    [0.0, 'rgba(255,255,255,0)'],
                    [0.1, 'rgba(204,255,204,0)'],
                    [0.25, 'rgb(102,255,51)'],
                    [0.5, 'rgb(0,153,0)'],
                    [0.75, 'rgb(51,102,0)'],
                    [1.0, 'rgb(0,51,0)']),
            x_nbins: int = 200,
            y_nbins: int = 100,
            smooth_sigma: int = 2,
            width: int = None,
            height: int = 600,
            plot_title: str = "Basecalled reads length vs reads %GC"):
        """
        Plot a 2D distribution of quality scores vs length of the reads
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/(Not recommanded to change)
        * x_nbins
            Number of bins to divide the read length values in(x-axis)
        * y_nbins
            Number of bins to divide the read quality values in(y-axis)
        * smooth_sigma
            standard deviation for 2D Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignment information are available
        if not self.has_gc_percent:
            raise pycoQCError("No %GC information available")

        fig = self._2D_density_plot(
            x_field_name="read_len",
            y_field_name="%GC",
            x_lab="Basecalled length",
            y_lab="%GC",
            x_scale="log",
            y_scale="linear",
            x_nbins=x_nbins,
            y_nbins=y_nbins,
            colorscale=colorscale,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height,
            plot_title=plot_title)

        return fig

    def read_len_align_len_2D(self,
                              colorscale=(
                                      [0.0, 'rgba(255,255,255,0)'],
                                      [0.1, 'rgba(255,150,0,0)'],
                                      [0.25, 'rgb(255,100,0)'],
                                      [0.5, 'rgb(200,0,0)'],
                                      [0.75, 'rgb(120,0,0)'],
                                      [1.0, 'rgb(70,0,0)']),
                              x_nbins: int = 200,
                              y_nbins: int = 100,
                              smooth_sigma: int = 1,
                              width: int = None,
                              height: int = 600,
                              plot_title: str = "Basecalled reads length vs alignments length"):
        """
        Plot a 2D distribution of length of the reads vs length of the alignments
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/(Not recommanded to change)
        * x_nbins
            Number of bins to divide the read length values in(x-axis)
        * y_nbins
            Number of bins to divide the read quality values in(y-axis)
        * smooth_sigma
            standard deviation for 2D Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignment information are available
        if not self.has_alignment:
            raise pycoQCError("No Alignment information available")

        fig = self._2D_density_plot(
            x_field_name="read_len",
            y_field_name="align_len",
            x_lab="Basecalled length",
            y_lab="Alignment length",
            x_scale="log",
            y_scale="log",
            x_nbins=x_nbins,
            y_nbins=y_nbins,
            colorscale=colorscale,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height,
            plot_title=plot_title)
        return fig

    def align_len_identity_freq_2D(
            self,
            colorscale=(
                    [0.0, 'rgba(255,255,255,0)'],
                    [0.1, 'rgba(255,150,0,0)'],
                    [0.25, 'rgb(255,100,0)'],
                    [0.5, 'rgb(200,0,0)'],
                    [0.75, 'rgb(120,0,0)'],
                    [1.0, 'rgb(70,0,0)']),
            x_nbins: int = 200,
            y_nbins: int = 100,
            smooth_sigma: int = 2,
            width: int = None,
            height: int = 600,
            plot_title: str = "Aligned reads length vs alignments identity"):
        """
        Plot a 2D distribution of alignments length vs alignments identity
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/(Not recommanded to change)
        * x_nbins
            Number of bins to divide the read length values in(x-axis)
        * y_nbins
            Number of bins to divide the read quality values in(y-axis)
        * smooth_sigma
            standard deviation for 2D Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignment information are available
        if not self.has_identity_freq:
            raise pycoQCError("No identity frequency information available")

        fig = self._2D_density_plot(
            x_field_name="align_len",
            y_field_name="identity_freq",
            x_lab="Alignment length",
            y_lab="Identity frequency",
            x_scale="log",
            y_scale="linear",
            x_nbins=x_nbins,
            y_nbins=y_nbins,
            colorscale=colorscale,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height,
            plot_title=plot_title)
        return fig

    def read_qual_identity_freq_2D(
            self,
            colorscale=(
                    [0.0, 'rgba(255,255,255,0)'],
                    [0.1, 'rgba(255,150,0,0)'],
                    [0.25, 'rgb(255,100,0)'],
                    [0.5, 'rgb(200,0,0)'],
                    [0.75, 'rgb(120,0,0)'],
                    [1.0, 'rgb(70,0,0)']),
            x_nbins: int = 200,
            y_nbins: int = 100,
            smooth_sigma: int = 1,
            width: int = None,
            height: int = 600,
            plot_title: str = "Reads PHRED quality vs alignments identity"):
        """
        Plot a 2D distribution of read quality vs alignments identity
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/(Not recommanded to change)
        * x_nbins
            Number of bins to divide the read length values in(x-axis)
        * y_nbins
            Number of bins to divide the read quality values in(y-axis)
        * smooth_sigma
            standard deviation for 2D Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignment information are available
        if not self.has_identity_freq:
            raise pycoQCError("No identity frequency information available")

        fig = self._2D_density_plot(
            x_field_name="mean_qscore",
            y_field_name="identity_freq",
            x_lab="PHRED quality",
            y_lab="Identity frequency",
            x_scale="linear",
            y_scale="linear",
            x_nbins=x_nbins,
            y_nbins=y_nbins,
            colorscale=colorscale,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height,
            plot_title=plot_title)
        return fig

    def _2D_density_plot(
            self,
            x_field_name, y_field_name, x_lab, y_lab, x_scale, y_scale, x_nbins, y_nbins,
            colorscale, smooth_sigma, width, height, plot_title):
        """Private function generating density plots for all 2D distribution functions"""
        self.logger.info("\t\tComputing plot")

        # Prepare all data
        prepared = [self._2D_density_data(level, x_field_name, y_field_name, x_nbins, y_nbins,
                                          x_scale, y_scale, smooth_sigma)
                    for level in LEVELS]
        _, dd1 = prepared[0]
        if not len(dd1["x"][0]):
            raise pycoQCError("No reads with both {} and {} to plot".format(
                x_field_name, y_field_name))

        # Plot initial data
        data = [
            go.Contour(x=dd1["x"][0], y=dd1["y"][0], z=dd1["z"][0], contours=dd1["contours"][0],
                       name="Density", hoverinfo="name+x+y", colorscale=colorscale, showlegend=True,
                       connectgaps=True, line={"width": 0}),
            go.Scatter(x=dd1["x"][1], y=dd1["y"][1],
                       mode='markers', name='Median', hoverinfo="name+x+y",
                       marker={"size": 12, "color": 'black', "symbol": "x"})]

        # tweak plot layout
        layout = go.Layout(
            hovermode="closest",
            plot_bgcolor="whitesmoke",
            legend={"x": -0.2, "y": 1, "xanchor": 'left', "yanchor": 'top'},
            updatemenus=self._level_menu(prepared, "restyle", x=-0.2, y=0,
                                         xanchor='left', yanchor='bottom'),
            width=width,
            height=height,
            title={"text":plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"},
            xaxis={"title": x_lab, "showgrid": True, "zeroline": False, "showline": True, "type": x_scale},
            yaxis={"title": y_lab, "showgrid": True, "zeroline": False, "showline": True, "type": y_scale})

        return go.Figure(data=data, layout=layout)

    def _2D_density_data(self, df_level, x_field_name, y_field_name, x_nbins, y_nbins, x_scale, y_scale, smooth_sigma):
        """ Private function preparing data for 2D_density_plot """

        self.logger.debug("\t\tPreparing data for {} reads".format(df_level))

        df, _ = self._level_sample(df_level)
        df = df[[x_field_name, y_field_name]].dropna()
        label = self._level_label(df_level)

        # A complete, empty update for an empty read set, for the same
        # reason as the 1D plots: an empty dict left the last set drawn.
        if df.empty:
            return label, dict(x=[[], []], y=[[], []], z=[[], None], contours=[None, None])

        # Prepare data for x and y
        x_data = df[x_field_name].to_numpy(dtype=np.float64)
        y_data = df[y_field_name].to_numpy(dtype=np.float64)
        x_min, x_med, x_max = np.percentile(x_data, (0, 50, 100))
        y_min, y_med, y_max = np.percentile(y_data, (0, 50, 100))
        x_bins = self._hist_bins(x_min, x_max, x_nbins, x_scale)
        y_bins = self._hist_bins(y_min, y_max, y_nbins, y_scale)

        # Compute 2D histogram
        z, y, x = np.histogram2d(x=y_data, y=x_data, bins=[y_bins, x_bins])
        if smooth_sigma:
            z = gaussian_filter(z, sigma=smooth_sigma)

        # histogram2d returns edges, one more than there are cells, and
        # a contour trace takes the coordinates it is given as cell
        # centres. Handing it the edges therefore drew the whole density
        # half a bin left and down, away from the exact median marker
        # drawn on the same axes.
        x = (x[:-1] + x[1:]) / 2
        y = (y[:-1] + y[1:]) / 2
        # Single precision is ample for a colour and halves the report's
        # share of each density, which it carries four times over.
        z = z.astype(np.float32)
        z_min, z_max = (float(v) for v in np.percentile(z, (0, 100)))
        # A flat density, as from a handful of identical reads, gave a
        # contour step of zero, which plotly cannot draw levels with.
        step = (z_max - z_min) / 15 or 1

        # Extract label and values
        data_dict = dict(
            x=[x, [x_med]], y=[y, [y_med]], z=[z, None],
            contours=[dict(start=z_min, end=z_max, size=step), None])

        return label, data_dict

    #~~~~~~~OUTPUT_OVER_TIME METHODS AND HELPER~~~~~~~#
    def output_over_time(
            self,
            cumulative_color: str = "rgb(204,226,255)",
            interval_color: str = "rgb(102,168,255)",
            time_bins: int = 500,
            width: int = None,
            height: int = 550,
            plot_title: str = "Output over experiment time"):
        """
        Plot a yield over time
        * cumulative_color
            Color of cumulative yield area(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * interval_color
            Color of interval yield line(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * time_bins
            Number of bins to divide the time values in(x-axis)
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self.logger.info("\t\tComputing plot")
        if "start_time" not in self.all_df:
            raise pycoQCError("No start_time information available")

        # Prepare all data. Counting bases needs a sequence, so a run that
        # has not been basecalled gets the read counts on their own.
        count_levels = ["reads", "bases"] if self.has_basecalls else ["reads"]
        prepared = [self._output_over_time_data(df_level=level,
                                                count_level=count,
                                                time_bins=time_bins)
                    for count in count_levels
                    for level in LEVELS]
        _, dd1, ld1 = prepared[0]
        if not len(dd1["x"][0]):
            raise pycoQCError("No reads with a start time to plot")

        # Plot initial data
        common = {
            "mode": "lines+text",
            "hoverinfo": "skip",
            "textposition": 'top center',
            "line":  {'color': 'gray', 'width': 1, 'dash': 'dot'}}
        data = [
            go.Scatter(x=dd1["x"][0], y=dd1["y"][0], name=dd1["name"][0], fill='tozeroy',
                       fillcolor=cumulative_color, mode='none'),
            go.Scatter(x=dd1["x"][1], y=dd1["y"][1], name=dd1["name"][1], mode='lines',
                       line={'color': interval_color, 'width': 2})]
        data += [go.Scatter(x=dd1["x"][i], y=dd1["y"][i], name=dd1["name"][i], text=dd1["text"][i], **common)
                 for i in range(2, len(dd1["x"]))]

        # tweak plot layout
        layout = go.Layout(
            plot_bgcolor="whitesmoke",
            width=width,
            height=height,
            updatemenus=self._level_menu(prepared, "update", x=-0.08, y=0,
                                         xanchor='right', yanchor='bottom'),
            legend={"x": -0.05, "y": 1, "xanchor": 'right', "yanchor": 'top'},
            title={"text": plot_title, "xref": "paper", "x": .5, "xanchor": "center"},
            xaxis={"title": "Experiment time (h)", "zeroline": False, "showline": True},
            yaxis={"title": "Count", "zeroline": False, "showline": True, "fixedrange": True,
                   "range": ld1["yaxis.range"]})

        return go.Figure(data=data, layout=layout)

    def _output_over_time_data(self, df_level, count_level, time_bins=500):
        """Private function preparing data for output_over_time"""
        if "start_time" not in self.all_df:
            raise pycoQCError("No start_time information available")
        self.logger.debug("\t\tPreparing data for {} {}".format(df_level, count_level))

        percents = (50, 75, 90, 99, 100)
        names = ["Cumulative", "Interval"] + ["{}%".format(lab) for lab in percents]
        label = self._level_label(df_level, count_level.capitalize())

        # Reads and bases are counted over the whole run. Scaling a sample
        # by the ratio of read counts is right for reads and wrong for
        # bases, because the sampled reads do not hold that same share of
        # the bases: the lambda run's yield was reported as 14,419,599,424
        # bases against the 14,338,032,127 in the summary table on the same
        # page. Binning the real column is one pass and settles it.
        df, _ = self._level_dfs(df_level)

        # Reads whose start time is missing cannot be placed on the axis;
        # the rest of the run still can. A read set with none left, such as
        # a run with no failing reads, gets a complete empty update.
        df = df.dropna(subset=["start_time"])
        if df.empty:
            data_dict = dict(x=[[] for _ in names], y=[[] for _ in names], name=names,
                             text=["", ""] + [["", ""] for _ in percents])
            return label, data_dict, {"yaxis.range": [0, 1]}

        # Divide the run into time_bins real intervals and label each by
        # its centre, as the channel heatmap does. Digitising against the
        # labels themselves made the first bin hold only the reads at
        # exactly t.min and drew every other bin's count half an interval
        # late.
        t = df["start_time"].to_numpy(dtype=np.float64) / 3600
        x, bin_idx = self._time_bins(t, time_bins)

        # Count reads or bases per categories. Counting reads gives
        # integers, and the smoothing below would truncate them whenever
        # the scaling factor is 1, which is any run small enough not to
        # be sampled, showing gaps in the trace that are not there.
        weights = df["read_len"].to_numpy(dtype=np.float64) if count_level == "bases" else None
        y = np.bincount(bin_idx, weights=weights, minlength=time_bins).astype(np.float64)

        # Transform to cumulative distribution
        y_cum = np.cumsum(y)
        y_cum_max = y_cum[-1]

        # Smooth and rescale interval trace. A read set small enough
        # that smoothing flattens every bin to zero leaves nothing to
        # rescale against, so the trace is left as it is.
        y = gaussian_filter1d(y, sigma=1)
        y_max = y.max() if y.size else 0
        if y_max:
            y = y * y_cum_max / y_max

        # Find percentages of data generated
        lab_text = []
        lab_x = []
        lab_y = []
        for n, lab in enumerate(percents):
            val = y_cum_max * lab / 100
            idx = (np.abs(y_cum-val)).argmin()
            # x comes from a float32 column, and rounding one keeps its
            # full binary expansion, so 16.84 printed as 16.84000015258789.
            lab_text.append(["", '{}%<br>{}h<br>{:,} {}'.format(
                lab, round(float(x[idx]), 2), int(y_cum[idx]), count_level)])
            lab_x.append([x[idx], x[idx]])
            # The last markers sit close together on a long run, so drop
            # alternate labels down the line to keep them legible.
            top = y_cum_max if n % 2 == 0 else y_cum_max * 0.82
            lab_y.append([0, top])

        # make data dict
        data_dict = dict(
            x=[x, x]+lab_x,
            y=[y_cum, y]+lab_y,
            name=names,
            text=["", ""] + lab_text)

        # Make layout dict = offset for labels on top
        layout_dict = {"yaxis.range": [0, (y_cum_max + y_cum_max / 6) or 1]}
        return label, data_dict, layout_dict

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

    #~~~~~~~QUAL_OVER_TIME METHODS AND HELPER~~~~~~~#
    def read_len_over_time(
            self,
            median_color: str = "rgb(102,168,255)",
            quartile_color: str = "rgb(153,197,255)",
            extreme_color: str = "rgba(153,197,255,0.5)",
            smooth_sigma: int = 1,
            time_bins: int = 500,
            width: int = None,
            height: int = 550,
            plot_title: str = "Read length over experiment time"):
        """
        Plot a read length over time
        * median_color
            Color of median line color(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-col
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in(x-axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self._require_basecalls()
        fig = self._over_time_plot(
            field_name="read_len",
            plot_title=plot_title,
            y_lab="Read length",
            y_scale="log",
            median_color=median_color,
            quartile_color=quartile_color,
            extreme_color=extreme_color,
            smooth_sigma=smooth_sigma,
            time_bins=time_bins,
            width=width,
            height=height)
        return fig

    def read_qual_over_time(
            self,
            median_color: str = "rgb(250,128,114)",
            quartile_color: str = "rgb(250,170,160)",
            extreme_color: str = "rgba(250,170,160,0.5)",
            smooth_sigma: int = 1,
            time_bins: int = 500,
            width: int = None,
            height: int = 500,
            plot_title: str = "Read quality over experiment time"):
        """
        Plot a mean quality over time
        * median_color
            Color of median line color(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-col
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in(x-axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self._require_basecalls()

        fig = self._over_time_plot(
            field_name="mean_qscore",
            plot_title=plot_title,
            y_lab="Mean read PHRED quality",
            y_scale="linear",
            median_color=median_color,
            quartile_color=quartile_color,
            extreme_color=extreme_color,
            smooth_sigma=smooth_sigma,
            time_bins=time_bins,
            width=width,
            height=height)
        return fig

    def read_gc_over_time(
            self,
            median_color: str = "rgb(51,153,0)",
            quartile_color: str = "rgb(85,255,0)",
            extreme_color: str = "rgba(153,255,102,0.5)",
            smooth_sigma: int = 1,
            time_bins: int = 500,
            width: int = None,
            height: int = 500,
            plot_title: str = "%GC over experiment time"):
        """
        Plot a mean quality over time
        * median_color
            Color of median line color(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-col
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in(x-axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Check data contains %GC information
        if not self.has_gc_percent:
            raise pycoQCError("No %GC information available")

        fig = self._over_time_plot(
            field_name="%GC",
            plot_title=plot_title,
            y_lab="Mean read %GC",
            y_scale="linear",
            median_color=median_color,
            quartile_color=quartile_color,
            extreme_color=extreme_color,
            smooth_sigma=smooth_sigma,
            time_bins=time_bins,
            width=width,
            height=height)

        return fig

    def align_len_over_time(
            self,
            median_color: str = "rgb(102,168,255)",
            quartile_color: str = "rgb(153,197,255)",
            extreme_color: str = "rgba(153,197,255,0.5)",
            smooth_sigma: int = 1,
            time_bins:int = 500,
            width: int = None,
            height: int = 500,
            plot_title: str = "Aligned reads length over experiment time"):
        """
        Plot a aligned reads length over time
        * median_color
            Color of median line color(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-col
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in(x-axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignment information are available
        if not self.has_alignment:
            raise pycoQCError("No Alignment information available")

        fig = self._over_time_plot(
            field_name="align_len",
            plot_title=plot_title,
            y_lab="Aligned reads length",
            y_scale="log",
            median_color=median_color,
            quartile_color=quartile_color,
            extreme_color=extreme_color,
            smooth_sigma=smooth_sigma,
            time_bins=time_bins,
            width=width,
            height=height)
        return fig

    def identity_freq_over_time(
            self,
            median_color: str = "rgb(250,128,114)",
            quartile_color: str = "rgb(250,170,160)",
            extreme_color: str = "rgba(250,170,160,0.5)",
            smooth_sigma: int = 1,
            time_bins: int = 500,
            width: int = None,
            height: int = 500,
            plot_title: str = "Aligned reads identity over experiment time"):
        """
        Plot the alignment identity scores over time
        * median_color
            Color of median line color(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-col
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in(x-axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """

        # Verify that alignment information are available
        if not self.has_identity_freq:
            raise pycoQCError("No identity frequency information available")

        fig = self._over_time_plot(
            field_name="identity_freq",
            plot_title=plot_title,
            y_lab="Identity frequency",
            y_scale="linear",
            median_color=median_color,
            quartile_color=quartile_color,
            extreme_color=extreme_color,
            smooth_sigma=smooth_sigma,
            time_bins=time_bins,
            width=width,
            height=height)
        return fig

    def _over_time_plot(self, field_name, plot_title, y_lab, y_scale, median_color, quartile_color, extreme_color,
                        smooth_sigma, time_bins, width, height):
        """Private function generating density plots for all over_time functions"""
        self.logger.info("\t\tComputing plot")
        if "start_time" not in self.all_df:
            raise pycoQCError("No start_time information available")

        prepared = [self._over_time_data(df_level=level, field_name=field_name,
                                         smooth_sigma=smooth_sigma, time_bins=time_bins)
                    for level in LEVELS]
        _, dd1 = prepared[0]
        if not len(dd1["x"][0]):
            raise pycoQCError("No reads with both a start time and {} to plot".format(field_name))

        # Plot initial data
        common = {
            "mode": "lines",
            "connectgaps": True}
        data = [
            go.Scatter(x=dd1["x"][0], y=dd1["y"][0], name=dd1["name"][0], line={"color": extreme_color},
                       legendgroup="Extreme", **common),
            go.Scatter(x=dd1["x"][1], y=dd1["y"][1], name=dd1["name"][1], fill="tonexty",
                       line={"color": extreme_color}, legendgroup="Extreme", **common),
            go.Scatter(x=dd1["x"][2], y=dd1["y"][2], name=dd1["name"][2], line={"color": quartile_color},
                       legendgroup="Quartiles", **common),
            go.Scatter(x=dd1["x"][3], y=dd1["y"][3], name=dd1["name"][3], fill="tonexty",
                       line={"color": quartile_color}, legendgroup="Quartiles", **common),
            go.Scatter(x=dd1["x"][4], y=dd1["y"][4], name=dd1["name"][4], line={"color": median_color}, **common)]

        # tweak plot layout
        layout = go.Layout(
            plot_bgcolor="whitesmoke",
            width=width,
            height=height,
            updatemenus=self._level_menu(prepared, "restyle", x=-0.07, y=0,
                                         xanchor='right', yanchor='bottom'),
            legend={"x": -0.07, "y": 1, "xanchor": 'right', "yanchor": 'top'},
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"},
            yaxis={"title": y_lab, "zeroline": False, "type": y_scale, "showline": True,
                   "rangemode": 'nonnegative', "fixedrange": True},
            xaxis={"title": "Experiment time (h)", "zeroline": False, "showline": True, "rangemode": 'nonnegative'})

        return go.Figure(data=data, layout=layout)

    def _over_time_data(self, df_level, field_name="read_len", smooth_sigma=1.5, time_bins=500):
        """Private function preparing data for qual_over_time"""
        if "start_time" not in self.all_df:
            raise pycoQCError("No start_time information available")
        self.logger.debug("\t\tPreparing data for {} reads and {}".format(df_level, field_name))

        val_name = ["Min", "Max", "25%", "75%", "Median"]
        quantiles = [0.0, 1.0, 0.25, 0.75, 0.5]
        label = self._level_label(df_level)

        # Drop rows missing either value together, so that each value stays
        # paired with its own start time. Dropping only the value column left
        # the two arrays different lengths, and zip() then silently paired
        # every value with the wrong time.
        df, _ = self._level_sample(df_level)
        subset = df[["start_time", field_name]].dropna()
        if subset.empty:
            return label, dict(x=[[] for _ in val_name], y=[[] for _ in val_name], name=val_name)

        # Same binning as the other over time plots: real intervals,
        # each labelled by its centre.
        t = subset["start_time"].to_numpy(dtype=np.float64) / 3600
        x, bin_idx = self._time_bins(t, time_bins)

        # Aggregate values per category. Grouping in pandas avoids building
        # one Python list per bin and keeps the percentiles vectorised.
        grouped = (pd.Series(subset[field_name].to_numpy(dtype=np.float64))
                   .groupby(bin_idx)
                   .quantile(quantiles)
                   .unstack()
                   .reindex(np.arange(x.size)))
        stat_dict = {name: grouped[q].to_numpy()
                     for name, q in zip(val_name, quantiles)}

        # Values smoothing
        if smooth_sigma:
            for val in val_name:
                stat_dict[val] = self._smooth_with_gaps(stat_dict[val], smooth_sigma)

        # make data dict
        data_dict = dict(
            x=[x for _ in val_name],
            y=[stat_dict[val] for val in val_name],
            name=val_name)
        return label, data_dict

    #~~~~~~~BARCODE_COUNT METHODS AND HELPER~~~~~~~#
    def barcode_counts(
            self,
            colors: list = ("#f8bc9c", "#f6e9a1", "#f5f8f2", "#92d9f5", "#4f97ba"),
            width: int = None,
            height: int = None,
            plot_title: str = "Percentage of reads per barcode"):
        """
        Plot a mean quality over time
        * colors
            List of colors(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel. Grows with the number of
            barcodes unless set.
        * plot_title
            Title to display on top of the plot
        """
        # Verify that barcode information are available
        if not self.has_barcodes:
            raise pycoQCError("No barcode information available")
        self.logger.info("\t\tComputing plot")

        # Prepare all data
        levels = [self._barcode_counts_data(df_level=level) for level in LEVELS]
        _, dd1 = levels[0]

        n_barcodes = len(dd1["labels"][0])

        # A 96 barcode kit needs more room than a 12 barcode one: at a fixed
        # height the legend and the table were simply cut off partway down.
        if height is None:
            height = max(500, 290 + 30 * n_barcodes)

        # A pie divided 96 ways is a ring of slivers under a thicket of
        # leader lines, and its legend is longer than the plot is tall. Past
        # a couple of dozen barcodes a bar for each one is readable, and
        # stays readable for the whole 96 of the largest kit.
        as_pie = n_barcodes <= PIE_MAX_BARCODES

        # One table and one pie per read set, with only the chosen pair
        # visible. The buttons used to restyle a single pair in place, and
        # because each carried both args and args2 plotly treated it as a
        # toggle: the first click undid the previous button rather than
        # applying the one just pressed, so every switch took two clicks.
        # They also only ever restyled the table, leaving the pie showing
        # the read set that had been selected before.

        data = []
        for i, (_, dd) in enumerate(levels):
            data.append(go.Table(
                header={
                    "values": ["Barcode", "Count"],
                    "align": "center", "fill": {"color": "grey"},
                    "font": {"size": 14, "color": "white"},
                    "height": 40},
                cells={
                    "values": [dd['labels'][0], dd['values'][0]],
                    "format": ["", ","],
                    "align": "center",
                    "fill": {"color": "whitesmoke"},
                    "font": {"size": 12}, "height": 30},
                columnwidth=[3, 2],
                domain={"x": [0.5, 1.0], "y": [0, 1.0]},
                visible=(i == 0)))
            if as_pie:
                data.append(go.Pie(
                    labels=dd["labels"][0],
                    values=dd["values"][0],
                    sort=False,
                    domain={"x": [0, 0.4], "y": [0, 1.0]},
                    visible=(i == 0)))
            else:
                data.append(go.Bar(
                    y=list(dd["labels"][0]),
                    x=list(dd["values"][0]),
                    orientation="h",
                    marker={"color": colors[0]},
                    hovertemplate="%{y}: %{x:,} reads<extra></extra>",
                    visible=(i == 0)))

        buttons = []
        for i, (label, _) in enumerate(levels):
            visible = [False] * len(data)
            visible[2 * i] = visible[2 * i + 1] = True
            buttons.append(dict(label=label, method="restyle",
                                args=[{"visible": visible}]))

        # Laid out in a row under the pie. Down the left hand side they sat
        # on top of the legend, which grows with the number of barcodes.
        updatemenus = [
            dict(type="buttons", active=0, direction="right",
                 x=0.2, y=-0.02, xanchor="center", yanchor="top",
                 pad={"t": 4}, buttons=buttons)]

        # tweak plot layout
        layout = go.Layout(
            plot_bgcolor="whitesmoke",
            legend={"x": -0.2, "y": 1, "xanchor": 'left', "yanchor": 'top'},
            updatemenus=updatemenus,
            width=width,
            height=height,
            # Room under the plot for the row of buttons.
            margin={"b": 70},
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"})

        if not as_pie:
            # Bars carry their own labels, so the legend would only repeat
            # them down the side.
            layout.update(
                showlegend=False,
                bargap=0.2,
                xaxis={"domain": [0.06, 0.42], "title": "Reads",
                       "zeroline": False, "showline": True},
                yaxis={"domain": [0, 1], "autorange": "reversed",
                       "tickfont": {"size": 10}})

        return go.Figure(data=data, layout=layout)

    def _barcode_counts_data(self, df_level):
        """Private function preparing data for barcode_counts"""
        self.logger.debug("\t\tPreparing data for {} reads".format(df_level))

        # Counting is cheap, so these are the run's own counts rather than
        # a sample scaled back up, which was unbiased but noisy and did not
        # match the barcode table on the same page.
        df, _ = self._level_dfs(df_level)
        counts = df["barcode"].value_counts()
        counts = counts.reindex(self._barcode_order(counts.index))

        # Extract label and values
        data_dict = dict(
            labels=[counts.index],
            values=[counts.values])

        return self._level_label(df_level), data_dict

    #~~~~~~~RUN HEALTH METHODS AND HELPER~~~~~~~#
    def read_speed_1D(self,
                      color: str = "plum",
                      nbins: int = 200,
                      smooth_sigma: int = 2,
                      width: int = None,
                      height: int = 500,
                      plot_title: str = "Translocation speed"):
        """
        Plot a distribution of how fast each strand passed through its pore,
        in bases per second. A kit is designed for one speed, 400 bases per
        second for most current DNA kits, and a run well below it points at
        the motor protein, the temperature or the flow cell.
        * color
            Color of the area(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * nbins
            Number of bins to devide the x-axis in
        * smooth_sigma
            standard deviation for Gaussian kernel
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self._require_speed()
        return self._1D_density_plot(
            field_name="speed",
            plot_title=plot_title,
            x_lab=self._speed_label,
            color=color,
            x_scale="linear",
            nbins=nbins,
            smooth_sigma=smooth_sigma,
            width=width,
            height=height)

    def read_speed_over_time(
            self,
            median_color: str = "rgb(153,102,204)",
            quartile_color: str = "rgb(204,170,230)",
            extreme_color: str = "rgba(204,170,230,0.5)",
            smooth_sigma: int = 1,
            time_bins: int = 500,
            width: int = None,
            height: int = 500,
            plot_title: str = "Translocation speed over experiment time"):
        """
        Plot translocation speed over experiment time. A slow decline is
        normal as the run ages; a sharp drop is worth investigating.
        * median_color
            Color of median line color(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * quartile_color
            Color of inter quartile area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * extreme_color
            Color of inter extreme area and lines(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in(x-axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        self._require_speed()
        return self._over_time_plot(
            field_name="speed",
            plot_title=plot_title,
            y_lab=self._speed_label,
            y_scale="linear",
            median_color=median_color,
            quartile_color=quartile_color,
            extreme_color=extreme_color,
            smooth_sigma=smooth_sigma,
            time_bins=time_bins,
            width=width,
            height=height)

    def _require_speed(self):
        if not self.has_speed:
            raise pycoQCError(
                "No translocation speed available: it needs a read length "
                "and a duration for each read")

    # Colours for the end reasons MinKNOW records, so that one reason looks
    # the same in every report. Anything else takes the next colour along.
    END_REASON_COLORS = {
        "signal_positive": "#4f97ba",
        "signal_negative": "#92d9f5",
        "unblock_mux_change": "#f44f39",
        "data_service_unblock_mux_change": "#fc8161",
        "mux_change": "#f6c344",
        "api_request": "#9b6fc3",
        "device_data_error": "#5a5a5a",
        "analysis_config_change": "#b5b5b5",
        "paused": "#d9d9d9",
        "partial": "#c7a27c",
        "unknown": "#828282"}
    OTHER_COLORS = ("#7fbf7b", "#e7a1c4", "#a6761d", "#1b9e77", "#666666")

    def end_reason_over_time(
            self,
            time_bins: int = 100,
            width: int = None,
            height: int = 550,
            plot_title: str = "Read end reason over experiment time"):
        """
        Plot how reads ended, as the share of reads per time interval.
        signal_positive is a strand that finished normally. unblock_mux_change
        is MinKNOW ejecting a blocked strand, and a share that climbs through
        the run means pores are clogging. data_service_unblock_mux_change is
        adaptive sampling rejecting a read, and mux_change a read cut short by
        a mux scan.
        * time_bins
            Number of bins to divide the time values in(x-axis)
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        if not self.has_end_reason:
            raise pycoQCError("No end_reason information available")
        if "start_time" not in self.all_df:
            raise pycoQCError("No start_time information available")
        self.logger.info("\t\tComputing plot")

        # The same reasons, in the same order and colours, for every read
        # set, most common first.
        reasons = list(self.all_df["end_reason"].fillna("unknown").value_counts().index)
        others = iter(self.OTHER_COLORS * len(reasons))
        colors = [self.END_REASON_COLORS.get(r) or next(others) for r in reasons]

        prepared = [self._end_reason_data(level, reasons, time_bins) for level in LEVELS]
        _, dd1 = prepared[0]
        if not len(dd1["x"][0]):
            raise pycoQCError("No reads with a start time to plot")

        data = [go.Scatter(x=dd1["x"][i], y=dd1["y"][i], customdata=dd1["customdata"][i],
                           name=dd1["name"][i], mode="lines", stackgroup="reasons",
                           line={"width": 0.5, "color": color}, fillcolor=color,
                           hovertemplate="%{y:.1f}% (%{customdata:,} reads)")
                for i, color in enumerate(colors)]

        layout = go.Layout(
            plot_bgcolor="whitesmoke",
            width=width,
            height=height,
            hovermode="x unified",
            updatemenus=self._level_menu(prepared, "restyle", x=-0.07, y=0,
                                         xanchor='right', yanchor='bottom'),
            # Under the plot, in a row: the stack fills the whole area, and
            # the reasons' names are long enough to run into the axis title
            # down the side.
            legend={"orientation": "h", "x": 0.5, "y": -0.18, "xanchor": "center",
                    "yanchor": "top", "traceorder": "normal"},
            margin={"b": 110},
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"},
            xaxis={"title": "Experiment time (h)", "zeroline": False, "showline": True},
            yaxis={"title": "Share of reads (%)", "range": [0, 100], "zeroline": False,
                   "showline": True, "fixedrange": True})

        return go.Figure(data=data, layout=layout)

    def _end_reason_data(self, df_level, reasons, time_bins=100):
        """Share and count of reads per end reason and time interval.

        Counted over the whole read set: a rare reason can be absent from a
        sample, and counting is cheap.
        """
        self.logger.debug("\t\tPreparing data for {} reads".format(df_level))
        label = self._level_label(df_level)

        df, _ = self._level_dfs(df_level)
        df = df.dropna(subset=["start_time"])
        if df.empty:
            empty = [[] for _ in reasons]
            return label, dict(x=empty, y=empty, customdata=empty, name=list(reasons))

        x, idx = self._time_bins(df["start_time"].to_numpy(dtype=np.float64) / 3600, time_bins)
        codes = pd.Categorical(df["end_reason"].fillna("unknown"), categories=reasons).codes
        counts = np.zeros((len(reasons), time_bins), dtype=np.int64)
        np.add.at(counts, (codes, idx), 1)

        # An interval without reads has no shares, rather than shares of 0
        # that would draw the stack down to nothing.
        total = counts.sum(axis=0)
        with np.errstate(invalid="ignore", divide="ignore"):
            share = np.where(total > 0, 100 * counts / total, np.nan)
        overall = 100 * counts.sum(axis=1) / counts.sum()

        return label, dict(
            x=[x for _ in reasons],
            y=[row for row in share],
            customdata=[row for row in counts],
            name=["{} ({:.1f}%)".format(r, o) for r, o in zip(reasons, overall)])

    def active_channels_over_time(
            self,
            total_color: str = "rgb(102,168,255)",
            interval_hours: float = 1.5,
            width: int = None,
            height: int = 500,
            plot_title: str = "Active channels over experiment time"):
        """
        Plot how many channels produced at least one read in each time
        interval, against the number the flow cell has. The decline over a
        run is how fast its pores are being used up. Where the mux is
        recorded, each channel is also counted under the mux its read came
        from: MinKNOW moves channels between their four pores, and one group
        dying early shows as one line falling away from the others.
        * total_color
            Color of the total line and area(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * interval_hours
            Length of the intervals channels are counted in. A channel
            counts once per interval it produced a read in, so the count
            depends on it. MinKNOW stops every channel for a mux scan every
            hour and a half by default, and intervals that do not match that
            cycle catch one scan in some and two in others: on a real run,
            0.72 hour intervals swung by half the count from one to the next,
            and 1.5 hour ones by a tenth. A run shorter than ten intervals is
            divided into ten.
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        if "channel" not in self.all_df or "start_time" not in self.all_df:
            raise pycoQCError("No channel or start_time information available")
        self.logger.info("\t\tComputing plot")

        muxes = []
        if self.has_mux:
            muxes = sorted(int(m) for m in self.all_df["mux"].dropna().unique())
        n_channels = self._flow_cell_channels

        # One set of intervals for every read set, from the whole run, so
        # switching between them moves the lines and not the axis.
        t = self.all_df["start_time"].dropna().to_numpy(dtype=np.float64) / 3600
        if t.size == 0:
            raise pycoQCError("No reads with both a channel and a start time to plot")
        span = float(t.max() - t.min())
        step = min(interval_hours, span / 10) if span > 0 else interval_hours
        intervals = (float(t.min()), step, max(1, int(round(span / step))))

        prepared = [self._active_channels_data(level, muxes, n_channels, intervals)
                    for level in LEVELS]
        _, dd1 = prepared[0]
        if not len(dd1["x"][0]):
            raise pycoQCError("No reads with both a channel and a start time to plot")

        data = [go.Scatter(x=dd1["x"][0], y=dd1["y"][0], name=dd1["name"][0], mode="lines",
                           fill="tozeroy", line={"color": total_color, "width": 2})]
        palette = ("#f44f39", "#f6c344", "#1b9e77", "#9b6fc3", "#828282")
        data += [go.Scatter(x=dd1["x"][i], y=dd1["y"][i], name=dd1["name"][i], mode="lines",
                            line={"color": palette[(i - 1) % len(palette)], "width": 1})
                 for i in range(1, len(muxes) + 1)]
        data.append(go.Scatter(x=dd1["x"][-1], y=dd1["y"][-1], name=dd1["name"][-1],
                               mode="lines", hoverinfo="skip",
                               line={"color": "gray", "width": 1, "dash": "dot"}))

        layout = go.Layout(
            plot_bgcolor="whitesmoke",
            width=width,
            height=height,
            hovermode="x unified",
            updatemenus=self._level_menu(prepared, "restyle", x=-0.07, y=0,
                                         xanchor='right', yanchor='bottom'),
            # Inside the plot, top right, where the count has fallen away
            # from the flow cell line by then. To the left of the axis, the
            # longer names ran into its title.
            legend={"x": 0.99, "y": 0.97, "xanchor": 'right', "yanchor": 'top',
                    "bgcolor": "rgba(255,255,255,0.7)"},
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"},
            xaxis={"title": "Experiment time (h)", "zeroline": False, "showline": True},
            yaxis={"title": "Channels producing reads per {:g} h".format(intervals[1])
                   if intervals[1] >= 0.1 else "Channels producing reads",
                   "rangemode": "tozero",
                   "zeroline": False, "showline": True, "fixedrange": True})

        return go.Figure(data=data, layout=layout)

    def _active_channels_data(self, df_level, muxes, n_channels, intervals):
        """Distinct channels with a read in each interval, in total and per
        mux, and the flow cell's channel count to compare them with.

        Counted over the whole read set: a sample would miss the channels
        that produced only a read or two, which is exactly the tail that
        matters late in a run.
        """
        self.logger.debug("\t\tPreparing data for {} reads".format(df_level))
        label = self._level_label(df_level)
        names = (["All channels"] + ["Mux {}".format(m) for m in muxes]
                 + ["Flow cell ({:,} channels)".format(n_channels)])

        df, _ = self._level_dfs(df_level)
        df = df.dropna(subset=["start_time", "channel"])
        if df.empty:
            empty = [[] for _ in names]
            return label, dict(x=empty, y=empty, name=names)

        # Each interval labelled by its centre; the end of the run that is
        # left over joins the last one rather than being drawn as a short
        # interval that looks like a collapse.
        start, step, time_bins = intervals
        x = start + step * (np.arange(time_bins) + 0.5)
        t = df["start_time"].to_numpy(dtype=np.float64) / 3600
        idx = np.clip(np.floor((t - start) / step).astype(np.int64), 0, time_bins - 1)
        channel = df["channel"].to_numpy(dtype=np.int64)

        def distinct(mask):
            # One entry per (interval, channel) pair, then counted per interval.
            pairs = np.unique(idx[mask] * (n_channels + 1) + channel[mask])
            return np.bincount(pairs // (n_channels + 1), minlength=time_bins)

        ys = [distinct(np.ones(len(df), dtype=bool))]
        if muxes:
            mux = df["mux"].to_numpy(dtype=np.float64, na_value=np.nan)
            ys += [distinct(mux == m) for m in muxes]
        ys.append(np.full(2, n_channels))
        xs = [x for _ in ys[:-1]] + [[x[0], x[-1]]]

        return label, dict(x=xs, y=ys, name=names)

    #~~~~~~~CHANNEL ACTIVITY METHODS AND HELPER~~~~~~~#
    def channels_activity(
            self,
            colorscale: list = (
                    [0.0, 'rgba(255,255,255,0)'],
                    [0.01, 'rgb(255,255,200)'],
                    [0.25, 'rgb(255,200,0)'],
                    [0.5, 'rgb(200,0,0)'],
                    [0.75, 'rgb(120,0,0)'],
                    [1.0, 'rgb(0,0,0)']),
            smooth_sigma: int = 1,
            time_bins: int = 100,
            width: int = None,
            height: int = 600,
            plot_title: str = "Output per channel over experiment time"):
        """
        Plot a yield over time
        * colorscale
            a valid plotly color scale https://plot.ly/python/colorscales/(Not recommended to change)
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * time_bins
            Number of bins to divide the time values in(y-axis)
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        if "channel" not in self.all_df or "start_time" not in self.all_df:
            raise pycoQCError("No channel or start_time information available")

        self.logger.info("\t\tComputing plot")

        # A channel numbered past the platform's count raised IndexError
        # and took the report down; the map is widened to fit instead.
        n_channels = self._flow_cell_channels

        # Prepare all data
        count_levels = ["reads", "bases"] if self.has_basecalls else ["reads"]
        prepared = [self._channels_activity_data(
                        df_level=level, count_level=count,
                        n_channels=n_channels, smooth_sigma=smooth_sigma,
                        time_bins=time_bins)
                    for count in count_levels
                    for level in LEVELS]
        _, dd1 = prepared[0]
        if not len(dd1["y"][0]):
            raise pycoQCError("No reads with both a channel and a start time to plot")

        # Plot initial data
        data = [go.Heatmap(x=dd1["x"][0], y=dd1["y"][0], z=dd1["z"][0], xgap=0.5, colorscale=colorscale,
                           hoverinfo="x+y+z")]

        # tweak plot layout
        layout = go.Layout(
            plot_bgcolor="whitesmoke",
            width=width,
            height=height,
            updatemenus=self._level_menu(prepared, "restyle", x=-0.06, y=0,
                                         xanchor='right', yanchor='bottom'),
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"},
            xaxis={"title": "Channel id", "zeroline": False, "showline": False, "nticks": 20, "showgrid": False},
            yaxis={"title": "Experiment time (h)", "zeroline": False, "showline": False,
                   "hoverformat": ".2f", "fixedrange": True})

        return go.Figure(data=data, layout=layout)

    def _channels_activity_data(self, df_level, count_level="bases", n_channels=512, smooth_sigma=2, time_bins=150):
        """Private function preparing data for channels_activity"""
        if "channel" not in self.all_df or "start_time" not in self.all_df:
            raise pycoQCError("No channel or start_time information available")
        self.logger.debug("\t\tPreparing data for {} {}".format(df_level, count_level))

        df, sf = self._level_sample(df_level)
        df = df.dropna(subset=["start_time", "channel"])
        label = self._level_label(df_level, count_level.capitalize())
        if df.empty:
            return label, dict(x=[[]], y=[[]], z=[[]])

        # Bin data in categories. Binning against the bin values left a
        # first row reachable only by reads at exactly t.min, and
        # labelling each row with the end of its interval shifted the
        # whole picture half a bin, since plotly draws bricks centred on
        # the y values. Dividing the span into time_bins real intervals
        # and labelling each by its centre avoids both.
        y, t_idx = self._time_bins(df["start_time"].to_numpy(dtype=np.float64) / 3600, time_bins)

        # Channels are numbered from 1. A channel 0 indexed column -1, so
        # its reads were silently drawn on the last channel of the flow cell.
        channel = df["channel"].to_numpy(dtype=np.int64)
        valid = (channel >= 1) & (channel <= n_channels)
        if count_level == "bases":
            weights = df["read_len"].to_numpy(dtype=np.float64)
        else:
            weights = np.ones(len(df), dtype=np.float64)

        # Count values per categories. This started from ones, so after
        # scaling every empty cell reported sf reads or bases on hover
        # rather than none. The colour scale already makes zero
        # transparent, so nothing is gained by the offset.
        z = np.zeros((time_bins, n_channels), dtype=np.float64)
        np.add.at(z, (t_idx[valid], channel[valid] - 1), weights[valid])

        # Scale counts in case of downsampling
        z *= sf

        # Time series smoothing
        if smooth_sigma:
            z = gaussian_filter1d(z, sigma=smooth_sigma, axis=0)

        # Define x and y-axis. digitize(right=True) returns 0 for reads
        # at t.min, so row i holds the reads at bins[i] and y needs every
        # bin. Dropping the first left y one shorter than z, and plotly
        # then treats the values as brick centres and extrapolates,
        # drawing every row one bin later than its data.
        x = ["c {}".format(i) for i in range(1, n_channels+1)]

        # Single precision is ample for a colour, and halves what the
        # report has to carry for a 3000 channel flow cell.
        return label, dict(x=[x], y=[y], z=[z.astype(np.float32)])

    #~~~~~~~ALIGNMENT_SUMMARY METHOD~~~~~~~#
    def alignment_reads_status(
            self,
            colors: list = ("#f44f39", "#fc8161", "#fcaf94", "#fdd5c4", "#828282"),
            width: int = None,
            height: int = 500,
            plot_title: str = "Summary of reads alignment status"):
        """
        Plot a basic alignment summary
        * colors
            List of colors(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignment information are available
        if not self.has_alignment:
            raise pycoQCError("No Alignment information available")
        self.logger.info("\t\tComputing plot")

        df = self.alignments_df
        # Create empty multiplot figure
        fig = make_subplots(rows=1, cols=2, column_widths=[0.4, 0.6], specs=[[{"type": "table"}, {"type": "pie"}]])

        # plot Table
        data = go.Table(
            columnwidth=[3, 2, 2],
            header={"values": list(df.columns), "align": "center", "fill_color": "grey", "font_size": 14,
                    "font_color": "white", "height": 40},
            cells={"values": df.values.T, "align": "center", "fill_color": "whitesmoke", "font_size": 12, "height": 30})
        fig.add_trace(data, row=1, col=1)

        # plot Pie plot. Each colour belongs to a category, in the parser's
        # fixed order, so a category keeps its colour when another one is
        # absent from the run rather than taking the next one along.
        order = ["Primary", "Secondary", "Supplementary", "Duplicated", "Unmapped"]
        slice_colors = [colors[order.index(c) % len(colors)] if c in order else None
                        for c in df["Alignments"]]
        data = go.Pie(
            labels=df["Alignments"],
            values=df["Counts"],
            sort=False,
            marker={"colors": slice_colors},
            name="Pie plot",
            textinfo='label+percent')
        fig.add_trace(data, row=1, col=2)

        # Change the layout
        fig.update_layout(
            width=width,
            height=height,
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"})

        return fig

    #~~~~~~~ALIGNMENT RATE METHOD AND HELPER~~~~~~~#
    def alignment_rate(
            self,
            colors: list = ("#fcaf94", "#828282", "#fc8161", "#828282", "#f44f39", "#d52221", "#828282",
                            "#828282", "#828282", "#828282"),
            width: int = None,
            height: int = 600,
            plot_title: str = "Bases alignment rate"):
        """
        Plot a basic alignment summary
        * colors
            List of colors(hex, rgb, rgba, hsl, hsv or any CSS named colors https://www.w3.org/TR/css-color-3/#svg-color
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """

        # Verify that alignment information are available
        if not self.has_identity_freq:
            raise pycoQCError("No identity frequency information available")
        self.logger.info("\t\tComputing plot")

        # Extract Data
        bc_bases = self.all_df["read_len"].sum()
        s = self.all_df[["read_len", "align_len", "insertion", "deletion", "soft_clip", "mismatch"]].dropna().sum()
        total_error = s["insertion"] + s["deletion"] + s["mismatch"]
        matching = s["align_len"] - total_error
        unmapped = bc_bases - s["read_len"]
        ct = namedtuple("ct", ["Bases", "Counts", "Total_freq", "Aligned_freq"])
        l = [
            ct("Basecalled", bc_bases, 1, 1),
            ct("Unmapped reads", unmapped, unmapped / bc_bases, 1),
            ct("Mapped reads", s["read_len"], s["read_len"] / bc_bases, 1),
            ct("Softclip", s["soft_clip"], s["soft_clip"] / bc_bases, 1),
            ct("Aligned", s["align_len"], s["align_len"] / bc_bases, 1),
            ct("Matching", matching, matching / bc_bases, matching / s["align_len"]),
            ct("Non-matching", total_error, total_error / bc_bases, total_error / s["align_len"]),
            ct("Insertions", s["insertion"], s["insertion"]/bc_bases, s["insertion"] / s["align_len"]),
            ct("Deletions", s["deletion"], s["deletion"] / bc_bases, s["deletion"] / s["align_len"]),
            ct("Mismatches", s["mismatch"], s["mismatch"] / bc_bases, s["mismatch"] / s["align_len"])]

        # Cast to df and compute percentage
        df = pd.DataFrame(l)

        # plot Table
        data1 = go.Table(
            columnwidth=[3, 2, 2, 2],
            header={
                "values": ["Bases", "Bases Count", "% Total", "% Aligned"],
                "align": "center", "fill_color": ["grey"], "font_size": 14, "font_color": "white", "height": 40},
            cells={
                "values": df.values.T, "format": ["", ".3e", ".3p", ".3p"],
                "align": "center", "fill_color": "whitesmoke", "font_size": 12, "height": 30})

        data2 = go.Sankey(
            arrangement="freeform",
            node=go.sankey.Node(
                pad=20,
                thickness=30,
                line={"width": 0},
                label=df["Bases"].values,
                x=[0, 0.25, 0.25, 0.5, 0.5, 1, 0.75],
                y=[0, 1, 0, 0.9, 0, 0, 0.8],
                color=colors),
            link=go.sankey.Link(
                source=[0, 0, 2, 2, 4, 4, 6, 6, 6],
                target=[1, 2, 3, 4, 5, 6, 7, 8, 9],
                value=df["Counts"].values[1:],
                hoverinfo="none"))

        # Create multipanel figure
        fig = make_subplots(rows=1, cols=2, column_widths=[0.4, 0.6], specs=[[{"type": "table"}, {"type": "sankey"}]])
        fig.add_trace(data1, row=1, col=1)
        fig.add_trace(data2, row=1, col=2)
        fig.update_layout(
            width=width,
            height=height,
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"})

        return fig

    #~~~~~~~ALIGNMENT COVERAGE METHOD AND HELPER~~~~~~~#
    def alignment_coverage(
            self,
            nbins: int = 500,
            color: str = 'rgba(70,130,180,0.70)',
            smooth_sigma: int = 1,
            width: int = None,
            height: int = 500,
            plot_title: str = "Coverage overview"):
        """
        Plot coverage over all the references
        * nbins
            Number of bins to divide the coverage into.
        * smooth_sigma
            sigma parameter for the Gaussian filter line smoothing
        * width
            With of the plotting area in pixel
        * height
            height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """
        # Verify that alignment information are available
        if not self.has_alignment:
            raise pycoQCError("No Alignment information available")
        self.logger.info("\t\tComputing plot")

        ref_offset_dict = self._ref_offset(self.ref_len_dict, "left", ret_type="dict")
        df = self.all_df[["ref_id", "ref_start", "ref_end"]].dropna()

        # Lay the references end to end and work in that coordinate.
        offsets = df["ref_id"].map(ref_offset_dict).to_numpy()
        starts = (offsets + df["ref_start"].to_numpy()).astype(np.float64)
        ends = (offsets + df["ref_end"].to_numpy()).astype(np.float64)

        # The dotted line has to be the mean of the area drawn beneath it,
        # so it measures the same thing: how much reference each alignment
        # spans. align_len counts aligned query bases instead, which
        # insertions inflate and deletions deflate, and the line then sat
        # off the curve by however indel rich the alignments were.
        mean_cov = round((ends - starts).sum() / self.total_ref_len, 2)

        y = self._binned_coverage(starts, ends, self.total_ref_len, nbins)

        # Time series smoothing
        if smooth_sigma:
            y = gaussian_filter1d(y, sigma=smooth_sigma)

        # Plot coverage area
        data1 = go.Scatter(
            x=list(range(nbins)),
            y=y,
            name="Mean coverage",
            hoveron="points",
            hoverinfo="y",
            fill='tozeroy',
            fillcolor=color,
            mode='none',
            showlegend=True,
            connectgaps=True)

        # Plot mean coverage
        data2 = go.Scatter(
            x=[0, nbins],
            y=[mean_cov, mean_cov],
            name=f"Overall coverage<br>{mean_cov}X",
            mode="lines",
            hoverinfo="skip",
            line={'color': 'gray', 'width': 2, 'dash': 'dot'})

        updatemenus = [
            dict(type="buttons", x=-0.2, y=0, xanchor='left', yanchor='bottom',  buttons=[
                dict(label="log", method='relayout', args=[{"yaxis": {"title": "Mean Coverage", "type": "log",
                                                                      "zeroline": False, "fixedrange": True}}]),
                dict(label="linear", method='relayout', args=[{"yaxis": {"title": "Mean Coverage", "type": "linear",
                                                                         "zeroline": False, "fixedrange": True}}])])]

        # Add chromosome shading and labels
        scale = nbins / self.total_ref_len
        shapes = []
        if len(self.ref_len_dict) <= COVERAGE_MAX_LABELLED_REFS:
            x_lab_coord = np.array(self._ref_offset(
                self.ref_len_dict, coordinates="middle", ret_type="list")) * scale
            x_lab = list(self.ref_len_dict.keys())
            x_title = None
            # Every other reference is shaded, the second, the fourth and so
            # on. The bounds were taken from the left edges alone, which have
            # no entry for where the last reference ends, so with an even
            # number of references the last one was never shaded.
            left = self._ref_offset(self.ref_len_dict, coordinates="left", ret_type="list")
            right = self._ref_offset(self.ref_len_dict, coordinates="right", ret_type="list")
            for i in range(1, len(self.ref_len_dict), 2):
                shapes.append(
                    go.layout.Shape(
                        type="rect", x0=left[i] * scale, x1=right[i] * scale,
                        y0=0, y1=1, yref="paper", opacity=0.5, layer="below", fillcolor="lightgrey",
                        line_width=0))
        else:
            x_lab_coord, x_lab = [], []
            x_title = "{:,} references laid end to end".format(len(self.ref_len_dict))

        # Tweak plot layout
        layout = go.Layout(
            width=width,
            height=height,
            plot_bgcolor="whitesmoke",
            updatemenus=updatemenus,
            shapes=shapes,
            hovermode="closest",
            legend={"x": -0.2, "y": 1, "xanchor": 'left', "yanchor": 'top'},
            xaxis={"zeroline": False, "showline": True, "ticktext": x_lab, "tickvals": x_lab_coord,
                   "tickangle": -45, "showgrid": False, "title": x_title},
            yaxis={"title": "Mean Coverage", "type": "log", "zeroline": False, "fixedrange": True},
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"})

        return go.Figure(data=[data1, data2], layout=layout)

    def _ref_offset(self, rlen, coordinates="left", ret_type="dict"):
        offset = [] if ret_type == "list" else OrderedDict()
        cumsum = 0
        for ref, rlen in rlen.items():
            # Define return val
            if coordinates == "left":
                v = cumsum
            elif coordinates == "middle":
                v = cumsum + rlen / 2
            else:
                v = cumsum + rlen
            # Add to appropriate collection
            if ret_type == "list":
                offset.append(v)
            else:
                offset[ref] = v
            cumsum += rlen
        return offset

    # ~~~~~~~BOXPLOT METHOD AND HELPER~~~~~~~#
    def length_boxplot_per_barcode(self, width: int = None, height: int = 500,
                                   plot_title: str = "Read length distribution per sample"):
        """
        Plot the distribution of read length per barcode, as one box per barcode
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """

        fig = self._per_barcode_boxplot(
            field_name="read_len",
            plot_title=plot_title,
            y_lab="Read length",
            line_color='rgb(102,168,255)',
            fill_color='rgb(204, 226, 255)',
            width=width,
            height=height,
            y_scale='log')

        return fig

    def gc_boxplot_per_barcode(self, width: int = None, height: int = 500,
                               plot_title: str = "%GC distribution per sample"):
        """
        Plot the distribution of %GC per barcode, as one box per barcode
        * width
            With of the plotting area in pixel
        * height
            Height of the plotting area in pixel
        * plot_title
            Title to display on top of the plot
        """

        # Check data contains %GC information
        if not self.has_gc_percent:
            raise pycoQCError("No %GC information available")

        fig = self._per_barcode_boxplot(
            field_name="%GC",
            plot_title=plot_title,
            y_lab="%GC",
            line_color='rgb(51, 204, 51)',
            fill_color='rgb(179, 255, 179)',
            width=width,
            height=height,
            y_scale='linear')

        return fig

    def _per_barcode_boxplot(self, field_name, plot_title, y_lab, line_color, fill_color, width, height, y_scale):
        """Private function generating boxplots plots for all per sample functions"""
        self._require_basecalls()

        # Verify that barcode information are available
        if not self.has_barcodes:
            raise pycoQCError("No barcode information available")

        self.logger.info("\t\tComputing plot")

        # Prepare all data
        prepared = [self._per_barcode_data(level, field_name) for level in LEVELS]
        _, dd1 = prepared[0]
        if not len(dd1["x"][0]):
            raise pycoQCError("No {} values to plot per barcode".format(field_name))

        # The boxes are drawn from statistics computed here. Handing plotly
        # every sampled read instead put the sample into the report four
        # times over, once for the plot and once per button: two boxplots
        # were 12 MB of a 100,000 read report, and the outlier points they
        # drew were most of what the browser had to render.
        data = [go.Box(x=dd1["x"][0], q1=dd1["q1"][0], median=dd1["median"][0], q3=dd1["q3"][0],
                       lowerfence=dd1["lowerfence"][0], upperfence=dd1["upperfence"][0],
                       mean=dd1["mean"][0], boxpoints=False, showlegend=False,
                       fillcolor=fill_color, line=dict(color=line_color))]

        # tweak plot layout
        layout = go.Layout(
            hovermode="closest",
            plot_bgcolor="whitesmoke",
            legend={"x": -0.2, "y": 1, "xanchor": 'left', "yanchor": 'top'},
            updatemenus=self._level_menu(prepared, "restyle", x=-0.2, y=0,
                                         xanchor='left', yanchor='bottom'),
            width=width,
            height=height,
            title={"text": plot_title, "xref": "paper", "x": 0.5, "xanchor": "center"},
            # In the order the rest of the report uses. Plotly's own
            # "category ascending" sorts as text, so barcode10 came between
            # barcode1 and barcode2, and Deepbinner's 10 between 1 and 2.
            xaxis={"categoryorder": "array",
                   "categoryarray": self._barcode_order(self.all_df["barcode"].unique())},
            yaxis={"title": y_lab, "type": y_scale, "zeroline": False, "showline": True, "fixedrange": True})

        return go.Figure(data=data, layout=layout)

    def _per_barcode_data(self, df_level, field_name="read_len"):
        """Box statistics per barcode for one read set.

        Quartiles are interpolated linearly, as plotly does, and the
        whiskers end at the furthest read within 1.5 times the interquartile
        range of the box, which is where plotly ends them too.
        """
        self.logger.debug("\t\tPreparing data for {} reads and {}".format(df_level, field_name))

        df, _ = self._level_sample(df_level)
        df = df[["barcode", field_name]].dropna()

        stats = {k: [] for k in ("x", "q1", "median", "q3", "lowerfence", "upperfence", "mean")}
        groups = {bc: values for bc, values in df.groupby("barcode", sort=False)[field_name]}
        for bc in self._barcode_order(groups):
            values = groups[bc].to_numpy(dtype=np.float64)
            q1, median, q3 = np.percentile(values, (25, 50, 75))
            iqr = q3 - q1
            stats["x"].append(bc)
            stats["q1"].append(q1)
            stats["median"].append(median)
            stats["q3"].append(q3)
            stats["lowerfence"].append(values[values >= q1 - 1.5 * iqr].min())
            stats["upperfence"].append(values[values <= q3 + 1.5 * iqr].max())
            stats["mean"].append(values.mean())

        # One trace, so each statistic is wrapped once for restyle.
        return self._level_label(df_level), {k: [v] for k, v in stats.items()}

    #~~~~~~~PRIVATE METHODS~~~~~~~#
    @staticmethod
    def _binned_coverage(starts, ends, total_len, nbins):
        """Mean depth per bin, counting every base each alignment covers

        An alignment used to be charged entirely to the bin its first base
        fell in. That is harmless when the bins are far wider than the reads,
        as on a human genome, but on the small references nanopore is
        routinely used for, a viral genome, a plasmid, an amplicon panel, it
        turned real coverage into a row of spikes with empty bins between
        them. The overall coverage line was right, so the picture looked
        plausible.

        Each alignment now contributes its bases to every bin it spans: the
        partial overlap at each end, and the whole bin width in between,
        accumulated through a difference array so the cost stays linear.

        * starts, ends
            Alignment bounds in the concatenated reference coordinate, with
            ends exclusive.
        * total_len
            Summed length of all references.
        * nbins
            Number of bins to divide them into.
        """
        width = total_len / nbins

        # Bins are indexed as position * nbins / total_len, the same mapping
        # the reference labels and the shading use, so the trace and the
        # annotations cannot drift apart. Integer floor division by a
        # truncated step did drift on short references.
        first = np.clip((starts / width).astype(np.int64), 0, nbins - 1)
        # The last bin an alignment touches is the one its final fraction
        # falls in. Deriving it from ends - 1 mixes a whole-base notion into
        # a continuous one: where a bin boundary is not an integer the tail
        # was charged to the bin before it, and once the bins are narrower
        # than a base, as on a reference shorter than nbins, each alignment
        # became a spike followed by empty bins.
        last = np.clip(np.ceil(ends / width).astype(np.int64) - 1,
                       0, nbins - 1)
        last = np.maximum(last, first)

        bases = np.zeros(nbins, dtype=np.float64)

        # Alignments confined to a single bin contribute their whole span.
        within = first == last
        if within.any():
            np.add.at(bases, first[within], (ends - starts)[within])

        spanning = ~within
        if spanning.any():
            head = first[spanning]
            tail = last[spanning]
            # The part of the first bin after the alignment starts, and the
            # part of the last bin before it ends.
            np.add.at(bases, head, width * (head + 1) - starts[spanning])
            np.add.at(bases, tail, ends[spanning] - width * tail)

            # Every bin strictly between them is covered end to end.
            interior = np.zeros(nbins + 1, dtype=np.float64)
            np.add.at(interior, head + 1, width)
            np.add.at(interior, tail, -width)
            bases += np.cumsum(interior)[:nbins]

        return bases / width

    @staticmethod
    def _smooth_with_gaps(values, sigma, min_weight=0.35):
        """Gaussian smoothing that tolerates empty bins

        Time bins with no reads in them are normal: a run pauses, output
        tails off, or the bins are simply finer than the data. Passing those
        gaps straight to gaussian_filter1d spreads the missing value over
        every neighbour, and a run with enough of them lost the whole series,
        leaving a plot with axes and no data.

        Each point is averaged over the neighbours that do have data, by
        filtering the values and a mask of where they are and dividing one by
        the other. Stretches with too little support nearby stay empty rather
        than being invented.

        * values
            One value per time bin, missing bins as nan.
        * sigma
            Width of the gaussian, in bins.
        * min_weight
            How much neighbouring data a point needs before it is drawn.
        """
        values = np.asarray(values, dtype=float)
        present = np.isfinite(values)
        if not present.any():
            return values

        filled = np.where(present, values, 0.0)
        smoothed = gaussian_filter1d(filled, sigma=sigma)
        weight = gaussian_filter1d(present.astype(float), sigma=sigma)

        out = np.full_like(values, np.nan)
        # A bin that was measured is always drawn. The weight threshold only
        # decides how far a value may be carried into surrounding gaps, and
        # on its own it discarded isolated measurements once sigma reached
        # about 1.5.
        usable = present | (weight > min_weight * weight.max())
        np.divide(smoothed, weight, out=out, where=usable)
        return out

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
                a=data, bins=pycoQC_plot._hist_bins(min_value, max_value, nbins, "log"))
        elif x_scale == "linear":
            count_y, bins = np.histogram(
                a=data, bins=pycoQC_plot._hist_bins(min_value, max_value, nbins))
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
