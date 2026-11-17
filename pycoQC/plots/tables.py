# -*- coding: utf-8 -*-
"""Summary tables and the statistics behind them and the JSON report.

Part of pycoQC_plot, which composes the plot families in pycoQC.plots into
one class. Methods here reach the read sets, the helpers and the other
families through self.
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library imports
from collections import OrderedDict
import datetime
import weakref

# Third party imports
import numpy as np
import plotly.graph_objs as go

# Local lib import
from pycoQC.common import pycoQCError
from pycoQC import __version__ as package_version
from pycoQC.plots.base import LEVELS

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class SummaryTables:
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

    def _metrics(self, df_level):
        """Every figure the summary tables and the JSON report give for one
        read set, computed once.

        Each table used to compute its own, through a method per figure,
        and the JSON report computed them all again: the barcode count was
        worked out five times per report and the reference length eight.
        """
        df = self._level_dfs(df_level)[0]
        # Each entry holds a weak reference to the frame it was computed
        # from and is used only while that same object is the read set, so
        # a replaced read set is computed again. Keying on id(df) could not
        # promise that, since an id is reused once its frame is freed, and a
        # strong reference kept a replaced frame in memory.
        cache = self.__dict__.setdefault("_metrics_cache", {})
        ref, metrics = cache.get(df_level, (None, None))
        if ref is None or ref() is not df:
            metrics = self._read_set_metrics(df)
            cache[df_level] = (weakref.ref(df), metrics)
        return metrics

    def _read_set_metrics(self, df):
        m = OrderedDict()
        m["run_duration"] = self._run_duration(df)
        m["active_channels"] = self._active_channels(df)
        m["runid_number"] = self._runid_number(df)
        m["barcodes_number"] = self._barcodes_number(df)

        if self.has_basecalls:
            m["reads"] = len(df)
            m["bases"] = int(df["read_len"].sum())
            m["N50"] = self._compute_N50(df["read_len"])
            m["median_read_len"] = self._median_or_nan(df["read_len"])
            m["median_qscore"] = self._median_or_nan(df["mean_qscore"])
            m["median_gc"] = (self._median_or_nan(df["%GC"])
                              if self.has_gc_percent else np.nan)

        if self.has_alignment:
            align_len = df["align_len"].dropna()
            m["aligned_reads"] = len(align_len)
            m["aligned_bases"] = int(align_len.sum())
            m["mean_coverage"] = self._alignment_mean_coverage(df)
            m["align_N50"] = self._compute_N50(align_len)
            m["align_median_len"] = self._median_or_nan(align_len)
            if self.has_identity_freq:
                m["median_identity"] = self._median_or_nan(df["identity_freq"])
                # Each error type per aligned base.
                for kind in ("insertion", "deletion", "mismatch"):
                    m[kind + "_rate"] = df[kind].dropna().sum() / m["aligned_bases"]
            else:
                m["median_identity"] = np.nan
        return m

    def summary_stats_dict(self):
        """
        Return a dictionnary containing exhaustive information about the run.
        """
        self.logger.info("\tCompute overall summary statistics")
        d = OrderedDict()

        d["pycoqc"] = OrderedDict()
        d["pycoqc"]["version"] = package_version
        d["pycoqc"]["date"] = datetime.datetime.now().strftime("%d/%m/%y")

        for level, lab in (("all", "All Reads"), ("pass", "Pass Reads")):
            d[lab] = self._compute_stats(level)
        return d

    def _compute_stats(self, df_level):
        df = self._level_dfs(df_level)[0]
        m = self._metrics(df_level)
        d = OrderedDict()
        # run information
        d["run"] = OrderedDict()
        for key in ("run_duration", "active_channels", "runid_number", "barcodes_number"):
            d["run"][key] = m[key]
        if self.has_end_reason:
            counts = df["end_reason"].fillna("unknown").value_counts()
            d["run"]["end_reasons"] = OrderedDict(
                (str(k), int(v)) for k, v in counts.items())
        # Everything under "basecall" describes a sequence, which a raw
        # signal run does not have yet. What it does have is how long each
        # read occupied its pore, which is reported in its place.
        if self.has_basecalls:
            d["basecall"] = OrderedDict()
            d["basecall"]["reads_number"] = m["reads"]
            d["basecall"]["bases_number"] = m["bases"]
            d["basecall"]["N50"] = m["N50"]
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
            d["alignment"]["reads_number"] = m["aligned_reads"]
            d["alignment"]["bases_number"] = m["aligned_bases"]
            d["alignment"]["mean_coverage"] = m["mean_coverage"]
            d["alignment"]["N50"] = m["align_N50"]
            d["alignment"]["len_percentiles"] = self._compute_percentiles(df["align_len"])
            x, y = self._compute_hist(data=df["align_len"], x_scale="log", smooth_sigma=2, nbins=100)
            d["alignment"]["len_hist"] = OrderedDict()
            d["alignment"]["len_hist"]["x"] = x
            d["alignment"]["len_hist"]["y"] = y

            if self.has_identity_freq:
                d["alignment"]["identity_freq_percentiles"] = self._compute_percentiles(df["identity_freq"])
                for key in ("insertion_rate", "deletion_rate", "mismatch_rate"):
                    d["alignment"][key] = m[key]
                x, y = self._compute_hist(data=df["identity_freq"], x_scale="linear", smooth_sigma=2, nbins=100)
                d["alignment"]["identity_freq_hist"] = OrderedDict()
                d["alignment"]["identity_freq_hist"]["x"] = x
                d["alignment"]["identity_freq_hist"]["y"] = y

        return d

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
        for level in LEVELS:
            m = self._metrics(level)
            data.append([self._level_label(level), m["run_duration"], m["active_channels"],
                         m["runid_number"], m["barcodes_number"]])

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
        for level in LEVELS:
            m = self._metrics(level)
            row = [self._level_label(level), m["reads"], m["bases"], m["N50"],
                   m["median_read_len"], m["median_qscore"]]
            if self.has_gc_percent:
                row.append(m["median_gc"])
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
        for level in LEVELS:
            m = self._metrics(level)
            data.append([self._level_label(level), m["aligned_reads"], m["aligned_bases"],
                         m["mean_coverage"], m["align_N50"], m["align_median_len"],
                         m["median_identity"]])

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

        # One pass over each frame, grouped by barcode. Selecting each
        # barcode with a mask scanned the whole run once per barcode, twice
        # over: a 96 barcode run read every row 192 times per read set.
        bc_list = self._barcode_order(df['barcode'].unique())
        totals = df.groupby('barcode', sort=False, observed=True)['read_len'].agg(['size', 'sum'])
        samples = dict(tuple(sample_df.groupby('barcode', sort=False, observed=True)))

        data = []
        for bc in bc_list:
            # A barcode with only a handful of reads can miss the sample
            # entirely, and its row then showed a real read count beside a
            # blank N50 and median. Such a barcode is small enough to
            # measure in full.
            sample1 = samples.get(bc)
            if sample1 is None or sample1.empty:
                sample1 = df[df['barcode'] == bc]
            row = [
                bc,
                int(totals.at[bc, 'size']),
                int(totals.at[bc, 'sum']),
                # Medians and N50 are distribution statistics, so the sample
                # estimates them directly.
                self._compute_N50(sample1['read_len']),
                self._median_or_nan(sample1['read_len']),
                self._median_or_nan(sample1['mean_qscore'])]
            if self.has_gc_percent:
                row.append(self._median_or_nan(sample1['%GC']))
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
