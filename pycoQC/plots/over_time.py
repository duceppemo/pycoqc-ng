# -*- coding: utf-8 -*-
"""How yield and the read distributions change over the run.

Part of pycoQC_plot, which composes the plot families in pycoQC.plots into
one class. Methods here reach the read sets, the helpers and the other
families through self.
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Third party imports
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
import plotly.graph_objs as go

# Local lib import
from pycoQC.common import pycoQCError
from pycoQC.plots.base import LEVELS

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class OverTime:
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
            updatemenus=self._level_menu(prepared, "update"),
            legend={"x": -0.05, "y": 1, "xanchor": 'right', "yanchor": 'top'},
            title=self._title(plot_title),
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
            top = y_cum_max if n % 2 == 0 else y_cum_max * 0.72
            lab_y.append([0, top])

        # make data dict
        data_dict = dict(
            x=[x, x]+lab_x,
            y=[y_cum, y]+lab_y,
            name=names,
            text=["", ""] + lab_text)

        # Headroom for the labels on top, which run to three lines. A sixth
        # of the range cut off their first line once the plot area had
        # given up some of its height to the row of buttons.
        layout_dict = {"yaxis.range": [0, (y_cum_max * 1.3) or 1]}
        return label, data_dict, layout_dict

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
            updatemenus=self._level_menu(prepared, "restyle"),
            legend={"x": -0.07, "y": 1, "xanchor": 'right', "yanchor": 'top'},
            title=self._title(plot_title),
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
