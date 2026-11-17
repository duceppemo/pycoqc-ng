# -*- coding: utf-8 -*-
"""The state of the flow cell: end reasons, pores and channels.

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
class RunHealth:
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
            updatemenus=self._level_menu(prepared, "restyle"),
            # Under the plot, in a row: the stack fills the whole area, and
            # the reasons' names are long enough to run into the axis title
            # down the side.
            legend={"orientation": "h", "x": 0.5, "y": -0.18, "xanchor": "center",
                    "yanchor": "top", "traceorder": "normal"},
            margin={"b": 110},
            title=self._title(plot_title),
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
            and 1.5 hour ones by a tenth. The length is adjusted so that a
            whole number of equal intervals spans the run, and a run shorter
            than ten intervals is divided into ten.
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
        # Equal intervals that span the run exactly. Fixed length ones leave
        # the end of the run over, and wherever it goes it distorts the last
        # count: as an interval of its own, 0.9 hours of reads at the end of
        # a 23.4 hour run were drawn as a 1.5 hour interval and read as a
        # collapse; joined to the last one, the end of a 47.99 hour run made
        # an 8 hour interval out of a 4 hour one. A run is at least ten
        # intervals, so the length moves by at most a twentieth, and its
        # drift against the mux scan cycle adds up to half an interval over
        # the whole run.
        n = max(1, int(round(span / step)))
        intervals = (float(t.min()), span / n if span > 0 else step, n)

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
            updatemenus=self._level_menu(prepared, "restyle"),
            # Inside the plot, top right, where the count has fallen away
            # from the flow cell line by then. To the left of the axis, the
            # longer names ran into its title.
            legend={"x": 0.99, "y": 0.97, "xanchor": 'right', "yanchor": 'top',
                    "bgcolor": "rgba(255,255,255,0.7)"},
            title=self._title(plot_title),
            xaxis={"title": "Experiment time (h)", "zeroline": False, "showline": True},
            yaxis={"title": "Channels producing reads per {:.3g} h".format(intervals[1])
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

        # Each interval labelled by its centre. The intervals span the run
        # exactly, so only its very last read, at the closing edge, needs
        # clipping into the last one.
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
            updatemenus=self._level_menu(prepared, "restyle"),
            title=self._title(plot_title),
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
