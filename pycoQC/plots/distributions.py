# -*- coding: utf-8 -*-
"""One and two dimensional distributions of the reads.

Part of pycoQC_plot, which composes the plot families in pycoQC.plots into
one class. Methods here reach the read sets, the helpers and the other
families through self.
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Third party imports
import numpy as np
from scipy.ndimage import gaussian_filter, gaussian_filter1d
import plotly.graph_objs as go

# Local lib import
from pycoQC.common import pycoQCError
from pycoQC.plots.base import LEVELS

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Distributions:
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
            updatemenus=self._level_menu(prepared, "update"),
            width=width,
            height=height,
            title=self._title(plot_title),
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

        # Percentile markers, each a vertical line labelled at its top. The
        # middle three sit close together on a narrow distribution, and on
        # the same line their labels printed over one another, so each label
        # is a step lower than the one before, over three steps.
        stat = np.percentile(data, [10, 25, 50, 75, 90])
        y_max = count_y.max()
        tops = [y_max * (1, 0.86, 0.72)[n % 3] for n in range(len(stat))]

        data_dict = dict(
            x=[count_x] + [[v, v] for v in stat],
            y=[count_y] + [[0, top] for top in tops],
            name=names,
            text=[""] + [["", "{}<br>{:,.2f}".format(n, v)] for n, v in zip(names[1:], stat)])

        # Make layout dict = Off set for labels on top
        layout_dict = {"yaxis.range": [0, (y_max + y_max / 6) or 1]}
        return label, data_dict, layout_dict

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
            updatemenus=self._level_menu(prepared, "restyle"),
            width=width,
            height=height,
            title=self._title(plot_title),
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
