# -*- coding: utf-8 -*-
"""Read counts and distributions per barcode.

Part of pycoQC_plot, which composes the plot families in pycoQC.plots into
one class. Methods here reach the read sets, the helpers and the other
families through self.
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Third party imports
import numpy as np
import plotly.graph_objs as go

# Local lib import
from pycoQC.common import pycoQCError
from pycoQC.plots.base import LEVELS, PIE_MAX_BARCODES

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Barcodes:
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
