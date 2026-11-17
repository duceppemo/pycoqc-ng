# -*- coding: utf-8 -*-
"""Alignment outcomes, rates and coverage.

Part of pycoQC_plot, which composes the plot families in pycoQC.plots into
one class. Methods here reach the read sets, the helpers and the other
families through self.
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library imports
from collections import OrderedDict, namedtuple

# Third party imports
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
import plotly.graph_objs as go
from plotly.subplots import make_subplots

# Local lib import
from pycoQC.common import pycoQCError
from pycoQC.plots.base import COVERAGE_MAX_LABELLED_REFS, MENU_ROW

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Alignment:
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
            title=self._title(plot_title))

        return fig

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
            title=self._title(plot_title))

        return fig

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
            dict(MENU_ROW, buttons=[
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
            title=self._title(plot_title))

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
