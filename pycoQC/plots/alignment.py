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
from pycoQC.plots.base import COVERAGE_MAX_LABELLED_REFS, LEGEND_RIGHT, MENU_ROW

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~FUNCTIONS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def _bases_text(n):
    """A base count with its unit, short enough for the table's column."""
    n = float(n)
    for size, unit in ((1e9, "Gb"), (1e6, "Mb"), (1e3, "kb")):
        if abs(n) >= size:
            return "{:.2f} {}".format(n / size, unit)
    return "{:,.0f} b".format(n)


def _share_text(fraction):
    """A share as a percentage, to three significant figures.

    plotly's ".3p" wrote every decimal place out, so the share of a few
    hundred bases in half a gigabase came out as 0.0000012800%, wider than
    its column, and ran into the next one.
    """
    percent = float(fraction) * 100
    if not np.isfinite(percent):
        return ""
    if percent == 0:
        return "0%"
    if abs(percent) < 0.001:
        # Not "<": plotly sizes a table cell holding one, escaped or not,
        # as if it were a tag, and the row came out taller than the others.
        return "≤0.001%"
    return "{:.3g}%".format(percent) if abs(percent) < 1000 else "{:,.0f}%".format(percent)


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
        # Create empty multiplot figure. More than half the width for the
        # table: with two fifths, a narrow window squeezed its headers and
        # counts into one another, and the pie needs less.
        fig = make_subplots(rows=1, cols=2, column_widths=[0.55, 0.45], specs=[[{"type": "table"}, {"type": "pie"}]])

        # plot Table
        data = go.Table(
            columnwidth=[3, 2, 2],
            header={"values": list(df.columns), "align": "center", "fill_color": "grey", "font_size": 12,
                    "font_color": "white", "height": 40},
            # Counts with thousands separators and shares with a % sign, as
            # the report's other tables write them; "58202" and "74.41"
            # were the bare numbers.
            cells={"values": [df["Alignments"],
                              ["{:,}".format(int(n)) for n in df["Counts"]],
                              ["{:.2f}%".format(p) for p in df["Percents"]]],
                   "align": "center", "fill_color": "whitesmoke", "font_size": 12, "height": 30})
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
            textinfo='label+percent',
            # Each slice is labelled, so a legend only repeated the names,
            # in room the table needed.
            showlegend=False)
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
            height: int = 1000,
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
            columnwidth=[3, 2.4, 2, 2],
            header={
                "values": ["Bases", "Bases Count", "% Total", "% Aligned"],
                "align": "center", "fill_color": ["grey"], "font_size": 14, "font_color": "white", "height": 40},
            cells={
                "values": [df["Bases"],
                           [_bases_text(n) for n in df["Counts"]],
                           [_share_text(f) for f in df["Total_freq"]],
                           [_share_text(f) for f in df["Aligned_freq"]]],
                "align": "center", "fill_color": "whitesmoke", "font_size": 12, "height": 30})

        data2 = self._alignment_flow(df["Bases"].tolist(), df["Counts"].tolist(), list(colors))

        # Create multipanel figure
        # The table above the diagram, each the full width. Side by side, the
        # table's cells wrapped onto a second line in a narrow figure, and
        # the taller rows pushed the last of them out of the bottom of it.
        fig = make_subplots(rows=2, cols=1, row_heights=[0.46, 0.54], vertical_spacing=0.05,
                            specs=[[{"type": "table"}], [{"type": "sankey"}]])
        fig.add_trace(data1, row=1, col=1)
        fig.add_trace(data2, row=2, col=1)
        fig.update_layout(
            width=width,
            height=height,
            title=self._title(plot_title))

        return fig

    @staticmethod
    def _alignment_flow(labels, counts, colors):
        """The Sankey diagram of where the basecalled bases went.

        Laid out by plotly. Each node used to be given a fixed position,
        which assumed every one was drawn, but plotly leaves out a node
        nothing flows into, and the positions then went to the wrong nodes:
        with no unmapped reads, "Mapped reads" was put at the foot of the
        plot, half off it. Flows of nothing are left out. The nodes are lined
        up from the left where plotly can, so that unmapped reads sit beside
        mapped reads rather than in the last column with everything else
        nothing flows out of. Older plotly has no such setting.
        """
        # Node i flows from parents[i].
        parents = {1: 0, 2: 0, 3: 2, 4: 2, 5: 4, 6: 4, 7: 6, 8: 6, 9: 6}
        links = [(parents[i], i, counts[i]) for i in parents if counts[i] > 0]
        node = dict(pad=15, thickness=20, line={"width": 0},
                    label=list(labels), color=list(colors))
        if "align" in go.sankey.Node._valid_props:
            node["align"] = "left"
        return go.Sankey(
            arrangement="snap",
            node=node,
            link=go.sankey.Link(
                source=[s for s, _, _ in links],
                target=[t for _, t, _ in links],
                value=[v for _, _, v in links],
                hoverinfo="none"))

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
            legend=LEGEND_RIGHT,
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
