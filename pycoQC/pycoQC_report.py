#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#~~~~~~~~~~~~~~IMPORTS~~~~~~~~~~~~~~#
# Standard library imports
import base64
import getpass
import io
import json
import logging
import math
import socket
from collections import OrderedDict
from importlib.resources import files
import datetime
import os

# Third party imports
import numpy as np
import plotly.graph_objs as go
import plotly.offline as py
import jinja2
from markupsafe import Markup

# Local imports
from pycoQC.common import get_logger, mkbasedir, pycoQCError
from pycoQC.pycoQC_parse import pycoQC_parse
from pycoQC.pycoQC_plot import pycoQC_plot
from pycoQC import __version__ as package_version
from pycoQC import __package_name__ as package_name

def _logo_data_uri():
    """The logo as a data URI, so the report stays one self-contained file

    Linking it from the repository would have made an offline report show a
    broken image, which defeats the point of embedding every plot.
    """
    try:
        svg = files("pycoQC").joinpath("templates/logo.svg").read_bytes()
    except (FileNotFoundError, OSError):
        return ""
    return "data:image/svg+xml;base64," + base64.b64encode(svg).decode("ascii")


def _spectre_css():
    """The vendored Spectre stylesheets, concatenated"""
    css = []
    for name in ("spectre.min.css", "spectre-exp.min.css",
                 "spectre-icons.min.css"):
        try:
            css.append(files("pycoQC").joinpath(
                "templates/vendor/" + name).read_text(encoding="utf-8"))
        except (FileNotFoundError, OSError):
            continue
    return "\n".join(css)


def _logo_png_uri():
    """The logo as a PNG data URI, for the PDF cover

    plotly's layout images go through the rendering engine, which does not
    take the SVG the HTML report embeds.
    """
    try:
        png = files("pycoQC").joinpath("templates/logo.png").read_bytes()
    except (FileNotFoundError, OSError):
        return ""
    return "data:image/png;base64," + base64.b64encode(png).decode("ascii")


def _current_user():
    """Who ran this, for the report's provenance line"""
    try:
        return getpass.getuser()
    except Exception:
        # getuser raises when there is no password entry, as in some
        # containers, and the report should not fail over a byline.
        return os.environ.get("USER") or os.environ.get("USERNAME") or "unknown"


# kaleido starts a headless browser for each image, and starting one can fail
# for reasons that do not recur on the next try (kaleido 1.5.1 often ran its
# script in the wrong page, failing with "kaleido_scopes is not defined").
# Each page gets a few attempts before the error is taken as a broken
# rendering engine.
PDF_RENDER_ATTEMPTS = 3


def _render_pdf_page(fig, width, height, logger=None):
    """One figure as PDF bytes, retrying kaleido's intermittent failures"""
    for attempt in range(1, PDF_RENDER_ATTEMPTS + 1):
        try:
            return fig.to_image(format="pdf", width=width, height=height)
        except Exception as E:
            if attempt == PDF_RENDER_ATTEMPTS:
                raise
            if logger is not None:
                logger.debug("\t\tRendering failed ({}), trying again".format(E))


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~MAIN CLASS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
def _json_safe(obj):
    """Plain Python values JSON can hold, with None for non finite floats.

    numpy scalars are converted first. np.float64 passes for a float, but
    np.float32 and the integer types do not, and the median of a float32
    column made json.dump fail on the whole report.
    """
    if isinstance(obj, np.generic):
        obj = obj.item()
    if isinstance(obj, float):
        return obj if math.isfinite(obj) else None
    if isinstance(obj, dict):
        return OrderedDict((k, _json_safe(v)) for k, v in obj.items())
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj


class pycoQC_report ():

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~INIT METHOD~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    def __init__(self,
                 parser: pycoQC_parse,
                 plotter: pycoQC_plot,
                 verbose: bool = False,
                 quiet: bool = False):
        """
        * parser
            A pycoQC_parse object
        * plotter
            A pycoQC_plot object
        * verbose
            Increase verbosity
        * quiet
            Reduce verbosity
        """
        # Set logging level
        self.logger = get_logger(name=__name__, verbose=verbose, quiet=quiet)

        # Check that parser is a valid instance of pycoQC_parse
        if not isinstance(parser, pycoQC_parse):
            raise pycoQCError("{} is not a valid pycoQC_parse object".format(parser))
        self.parser = parser

        # Check that plotter is a valid instance of pycoQC_plot
        if not isinstance(plotter, pycoQC_plot):
            raise pycoQCError("{} is not a valid pycoQC_plot object".format(plotter))
        self.plotter = plotter

    def __repr__(self):
        return "[{}]\n".format(self.__class__.__name__)

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PUBLIC METHODS~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    def pdf_report(self,
                   outfile: str,
                   config_file: str = "",
                   report_title: str = "PycoQC report",
                   skip_coverage_plot: bool = False,
                   width: int = 1000,
                   height: int = 600):
        """Write the same plots as the HTML report, one per page, as a PDF

        The HTML report embeds the plotly runtime and every point of every
        plot, which runs to several megabytes. This renders each plot once as
        a static image, giving a file small enough to attach to an email or
        drop into a lab notebook, at the cost of the interactivity.

        * outfile
            Path of the PDF to write.
        * config_file
            Same configuration as the HTML report, choosing the plots.
        * report_title
            Title printed on each page.
        * skip_coverage_plot
            Skip the coverage plot, as for the HTML report.
        * width, height
            Page size in pixels.
        """
        self.logger.info("Generating PDF report")

        # Imported here so the dependency is only needed when a PDF is asked
        # for; the HTML report does not require a rendering engine.
        try:
            from pypdf import PdfReader, PdfWriter
        except ImportError:
            raise pycoQCError(
                "Writing a PDF report needs the pdf extra: "
                "pip install 'pycoqc-ng[pdf]'")

        # kaleido and the browser driver beneath it log their internals,
        # including temporary directory cleanup, to stderr. That is noise in a
        # command line tool, so quiet them unless debug output was asked for.
        if not self.logger.isEnabledFor(logging.DEBUG):
            for name in ("kaleido", "choreographer", "logistro"):
                logging.getLogger(name).setLevel(logging.WARNING)

        self.logger.info("\tParsing config file")
        config_dict = self._get_config(config_file)

        writer = PdfWriter()
        n_plots = 0

        # A cover carrying the same title, run details and provenance as the
        # HTML header, so a PDF filed as a QA record stands on its own.
        cover = self._pdf_cover(report_title, width, height)
        if cover is not None:
            try:
                writer.append(PdfReader(io.BytesIO(
                    _render_pdf_page(cover, width, height, self.logger))))
            except Exception as E:
                self.logger.info("\t\tCould not render the cover: {}".format(E))
        for method_name, plot_title, fig in self._iter_figures(config_dict, skip_coverage_plot):
            # A plot that asked for more room than a page gets it: the
            # per barcode table and pie size themselves to the number of
            # barcodes, and forcing them back to the page height cut a
            # 96 barcode run off after the twelfth row.
            page_height = max(height, int(fig["layout"]["height"] or 0))

            fig.update_layout(title=dict(text=plot_title),
                              width=width, height=page_height)
            try:
                page = _render_pdf_page(fig, width, page_height, self.logger)
            except Exception as E:
                # Unlike a plot without data, a rendering engine that does
                # not work is fatal: every page would fail the same way.
                raise pycoQCError(
                    "Could not render {} to PDF. Writing a PDF report "
                    "needs the pdf extra: pip install 'pycoqc-ng[pdf]'. "
                    "Underlying error: {}".format(method_name, E))
            writer.append(PdfReader(io.BytesIO(page)))
            n_plots += 1

        if not n_plots:
            raise pycoQCError("No plot could be rendered for the PDF report")

        self.logger.info("\tWriting to PDF file")
        mkbasedir(outfile, exist_ok=True)
        with open(outfile, "wb") as fh:
            writer.write(fh)
        self.logger.info("\t\t{} pages written".format(n_plots))

    def _pdf_cover(self, report_title, width, height):
        """Build the PDF's first page: title, run details and provenance"""
        generated = datetime.datetime.now().astimezone()
        rows = []
        for key, value in (getattr(self.parser, "run_metadata", None) or {}).items():
            rows.append((key, str(value)))
        if rows:
            rows.append(("", ""))
        rows += [
            ("Generated", generated.strftime("%Y-%m-%d %H:%M %Z").strip()),
            ("Generated by", _current_user()),
            ("On host", socket.gethostname()),
            ("Software", "{} {}".format(package_name, package_version)),
        ]

        # Margins are zeroed and everything placed by domain, so paper
        # coordinates span the whole page. With a top margin, paper y=1 is
        # the top of the plotting area rather than of the figure, and the
        # logo landed on top of the table.
        table_top = 0.58
        fig = go.Figure(go.Table(
            domain=dict(x=[0.09, 0.91], y=[0.04, table_top]),
            columnwidth=[34, 66],
            header=dict(values=["", ""], height=1,
                        fill_color="rgba(0,0,0,0)", line_color="rgba(0,0,0,0)"),
            cells=dict(
                values=[[r[0] for r in rows], [r[1] for r in rows]],
                align="left", height=27,
                font=dict(size=12, color=["#64748b", "#1f2933"]),
                fill_color="rgba(0,0,0,0)",
                line_color="#e7e9ed")))

        fig.update_layout(
            margin=dict(l=0, r=0, t=0, b=0),
            paper_bgcolor="white", plot_bgcolor="white",
            annotations=[dict(
                text=report_title, x=0.5, y=table_top + 0.10,
                xref="paper", yref="paper", xanchor="center", yanchor="middle",
                showarrow=False, font=dict(size=23, color="#1f2933"))])

        logo = _logo_png_uri()
        if logo:
            fig.add_layout_image(dict(
                source=logo, xref="paper", yref="paper",
                x=0.5, y=0.96, sizex=0.46, sizey=0.18,
                xanchor="center", yanchor="top", layer="above"))
        return fig

    def html_report(self,
                    outfile: str,
                    config_file: str = "",
                    template_file: str = "",
                    report_title: str = "PycoQC report",
                    skip_coverage_plot: bool = False):
        """"""
        self.logger.info("Generating HTML report")

        # Parse configuration file
        self.logger.info("\tParsing html config file")
        config_dict = self._get_config(config_file)
        self.logger.debug(config_dict)

        # Loop over configuration file and run the pycoQC functions defined
        plots = list()
        titles = list()
        for method_name, plot_title, fig in self._iter_figures(
                config_dict, skip_coverage_plot, blank_title=True):
            plots.append(Markup(py.plot(
                fig,
                output_type='div',
                include_plotlyjs=False,
                image_width='',
                image_height='',
                auto_open=False)))
            titles.append(plot_title)

        # Load HTML template for Jinja
        self.logger.info("\tLoading HTML template")
        template = self._get_jinja_template(template_file)

        # Provenance for the report itself: who generated it, when and with
        # what, which is what a QA record needs alongside the run details.
        generated = datetime.datetime.now().astimezone()
        report_meta = OrderedDict((
            ("Generated", generated.strftime("%Y-%m-%d %H:%M %Z").strip()),
            ("Generated by", _current_user()),
            ("On host", socket.gethostname()),
            ("Software", "{} {}".format(package_name, package_version)),
        ))
        report_subtitle = "Generated on {} with {} {}".format(
            generated.strftime("%Y-%m-%d"), package_name, package_version)

        # Run details collected while parsing, before the unused columns were
        # dropped. Empty for a summary that carries none.
        run_meta = getattr(self.parser, "run_metadata", None) or OrderedDict()

        # Define source files list. The paths are escaped: a directory
        # name holding & or < broke the page's markup around it.
        source_files = OrderedDict()
        for files_list, name in (
                (self.parser.summary_files_list, "summary"),
                (self.parser.barcode_files_list, "barcode"),
                (self.parser.bam_file_list, "bam")):
            if files_list:
                source_files[name] = [os.path.abspath(f) for f in files_list]
        # Kept as ready made HTML for templates written against it.
        src_files = Markup("").join(
            Markup("<h4>Source {} files</h4><ul>{}</ul>").format(
                name, Markup("").join(Markup("<li>{}</li>").format(f) for f in paths))
            for name, paths in source_files.items())

        # Render plots. Everything generated here is marked safe; what came
        # from the input, the run details and the title, is escaped by the
        # default template.
        self.logger.info("\tRendering plots in d3js")
        rendering = template.render(
            plots=plots,
            titles=titles,
            plotlyjs=Markup(py.get_plotlyjs()),
            report_title=report_title,
            report_subtitle=report_subtitle,
            report_meta=report_meta,
            run_meta=run_meta,
            logo_url=_logo_data_uri(),
            spectre_css=Markup(_spectre_css()),
            source_files=source_files,
            src_files=src_files)

        # Write to HTML file, in the encoding the page declares rather than
        # whatever the locale's default happens to be.
        self.logger.info("\tWriting to HTML file")
        mkbasedir(outfile, exist_ok=True)
        with open(outfile, "w", encoding="utf-8") as fp:
            fp.write(rendering)

    def json_report(self,
                    outfile: str):
        """"""
        self.logger.info("Generating JSON report")
        self.logger.info("\tRunning summary_stats_dict method")
        res_dict = self.plotter.summary_stats_dict()

        self.logger.info("\tWriting to JSON file")
        mkbasedir(outfile, exist_ok=True)
        with open(outfile, "w", encoding="utf-8") as fp:
            # A run with no alignments, or a barcode with no passing reads,
            # leaves NaN in the stats. json.dump writes that as a bare NaN,
            # which Python reads back but which is not JSON: jq, R and
            # JavaScript all reject the file. null says the same thing in a
            # form every reader accepts.
            json.dump(_json_safe(res_dict), fp, indent=2, allow_nan=False)

    #~~~~~~~~~~~~~~PRIVATE FUNCTION~~~~~~~~~~~~~~#

    def _iter_figures(self, config_dict, skip_coverage_plot=False, blank_title=False):
        """Run each plot the configuration names, yielding the ones drawn.

        Shared by the HTML and PDF reports, which used to carry a copy each
        of this loop and of its error handling.

        Yields (method name, title, figure). A plot whose data the run does
        not have raises pycoQCError and is skipped, as is a name that is not
        a plotting method. Any other failure is reported as a warning and
        the rest of the report is still written: one plot tripping over an
        unusual run is no reason to lose all the others.

        * config_dict
            Plot names mapped to the arguments for each.
        * skip_coverage_plot
            Leave out alignment_coverage.
        * blank_title
            Draw the figure without its title, for a page that prints the
            title itself.
        """
        for method_name, method_args in config_dict.items():
            if skip_coverage_plot and method_name == "alignment_coverage":
                self.logger.info("\tSkipping method {}".format(method_name))
                continue

            # Only the public plotting methods. A config naming a private
            # helper, or summary_stats_dict, called it and then failed on
            # its return value.
            method = getattr(self.plotter, method_name, None)
            if method_name.startswith("_") or not callable(method):
                self.logger.info("\t\t{} is not a valid plotting method".format(method_name))
                continue

            self.logger.info("\tRunning method {}".format(method_name))
            self.logger.debug("\t{} ({})".format(method_name, method_args))

            # Store plot title for HTML title and remove from data
            # passed to plotly. A user written config entry need not
            # carry one, and requiring it raised an uncaught KeyError.
            method_args = dict(method_args or {})
            plot_title = method_args.pop("plot_title", method_name)
            if blank_title:
                method_args["plot_title"] = ""

            try:
                fig = method(**method_args)
            except pycoQCError as E:
                self.logger.info("\t\t{}".format(E))
                continue
            except Exception as E:
                self.logger.warning(
                    "WARNING: {} failed with {}: {}".format(
                        method_name, type(E).__name__, E))
                self.logger.debug("", exc_info=True)
                continue

            if not isinstance(fig, go.Figure):
                self.logger.info("\t\t{} is not a valid plotting method".format(method_name))
                continue
            # Do not display %GC graph if no data available
            if not fig["data"]:
                continue

            yield method_name, plot_title, fig

    def _get_config(self, config_file=None):
        """"""
        # First, try to read provided configuration file if given
        if config_file:
            self.logger.debug("\tTry to read provided config file")
            try:
                with open(config_file, 'r') as cf:
                    return json.load(cf)
            except (FileNotFoundError, IOError, json.JSONDecodeError) as E:
                # Falling back silently meant a typo in --config_file
                # produced a different report with no sign anything was
                # wrong. The user named this file, so say it was not used.
                self.logger.warning(
                    "WARNING: could not read the configuration file {} ({}). "
                    "Falling back to the default configuration.".format(
                        config_file, E))

        # Then the per user configuration the command line help has always
        # promised. Nothing looked in either location, so a file put there
        # was ignored and the report came out with the default plots.
        for candidate in (os.path.expanduser("~/.pycoQC"),
                          os.path.expanduser("~/.config/pycoQC/config")):
            if os.path.isfile(candidate):
                self.logger.debug("\tTry to read {}".format(candidate))
                try:
                    with open(candidate) as cf:
                        config = json.load(cf)
                    self.logger.info(
                        "\tUsing the configuration in {}".format(candidate))
                    return config
                except (IOError, json.JSONDecodeError) as E:
                    self.logger.warning(
                        "WARNING: could not read the configuration file {} "
                        "({}). Falling back to the default "
                        "configuration.".format(candidate, E))

        # Last use the default hardcoded config_dict. A run that has not
        # been basecalled has no read lengths or quality scores, so most of
        # the usual plots would be skipped one after another and the report
        # would be mostly empty; it gets its own set describing the run
        # itself instead.
        if self.plotter.has_basecalls:
            name = "templates/pycoQC_config.json"
        else:
            self.logger.debug("\tNo basecalls: using the raw signal configuration")
            name = "templates/pycoQC_config_signal.json"
        self.logger.debug("\tRead default configuration file")
        config_json = files("pycoQC").joinpath(name)
        return json.loads(config_json.read_text())

    def _get_jinja_template(self, template_file=None):
        """"""
        # First, try to read provided configuration file if given
        if template_file:
            self.logger.debug("\tTry to load provided jinja template file")
            try:
                with open(template_file) as fp:
                    template = jinja2.Template(fp.read())
                    return template
            except (FileNotFoundError, IOError, jinja2.exceptions.TemplateNotFound,
                    jinja2.exceptions.TemplateSyntaxError) as E:
                self.logger.warning(
                    "WARNING: could not read the template file {} ({}). "
                    "Falling back to the default template.".format(
                        template_file, E))

        # Last use the default harcoded config_dict
        self.logger.debug("\tRead default jinja template")
        # select_autoescape matches on the end of the name, and the
        # template is called spectre.html.j2, so ["html"] alone left it
        # unescaped: a sample name holding & or < went into the page as
        # markup.
        env = jinja2.Environment(
            loader=jinja2.PackageLoader('pycoQC', 'templates'),
            autoescape=jinja2.select_autoescape(["html", "html.j2"]))
        template = env.get_template('spectre.html.j2')
        return template
