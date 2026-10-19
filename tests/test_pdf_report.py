# -*- coding: utf-8 -*-
"""The static PDF report.

The HTML report embeds the plotly runtime and every data point, which runs
to megabytes. The PDF renders each plot once as an image instead.
"""

import pytest

from pycoQC.common import pycoQCError
from pycoQC.pycoQC import pycoQC

# Rendering needs a headless browser through kaleido, which is an optional
# extra, so skip rather than fail where it is unavailable.
pytest.importorskip("pypdf")
kaleido = pytest.importorskip("kaleido")


def can_render():
    import plotly.graph_objs as go
    try:
        go.Figure(go.Bar(x=[1], y=[1])).to_image(format="pdf", width=200,
                                                 height=200)
        return True
    except Exception:
        return False


pytestmark = pytest.mark.skipif(
    not can_render(), reason="kaleido cannot render here")


@pytest.fixture(scope="module")
def small_config(tmp_path_factory):
    """Two plots is enough for everything except the size comparison, and
    each rendered page costs a browser round trip."""
    import json
    path = tmp_path_factory.mktemp("cfg") / "config.json"
    path.write_text(json.dumps({
        "read_len_1D": {"plot_title": "Read length"},
        "read_qual_1D": {"plot_title": "Read quality"},
    }))
    return path


@pytest.fixture(scope="module")
def pdf_and_html(tmp_path_factory, guppy_summary):
    out = tmp_path_factory.mktemp("pdf")
    pdf = out / "report.pdf"
    html = out / "report.html"
    pycoQC(summary_file=str(guppy_summary), pdf_outfile=str(pdf),
           html_outfile=str(html), quiet=True)
    return pdf, html


class TestPdfReport:

    def test_pdf_is_written_with_one_page_per_plot(self, pdf_and_html):
        from pypdf import PdfReader
        pdf, _ = pdf_and_html
        assert pdf.exists()
        assert len(PdfReader(str(pdf)).pages) > 1

    def test_pdf_is_much_smaller_than_the_html(self, pdf_and_html):
        """The whole point of the PDF: something small enough to email."""
        pdf, html = pdf_and_html
        assert pdf.stat().st_size < html.stat().st_size / 5

    def test_pdf_alone_is_enough(self, tmp_path, guppy_summary, small_config):
        """A PDF on its own satisfies the "at least one output" rule."""
        pdf = tmp_path / "only.pdf"
        pycoQC(summary_file=str(guppy_summary), pdf_outfile=str(pdf),
               config_file=str(small_config), quiet=True)
        assert pdf.exists()

    def test_output_directory_is_created(self, tmp_path, guppy_summary,
                                         small_config):
        pdf = tmp_path / "nested" / "dir" / "report.pdf"
        pycoQC(summary_file=str(guppy_summary), pdf_outfile=str(pdf),
               config_file=str(small_config), quiet=True)
        assert pdf.exists()

    def test_config_selects_the_plots(self, tmp_path, guppy_summary,
                                      small_config):
        from pypdf import PdfReader
        pdf = tmp_path / "two.pdf"
        pycoQC(summary_file=str(guppy_summary), pdf_outfile=str(pdf),
               config_file=str(small_config), quiet=True)
        assert len(PdfReader(str(pdf)).pages) == 2

    def test_empty_config_is_reported(self, tmp_path, guppy_summary):
        import json
        config = tmp_path / "empty.json"
        config.write_text(json.dumps({}))
        with pytest.raises(pycoQCError, match="No plot"):
            pycoQC(summary_file=str(guppy_summary),
                   pdf_outfile=str(tmp_path / "x.pdf"),
                   config_file=str(config), quiet=True)
