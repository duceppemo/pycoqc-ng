# -*- coding: utf-8 -*-
"""pycoQC_plot: every plot in the report, as one object.

Each family of plots lives in its own module under pycoQC.plots and is a
mixin here; PlotBase, last in line, holds the read sets and the helpers
they share. A plot is a public method returning a plotly figure, which is
what the report configuration names.
"""

from pycoQC.plots.alignment import Alignment
from pycoQC.plots.barcodes import Barcodes
from pycoQC.plots.base import (  # noqa: F401  re-exported for callers
    COVERAGE_MAX_LABELLED_REFS, LEVELS, PIE_MAX_BARCODES, SEED, PlotBase)
from pycoQC.plots.distributions import Distributions
from pycoQC.plots.over_time import OverTime
from pycoQC.plots.run_health import RunHealth
from pycoQC.plots.tables import SummaryTables

__all__ = ["pycoQC_plot", "COVERAGE_MAX_LABELLED_REFS", "LEVELS",
           "PIE_MAX_BARCODES", "SEED"]


class pycoQC_plot(SummaryTables, Distributions, OverTime, RunHealth,
                  Barcodes, Alignment, PlotBase):
    __doc__ = PlotBase.__init__.__doc__
