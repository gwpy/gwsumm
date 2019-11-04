# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013)
#
# This file is part of GWSumm.
#
# GWSumm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GWSumm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWSumm.  If not, see <http://www.gnu.org/licenses/>.

"""A `Plot` is a representation of an image to be included in the HTML
output a :doc:`tab </tabs>`.

For simple purposes, a `Plot` is just a reference to an existing image file
that can be imported into an HTML page via the ``<img>`` tag.

For more complicated purposes, a number of data plot classes are provided to
allow users to generate images on-the-fly.
The available classes are:

.. autosummary::
   :toctree: api

   TimeSeriesDataPlot
   SpectrogramDataPlot
   SegmentDataPlot
   StateVectorDataPlot
   SpectrumDataPlot
   TimeSeriesHistogramPlot
   TriggerTimeSeriesDataPlot
   TriggerHistogramPlot
   TriggerRateDataPlot
"""

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

from .registry import (
    register_plot,
    get_plot,
)
from .utils import (
    get_column_label,
    get_column_string,
    hash,
)
from .core import (
    format_label,
    SummaryPlot,
    DataPlot,
    BarPlot,
    PiePlot,
)
from .builtin import (
    undo_demodulation,
    TimeSeriesDataPlot,
    SpectrogramDataPlot,
    CoherenceSpectrogramDataPlot,
    SpectrumDataPlot,
    CoherenceSpectrumDataPlot,
    TimeSeriesHistogramPlot,
    TimeSeriesHistogram2dDataPlot,
    SpectralVarianceDataPlot,
    RayleighSpectrogramDataPlot,
    RayleighSpectrumDataPlot,
)
from .segments import (
    tint_hex,
    common_limits,
    SegmentDataPlot,
    StateVectorDataPlot,
    DutyDataPlot,
    ODCDataPlot,
    SegmentPiePlot,
    NetworkDutyPiePlot,
    SegmentBarPlot,
    SegmentHistogramPlot,
)
from .triggers import (
    TriggerPlotMixin,
    TriggerDataPlot,
    TriggerTimeSeriesDataPlot,
    TriggerHistogramPlot,
    TriggerRateDataPlot,
)
from .range import (
    _get_params,
    RangePlotMixin,
    RangeDataPlot,
    RangeDataHistogramPlot,
    RangeSpectrogramDataPlot,
    RangeSpectrumDataPlot,
    RangeCumulativeSpectrumDataPlot,
    SimpleTimeVolumeDataPlot,
    GWpyTimeVolumeDataPlot,
)
from .noisebudget import (
    NoiseBudgetPlot,
    RelativeNoiseBudgetPlot,
)
from .guardian import GuardianStatePlot
from .sei import SeiWatchDogPlot
