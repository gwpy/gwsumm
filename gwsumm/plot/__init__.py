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

from matplotlib import __version__ as mpl_version

from gwpy.plotter import rcParams

from .registry import *
from .utils import *
from .core import *
from .builtin import *
from .segments import *
from .triggers import *
from .range import *

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

rcParams.update({
    'grid.color': 'gray',
    'svg.fonttype': 'none',
})

# use new viridis colormap by default
if '1.5' <= mpl_version < '2.0':
    rcParams['image.cmap'] = 'viridis'
