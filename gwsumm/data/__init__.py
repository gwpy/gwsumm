# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013-2016)
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

"""Methods and classes for loading and pre-processing data

Each of the sub-modules are designed to read or create the data requested
only once, with the containers from the `globalv` module used as a storage
buffer for each unique data type
"""

# read TimeSeries data
from .timeseries import *

# generate Spectrograms and FrequencySeries
from .spectral import *

# generate coherence Spectrograms and Spectra
from .coherence import *

# generate TimeSeries of sensitive distance (range)
from .range import *
