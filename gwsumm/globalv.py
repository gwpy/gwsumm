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


"""Set of global memory variables for GWSumm package
"""

import time

from gwpy.time import to_gps
from gwpy.segments import DataQualityDict
from gwpy.detector import ChannelList

CHANNELS = ChannelList()
STATES = {}

DATA = {}
SPECTROGRAMS = {}
SPECTRUM = {}
COHERENCE_COMPONENTS = {}
COHERENCE_SPECTRUM = {}
SEGMENTS = DataQualityDict()
TRIGGERS = {}

VERBOSE = False
PROFILE = False
START = time.time()

# run time variables
MODE = 0
WRITTEN_PLOTS = []
NOW = int(to_gps('now'))
HTMLONLY = False

# comments
IFO = None
HTML_COMMENTS_NAME = None
