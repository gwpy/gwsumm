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

"""Job modes
"""

from . import version

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

# set mode enum
SUMMARY_MODE_DAY = 0
SUMMARY_MODE_MONTH = 1
SUMMARY_MODE_YEAR = 2
SUMMARY_MODE_WEEK = 3
SUMMARY_MODE_GPS = 4

MODE_NAME = {SUMMARY_MODE_GPS: 'GPS',
             SUMMARY_MODE_DAY: 'DAY',
             SUMMARY_MODE_WEEK: 'WEEK',
             SUMMARY_MODE_MONTH: 'MONTH',
             SUMMARY_MODE_YEAR: 'YEAR'}

MODE_ENUM = dict((val, key) for (key, val) in MODE_NAME.iteritems())

MODE = None
