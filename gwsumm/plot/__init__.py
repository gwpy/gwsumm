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

"""Plotting utilities for GWSummary
"""

from gwpy.plotter import rcParams
rcParams.update({
    'figure.subplot.bottom': 0.17,
    'figure.subplot.left': 0.1,
    'figure.subplot.right': 0.9,
    'figure.subplot.top': 0.90,
})
rcParams['axes.labelsize'] = 20

from gwsumm import version

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

# import SuperClasses
from .core import *

# this import registers the builtins with their given names
from .builtin_plots import *

