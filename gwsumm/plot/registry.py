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

"""Registry for GWSumm output plot types

All plot types should be registered for easy identification from the
configuration INI files
"""

import re
from six import string_types

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

_PLOTS = {}


def register_plot(plot, name=None, force=False):
    """Register a new summary `Plot` to the given ``name``

    Parameters
    ----------
    name : `str`
        unique descriptive name for this type of plot, must not
        contain any spaces, e.g. 'timeseries'
    plotclass : `type`
        defining Class for this plot type
    force : `bool`
        overwrite existing registration for this type

    Raises
    ------
    ValueError
        if name is already registered and ``force`` not given as `True`
    """
    if name is None:
        name = plot.type
    if name not in _PLOTS or force:
        _PLOTS[name] = plot
    else:
        raise ValueError("Plot '%s' has already been registered to the %s "
                         "class" % (name, plot.__name__))


def get_plot(name):
    """Query the registry for the plot class registered to the given
    name
    """
    if isinstance(name, string_types):
        name = re.sub(r'[\'\"]', '', name)
    try:
        return _PLOTS[name]
    except KeyError:
        raise ValueError("No TabPlot registered with name '%s'" % name)
