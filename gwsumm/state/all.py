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

"""Definition of the 'All' state.

This is a special `SummaryState` that has valid and active segments spanning
the full analysis interval.
"""

from gwpy.segments import (Segment, SegmentList)

from ..globalv import NOW
from .core import SummaryState
from .registry import register_state

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

ALLSTATE = 'all'


def generate_all_state(start, end, register=True, **kwargs):
    """Build a new `SummaryState` for the given [start, end) interval.

    Parameters
    ----------
    start : `~gwpy.time.LIGOTimeGPS`, float
        the GPS start time of the current analysis
    end : `~gwpy.time.LIGOTimeGPS`, float
        the GPS end time of the current analysis
    register : `bool`, optional
        should the new `SummaryState` be registered, default `True`
    **kwargs
        other keyword arguments passed to the `SummaryState` constructor

    Returns
    -------
    allstate : `SummaryState`
        the newly created 'All' `SummaryState`
    """
    now = min(end, NOW)
    all_ = SummaryState(ALLSTATE,
                        known=SegmentList([Segment(start, end)]),
                        active=SegmentList([Segment(start, now)]),
                        **kwargs)
    all_.ready = True
    if register:
        register_state(all_)
    return all_
