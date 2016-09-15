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

import os.path

from . import globalv

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

# set mode enum
SUMMARY_MODE_STATIC = 0
SUMMARY_MODE_EVENT = 1
SUMMARY_MODE_GPS = 2
SUMMARY_MODE_DAY = 10
SUMMARY_MODE_WEEK = 11
SUMMARY_MODE_MONTH = 12
SUMMARY_MODE_YEAR = 13

MODE_NAME = {
    SUMMARY_MODE_STATIC: 'STATIC',
    SUMMARY_MODE_EVENT: 'EVENT',
    SUMMARY_MODE_GPS: 'GPS',
    SUMMARY_MODE_DAY: 'DAY',
    SUMMARY_MODE_WEEK: 'WEEK',
    SUMMARY_MODE_MONTH: 'MONTH',
    SUMMARY_MODE_YEAR: 'YEAR',
}

MODE_ENUM = dict((val, key) for (key, val) in MODE_NAME.iteritems())


def get_mode(mode=None):
    """Return the current mode.
    """
    if mode is None:
        return globalv.MODE
    elif isinstance(mode, int):
        try:
            MODE_NAME[mode]
        except KeyError:
            ValueError("%s is not a valid mode" % mode)
        else:
            return mode
    else:
        try:
            return MODE_ENUM[mode]
        except KeyError:
            raise ValueError("%s is not a valid mode" % mode)


def set_mode(m):
    """Set the current mode.
    """
    if not isinstance(m, int):
        m = MODE_ENUM[str(m).upper()]
    globalv.MODE = m


def get_base(date, mode=None):
    """Determine the correct base attribute for the given date and mode.

    Parameters
    ----------
    date : :class:`datetime.datetime`
        formatted date
    mode : `int`, `str`
        enumerated interger code (or name) for the required mode

    Returns
    -------
    base : `str`
        the recommended base URL to have a correctly linked calendar
    """
    if mode is None:
        mode = get_mode()
    elif not isinstance(mode, int):
        mode = MODE_ENUM[str(mode).upper()]
    if mode == SUMMARY_MODE_DAY:
        return os.path.join('day', date.strftime('%Y%m%d'))
    elif mode == SUMMARY_MODE_WEEK:
        return os.path.join('week', date.strftime('%Y%m%d'))
    elif mode == SUMMARY_MODE_MONTH:
        return os.path.join('month', date.strftime('%Y%m'))
    elif mode == SUMMARY_MODE_YEAR:
        return os.path.join('year', date.strftime('%Y'))
    raise ValueError("Cannot format base for unknown mode %r" % mode)
