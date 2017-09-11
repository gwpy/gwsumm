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
from enum import (Enum, unique)

from . import globalv

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


# -- operating Mode -----------------------------------------------------------

# https://docs.python.org/3/library/enum.html#orderedenum
class OrderedEnum(Enum):
    def __ge__(self, other):
        if self.__class__ is other.__class__:
            return self.value >= other.value
        return NotImplemented

    def __gt__(self, other):
        if self.__class__ is other.__class__:
            return self.value > other.value
        return NotImplemented

    def __le__(self, other):
        if self.__class__ is other.__class__:
            return self.value <= other.value
        return NotImplemented

    def __lt__(self, other):
        if self.__class__ is other.__class__:
            return self.value < other.value
        return NotImplemented


# set mode enum
@unique
class Mode(OrderedEnum):
    """Enumeration of valid processing 'modes'

    Each mode provides an association with a particular GPS interval
    """
    # no GPS associations
    static = 0
    # central GPS time (with duration)
    event = 1
    # arbitrary GPS [start, end) interval
    gps = 2
    # calendar epochs
    day = 10
    week = 11
    month = 12
    year = 13

    def dir_format(self):
        if self == Mode.day:
            return os.path.join('day', '%Y%m%d')
        elif self == Mode.week:
            return os.path.join('week', '%Y%m%d')
        elif self == Mode.month:
            return os.path.join('month', '%Y%m')
        elif self == Mode.year:
            return os.path.join('year', '%Y')
        raise ValueError("Cannot format base for Mode %s" % self)

    def is_calendar(self):
        if self >= Mode.day:
            return True
        return False


# -- Mode accessors -----------------------------------------------------------

def get_mode(m=None):
    """Return the enum for the given mode, defaults to the current mode.
    """
    if m is None:
        m = globalv.MODE
    if isinstance(m, (int, Enum)):
        return Mode(m)
    else:
        try:
            return Mode[str(m).lower()]
        except KeyError:
            raise ValueError("%s is not a valid Mode" % m)


def set_mode(m):
    """Set the current mode.
    """
    if isinstance(m, int):
        m = Mode(m)
    elif not isinstance(m, Mode):
        try:
            m = Mode[str(m).lower()]
        except KeyError:
            raise ValueError("%s is not a valid Mode" % m)
    globalv.MODE = m.value


# -- Mode utilities -----------------------------------------------------------

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
    mode = get_mode(mode)
    return date.strftime(mode.dir_format())
