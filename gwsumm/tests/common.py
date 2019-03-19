# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013)
#
# This file is part of GWpy.
#
# GWpy is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# GWpy is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWpy.  If not, see <http://www.gnu.org/licenses/>.

"""Compatibility module
"""

from functools import wraps

from gwsumm import globalv


# -- test decorators ----------------------------------------------------------

def empty_globalv_CHANNELS(f):
    @wraps(f)
    def wrapped_f(*args, **kwargs):
        _channels = globalv.CHANNELS
        globalv.CHANNELS = type(globalv.CHANNELS)()
        try:
            return f(*args, **kwargs)
        finally:
            globalv.CHANNELS = _channels
    return wrapped_f
