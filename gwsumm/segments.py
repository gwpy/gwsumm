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

"""Utilities for segment handling and display
"""

from .version import version as __version__
__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

try:
    from configparser import (ConfigParser, NoSectionError, NoOptionError)
except ImportError:
    from ConfigParser import (ConfigParser, NoSectionError, NoOptionError)
import warnings
import operator

from gwpy.segments import (DataQualityFlag, DataQualityDict, SegmentList)

import globalv
from .utils import *


def get_segments(flag, validity=None, config=ConfigParser(), cache=None,
                 query=True):
    """Retrieve the segments for a given flag

    Segments will be loaded from global memory if already defined,
    otherwise they will be loaded from the given
    :class:`~glue.lal.Cache`, or finally from the segment database

    Parameters
    ----------
    FIXME

    Returns
    -------
    FIXME
    """
    if isinstance(flag, basestring):
        flags = flag.split(',')
    else:
        flags = flag
    # check validity
    if not validity:
        start = config.get('general', 'gps-start-time')
        end = config.get('general', 'gps-end-time')
        validity = [(start, end)]
    # generate output object
    out = DataQualityDict()
    for f in flags:
        out[f] = DataQualityFlag(f, valid=validity, active=validity)
        globalv.SEGMENTS.setdefault(f, DataQualityFlag(f))
    # read segments from global memory and get the union of needed times
    try:
        old = reduce(operator.or_, (globalv.SEGMENTS.get(
                                        f, DataQualityFlag(f)).valid
                                    for f in flags))
    except TypeError:
        old = SegmentList()
    newsegs = validity - old
    # load new segments
    query &= abs(newsegs) != 0
    if query:
        if cache:
            raise NotImplementedError("Reading segments from cache has not "
                                      "been implemented yet.")
        else:
            # parse configuration for query
            kwargs = {}
            try:
                kwargs['url'] = config.get('segment-database', 'url')
            except (NoSectionError, NoOptionError):
                pass
            new = DataQualityDict.query(flags, newsegs, **kwargs)
            for f in flags:
                vprint("Downloaded %d new segments from database for %s.\n"
                       % (len(new[f].active), f))
        # record new segments
        globalv.SEGMENTS += new
    # return what was asked for
    for f in flags:
        out[f] &= globalv.SEGMENTS[f]
    if isinstance(flag, basestring):
        return out[flag]
    else:
        return out
