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

import re

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
    allflags = [f for cf in flags for f in re_flagdiv.split(cf) if f]

    # check validity
    if validity is None:
        start = config.get('general', 'gps-start-time')
        end = config.get('general', 'gps-end-time')
        validity = [(start, end)]
    elif isinstance(validity, DataQualityFlag):
        validity = validity.active

    # generate output object
    out = DataQualityDict()
    for f in flags:
        out[f] = DataQualityFlag(f, valid=validity, active=validity)
    for f in allflags:
        globalv.SEGMENTS.setdefault(f, DataQualityFlag(f))

    # read segments from global memory and get the union of needed times
    try:
        old = reduce(operator.and_, (globalv.SEGMENTS.get(
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
            new = DataQualityDict.query(allflags, newsegs, **kwargs)
            for f in new:
                vprint("    Downloaded %d new segments from database for %s.\n"
                       % (len(new[f].active), f))
        # record new segments
        globalv.SEGMENTS += new

    # return what was asked for
    for compound in flags:
        union, intersection, exclude, notequal = split_compound_flag(compound)
        for f in exclude:
            out[compound] -= globalv.SEGMENTS[f]
        for f in intersection:
            out[compound] &= globalv.SEGMENTS[f]
        for f in union:
            out[compound] |= globalv.SEGMENTS[f]
        for f in notequal:
            diff1 = out[compound] - globalv.SEGMENTS[f]
            diff2 = globalv.SEGMENTS[f] - out[compound]
            out[compound] &= (diff1 | diff2)
    if isinstance(flag, basestring):
        return out[flag]
    else:
        return out


def split_compound_flag(compound):
    """Parse the configuration for this state.

    Returns
    -------
    flags : `tuple`
        a 2-tuple containing lists of flags defining required ON
        and OFF segments respectively for this state
    """
    # find flags
    divs = re_flagdiv.findall(compound)
    keys = re_flagdiv.split(compound)
    # load flags and vetoes
    union = []
    intersection = []
    exclude = []
    notequal = []
    for i, key in enumerate(keys[::2]):
        if not key:
            continue
        # get veto bit
        if i != 0 and divs[i-1] == '!':
            exclude.append(key)
        elif i != 0 and divs[i-1] == '|':
            union.append(key)
        elif i != 0 and divs[i-1] == '!=':
            notequal.append(key)
        else:
            intersection.append(key)
    return union, intersection, exclude, notequal
