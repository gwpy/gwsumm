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

from __future__ import (division, print_function)
import sys

try:
    from configparser import (ConfigParser, NoSectionError, NoOptionError)
except ImportError:
    from ConfigParser import (ConfigParser, NoSectionError, NoOptionError)
import warnings
import operator

from gwpy.segments import (DataQualityFlag, DataQualityDict, SegmentList)

from . import (globalv, version)
from .utils import *

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


def get_segments(flag, validity=None, config=ConfigParser(), cache=None,
                 query=True, return_=True, segdb_error='raise', url=None):
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
    if isinstance(flag, (unicode, str)):
        flags = flag.split(',')
    else:
        flags = flag
    allflags = set([f for cf in flags for f in re_flagdiv.split(cf)[::2] if f])

    # check validity
    if validity is None:
        start = config.get('general', 'gps-start-time')
        end = config.get('general', 'gps-end-time')
        validity = [(start, end)]
    elif isinstance(validity, DataQualityFlag):
        validity = validity.active
    validity = SegmentList(validity)

    # generate output object
    out = DataQualityDict()
    for f in flags:
        out[f] = DataQualityFlag(f, known=validity, active=validity)
    for f in allflags:
        globalv.SEGMENTS.setdefault(f, DataQualityFlag(f))

    # read segments from global memory and get the union of needed times
    try:
        old = reduce(operator.and_, (globalv.SEGMENTS.get(
                                        f, DataQualityFlag(f)).known
                                    for f in flags))
    except TypeError:
        old = SegmentList()
    newsegs = validity - old
    # load new segments
    query &= abs(newsegs) != 0
    if query:
        if cache:
            try:
                new = DataQualityDict.read(cache, allflags)
            except Exception as e:
                if type(e) is not Exception:
                    raise
                if len(allflags) == 1:
                    f = list(allflags)[0]
                    new = DataQualityDict()
                    new[f] = DataQualityFlag.read(cache, f, coalesce=True)
            for flag in new:
                new[flag].known &= newsegs
                new[flag].active &= newsegs
        else:
            # parse configuration for query
            kwargs = {}
            if url is not None:
                kwargs['url'] = url
            else:
                try:
                    kwargs['url'] = config.get('segment-database', 'url')
                except (NoSectionError, NoOptionError):
                    pass
            try:
                if 'url' in kwargs and 'dqsegdb' in kwargs['url']:
                    new = DataQualityDict.query_dqsegdb(allflags, newsegs,
                                                        **kwargs)
                else:
                    new = DataQualityDict.query(allflags, newsegs, **kwargs)
            except Exception as e:
                # ignore error from SegDB
                if segdb_error in ['ignore', None]:
                    pass
                # convert to warning
                elif segdb_error in ['warn']:
                    print('%sWARNING: %sCaught %s: %s [gwsumm.segments]'
                          % (WARNC, ENDC, type(e).__name__, str(e)),
                          file=sys.stderr)
                    warnings.warn('%s: %s' % (type(e).__name__, str(e)))
                # otherwise raise as normal
                else:
                    raise
                new = DataQualityDict()
            for f in new:
                vprint("    Downloaded %d segments for %s (%.2f%% coverage).\n"
                       % (len(new[f].active), f,
                          float(abs(new[f].known))/float(abs(newsegs))*100))
        # record new segments
        globalv.SEGMENTS += new

    # return what was asked for
    if return_:
        for compound in flags:
            union, intersection, exclude, notequal = split_compound_flag(
                compound)
            for f in exclude:
                out[compound] -= globalv.SEGMENTS[f]
            for f in intersection:
                out[compound] &= globalv.SEGMENTS[f]
            for f in union:
                out[compound] |= globalv.SEGMENTS[f]
            for f in notequal:
                diff1 = out[compound] - globalv.SEGMENTS[f]
                diff2 = globalv.SEGMENTS[f] - out[compound]
                out[compound] = (diff1 | diff2)
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
