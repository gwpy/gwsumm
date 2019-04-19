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
import operator
import sys
import warnings
from collections import OrderedDict
from configparser import (
    DEFAULTSECT,
    ConfigParser,
    NoSectionError,
    NoOptionError,
)

from six import string_types
from six.moves import reduce

from astropy.io.registry import IORegistryError

from gwpy.segments import (DataQualityFlag, DataQualityDict,
                           SegmentList, Segment)

from . import globalv
from .utils import (
    re_flagdiv,
    vprint,
    WARNC,
    ENDC,
)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

SEGDB_URLS = [
    'https://segdb.ligo.caltech.edu',
    'https://metaserver.phy.syr.edu',
    'https://geosegdb.atlas.aei.uni-hannover.de',
    'http://10.20.50.30'  # geosegdb internal
]


def get_segments(flag, validity=None, config=ConfigParser(), cache=None,
                 query=True, return_=True, coalesce=True, padding=None,
                 segdb_error='raise', url=None):
    """Retrieve the segments for a given flag

    Segments will be loaded from global memory if already defined,
    otherwise they will be loaded from the given
    :class:`~glue.lal.Cache`, or finally from the segment database

    Parameters
    ----------
    flag : `str`, `list`
        either the name of one flag, or a list of names

    validity : `~gwpy.segments.SegmentList`
        the segments over which to search for other segments

    query : `bool`, optional, default: `True`
        actually execute a read/query operation (if needed), otherwise
        just retrieve segments that have already been cached

    config : `~configparser.ConfigParser`, optional
        the configuration for your analysis, if you have one. If
        present the ``[segment-database]`` section will be queried
        for the following options

        - ``gps-start-time``, and ``gps-end-time``, if ``validity`` is
          not given
        - ``url`` (the remote hostname for the segment database) if
          the ``url`` keyword is not given

    cache : :class:`glue.lal.Cache`, optional
        a cache of files from which to read segments, otherwise segments
        will be downloaded from the segment database

    coalesce : `bool`, optional, default: `True`
        coalesce all segmentlists before returning, otherwise just return
        segments as they were downloaded/read

    padding : `tuple`, or `dict` of `tuples`, optional
        `(start, end)` padding with which to pad segments that are
        downloaded/read

    segdb_error : `str`, optional, default: ``'raise'``
        how to handle errors returned from the segment database, one of

        - ``'raise'`` (default) : raise the exception as normal
        - ``'warn'`` : print the exception as a warning, but return no
          segments
        - ``'ignore'`` : silently ignore the error and return no segments

    url : `str`, optional
        the remote hostname for the target segment database

    return_ : `bool`, optional, default: `True`
        internal flag to enable (True) or disable (False) actually returning
        anything. This is useful if you want to download/read segments now
        but not use them until later (e.g. plotting)

    Returns
    -------
    flag : `~gwpy.segments.DataQualityFlag`
        the flag object representing segments for the given single flag, OR

    flagdict : `~gwpy.segments.DataQualityDict`
        the dict of `~gwpy.segments.DataQualityFlag` objects for multiple
        flags, if ``flag`` is given as a `list`, OR

    None
       if ``return_=False``
    """
    if isinstance(flag, string_types):
        flags = flag.split(',')
    else:
        flags = flag
    allflags = set([f for cf in flags for f in
                    re_flagdiv.split(str(cf))[::2] if f])

    if padding is None and isinstance(flag, DataQualityFlag):
        padding = {flag: flag.padding}
    elif padding is None:
        padding = dict((flag,
                        isinstance(flag, DataQualityFlag) and
                        flag.padding or None) for flag in flags)

    # check validity
    if validity is None:
        start = config.get(DEFAULTSECT, 'gps-start-time')
        end = config.get(DEFAULTSECT, 'gps-end-time')
        span = SegmentList([Segment(start, end)])
    elif isinstance(validity, DataQualityFlag):
        validity = validity.active
        try:
            span = SegmentList([validity.extent()])
        except ValueError:
            span = SegmentList()
    else:
        try:
            span = SegmentList([SegmentList(validity).extent()])
        except ValueError:
            span = SegmentList()
    validity = SegmentList(validity)

    # generate output object
    out = DataQualityDict()
    for f in flags:
        out[f] = DataQualityFlag(f, known=validity, active=validity)
    for f in allflags:
        globalv.SEGMENTS.setdefault(f, DataQualityFlag(f))

    # read segments from global memory and get the union of needed times
    try:
        old = reduce(
            operator.and_,
            (globalv.SEGMENTS.get(f, DataQualityFlag(f)).known for f in flags))
    except TypeError:
        old = SegmentList()
    newsegs = validity - old
    # load new segments
    query &= abs(newsegs) != 0
    query &= len(allflags) > 0
    if cache is not None:
        query &= len(cache) != 0
    if query:
        if cache is not None:
            try:
                new = DataQualityDict.read(cache, list(allflags))
            except IORegistryError as e:
                # can remove when astropy >= 1.2 is required
                if type(e) is not IORegistryError:
                    raise
                if len(allflags) == 1:
                    f = list(allflags)[0]
                    new = DataQualityDict()
                    new[f] = DataQualityFlag.read(cache, f, coalesce=False)
            for f in new:
                new[f].known &= newsegs
                new[f].active &= newsegs
                if coalesce:
                    new[f].coalesce()
                vprint("    Read %d segments for %s (%.2f%% coverage).\n"
                       % (len(new[f].active), f,
                          float(abs(new[f].known))/float(abs(newsegs))*100))
        else:
            if len(newsegs) >= 10:
                qsegs = span
            else:
                qsegs = newsegs
            # parse configuration for query
            kwargs = {}
            if url is not None:
                kwargs['url'] = url
            else:
                try:
                    kwargs['url'] = config.get('segment-database', 'url')
                except (NoSectionError, NoOptionError):
                    pass
            if kwargs.get('url', None) in SEGDB_URLS:
                query_func = DataQualityDict.query_segdb
            else:
                query_func = DataQualityDict.query_dqsegdb
            try:
                new = query_func(allflags, qsegs, on_error=segdb_error,
                                 **kwargs)
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
                new[f].known &= newsegs
                new[f].active &= newsegs
                if coalesce:
                    new[f].coalesce()
                vprint("    Downloaded %d segments for %s (%.2f%% coverage).\n"
                       % (len(new[f].active), f,
                          float(abs(new[f].known))/float(abs(newsegs))*100))
        # record new segments
        globalv.SEGMENTS += new
        for f in new:
            globalv.SEGMENTS[f].description = str(new[f].description)

    # return what was asked for
    if return_:
        for compound in flags:
            union, intersection, exclude, notequal = split_compound_flag(
                compound)
            if len(union + intersection) == 1:
                out[compound].description = globalv.SEGMENTS[f].description
                out[compound].padding = padding.get(f, (0, 0))
            for flist, op in zip([exclude, intersection, union, notequal],
                                 [operator.sub, operator.and_, operator.or_,
                                  not_equal]):
                for f in flist:
                    pad = padding.get(f, (0, 0))
                    segs = globalv.SEGMENTS[f].copy()
                    if isinstance(pad, (float, int)):
                        segs = segs.pad(pad, pad)
                    elif pad is not None:
                        segs = segs.pad(*pad)
                    if coalesce:
                        segs = segs.coalesce()
                    out[compound] = op(out[compound], segs)
            out[compound].known &= validity
            out[compound].active &= validity
            if coalesce:
                out[compound].coalesce()
        if isinstance(flag, string_types):
            return out[flag]
        else:
            return out


def not_equal(a, b, f):
    diff1 = a - b
    diff2 = b - a
    return diff1 | diff2


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


def format_padding(flags, padding):
    """Format an arbitrary collection of paddings into a `dict`
    """
    # parse string to start with
    if isinstance(padding, string_types):
        padding = list(eval(str))
    # zip list into dict
    if (isinstance(padding, (list)) or
            (isinstance(padding, tuple) and len(padding) and
             (any(isinstance(p, (list, tuple)) for p in padding) or
              len(padding) > 2))):
        return OrderedDict(list(zip(flags, padding)))
    # otherwise copy single padding param for all flags
    elif not isinstance(padding, dict):
        return OrderedDict((c, padding) for c in flags)
    else:
        return padding
