
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
