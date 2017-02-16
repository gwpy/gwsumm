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

"""Read and store transient event triggers
"""

import re
import warnings

from six import string_types

import numpy

from astropy.table import vstack as vstack_tables

from glue.lal import Cache

from gwpy.io.cache import cache_segments
from gwpy.table import (lsctables, EventTable)
from gwpy.table.utils import get_table_column
from gwpy.table.io.pycbc import filter_empty_files as filter_pycbc_live_files
from gwpy.segments import (DataQualityFlag, SegmentList)

try:
    import trigfind
except ImportError:
    from gwpy.table.io import trigfind

from . import globalv
from .utils import (re_cchar, vprint, count_free_cores, safe_eval)
from .config import (GWSummConfigParser, NoSectionError, NoOptionError)
from .channels import get_channel

# build list of format arguments for reading ETGs
ETG_FORMAT = {
    'cwb': 'root.cwb',
    'daily_ihope': 'ligolw.sngl_inspiral',
    'daily_ahope': 'ligolw.sngl_inspiral',
    'dmt_omega': 'ligolw.sngl_burst',
    'dmt_wsearch': 'ligolw.sngl_burst',
    'kleinewelle': 'ligolw.sngl_burst',
    'kw': 'ligolw.sngl_burst',
    'omega': 'ascii',
    'omegadq': 'ascii',
    'omicron': 'ligolw.sngl_burst',
    'pycbc_live': 'hdf5.pycbc_live',
}
for name in lsctables.TableByName:
    ETG_FORMAT[name] = 'ligolw.%s' % name


def get_etg_format(etg):
    return ETG_FORMAT[re_cchar.sub('_', etg).lower()]


def get_etg_table(etg):
    """Find which table should be used for the given etg

    Parameters
    ----------
    etg : `str`
        name of Event Trigger Generator for which to query

    Returns
    -------
    table : `~gwpy.table.Table`
        LIGO_LW table registered to the given ETG

    Raises
    ------
    KeyError
        if the ETG is not registered
    """
    try:
        format_ = get_etg_format(etg)
    except KeyError as e:
        e.args = ('No LIGO_LW table registered to etg %r' % etg,)
        raise
    if format_.startswith('ligolw.'):
        return lsctables.TableByName[format_[7:]]
    raise KeyError("No LIGO_LW table registered to etg %r" % etg)


def get_triggers(channel, etg, segments, config=GWSummConfigParser(),
                 cache=None, columns=None, format=None, query=True,
                 multiprocess=False, ligolwtable=None, filterstr=None,
                 return_=True):
    """Read a table of transient event triggers for a given channel.
    """
    key = '%s,%s' % (str(channel), etg.lower())
    # convert input segments to a segmentlist (for convenience)
    if isinstance(segments, DataQualityFlag):
        segments = segments.active
    segments = SegmentList(segments)

    # get processes
    if multiprocess is True:
        nproc = count_free_cores()
    elif multiprocess is False:
        nproc = 1
    else:
        nproc = multiprocess

    # find LIGO_LW table for this ETG
    try:
        TableClass = get_etg_table(etg)
    except KeyError:
        TableClass = None

    # work out columns
    if columns is None:
        try:
            columns = config.get(etg.lower(), 'columns').split(',')
        except (NoSectionError, NoOptionError):
            columns = None
        else:
            columns = [c.strip() for c in columns]

    # read segments from global memory
    try:
        havesegs = globalv.TRIGGERS[key].meta['segments']
    except KeyError:
        new = segments
    else:
        new = segments - havesegs

    # read new triggers
    if query and abs(new) != 0:
        ntrigs = 0
        vprint("    Grabbing %s triggers for %s" % (etg, str(channel)))

        # store read kwargs
        kwargs = get_etg_read_kwargs(config, etg, exclude=['columns'])
        kwargs['columns'] = columns
        if etg.lower().replace('-', '_') in ['cwb', 'pycbc_live']:
            kwargs['ifo'] = get_channel(channel).ifo
        if format is not None:
            kwargs['format'] = format
        if 'format' not in kwargs:
            try:
                kwargs['format'] = get_etg_format(etg)
            except KeyError:
                kwargs['format'] = etg.lower()

        # if single file
        if cache is not None and len(cache) == 1:
            trigs = read_cache(cache, new, etg, nproc=nproc,
                               filterstr=filterstr, **kwargs)
            if trigs is not None:
                add_triggers(trigs, key)
                ntrigs += len(trigs)
        # otherwise, loop over segments
        else:
            for segment in new:
                # find trigger files
                if cache is None and etg.lower() in ['kw', 'kleinewelle']:
                    kwargs['filt'] = lambda t: t.channel == str(channel)
                if cache is None and not etg.lower() == 'hacr':
                    try:
                        segcache = trigfind.find_trigger_urls(str(channel), etg,
                                                              segment[0],
                                                              segment[1])
                    except ValueError as e:
                        warnings.warn("Caught %s: %s" % (type(e).__name__, str(e)))
                        continue
                elif cache is not None:
                    segcache = cache
                # read table
                if etg.lower() == 'hacr':
                    from gwpy.table.io.hacr import get_hacr_triggers
                    trigs = get_hacr_triggers(channel, segment[0], segment[1],
                                              columns=columns)
                    trigs.meta['segments'] = SegmentList([segment])
                else:
                    trigs = read_cache(segcache, SegmentList([segment]), etg,
                                       nproc=nproc,
                                       filterstr=filterstr,
                                       **kwargs)
                if trigs is not None:
                    add_triggers(trigs, key)
                    ntrigs += len(trigs)
                vprint(".")
        vprint(" | %d events read\n" % ntrigs)

    # if asked to read triggers, but didn't actually read any,
    # create an empty table so that subsequent calls don't raise KeyErrors
    if query and key not in globalv.TRIGGERS:
        if columns is None and TableClass is not None:
            tab = EventTable(lsctables.New(TableClass), get_as_columns=True)
        else:
            tab = EventTable(names=columns)
        tab.meta['segments'] = SegmentList()
        add_triggers(tab, key)

    # work out time function
    if return_:
        return keep_in_segments(globalv.TRIGGERS[key], segments, etg)
    else:
        return


def add_triggers(table, key, segments=None):
    """Add a `EventTable` to the global memory cache
    """
    if segments is None:
        segments = table.meta['segments']
    try:
        old = globalv.TRIGGERS[key]
    except KeyError:
        globalv.TRIGGERS[key] = table
        globalv.TRIGGERS[key].meta['segments'] = segments
    else:
        segs = old.meta['segments']
        globalv.TRIGGERS[key] = vstack_tables((old, table))
        globalv.TRIGGERS[key].meta['segments'] = segs + segments
    globalv.TRIGGERS[key].meta['segments'].coalesce()


def time_in_segments(times, segmentlist):
    """Find which times lie inside a segmentlist

    Parameters
    ----------
    times : `numpy.ndarray`
        1-dimensional array of times
    segmentlist : :class:`~glue.segments.segmentlist`
        list of time segments to test

    Returns
    -------
    insegments : `numpy.ndarray`
        boolean array indicating which times lie inside the segmentlist

    Notes
    -----
    A time `t` lies inside a segment `[a..b]` if `a <= t <= b`
    """
    segmentlist = type(segmentlist)(segmentlist).coalesce()
    keep = numpy.zeros(times.shape[0], dtype=bool)
    j = 0
    try:
        a, b = segmentlist[j]
    except IndexError:  # no segments, return all False
        return keep
    i = 0
    while i < keep.size:
        t = times[i]
        # if before start, move to next trigger now
        if t < a:
            i += 1
            continue
        # if after end, find the next segment and check this trigger again
        if t > b:
            j += 1
            try:
                a, b = segmentlist[j]
                continue
            except IndexError:
                break
        # otherwise it must be in this segment, record and move to next
        keep[i] = True
        i += 1
    return keep


def keep_in_segments(table, segmentlist, etg=None):
    """Return a view of the table containing only those rows in the segmentlist
    """
    times = get_times(table, etg)
    order = times.argsort()
    t2 = table[order]
    times = times[order]
    keep = time_in_segments(times, segmentlist)
    out = t2[keep]
    out.meta['segments'] = segmentlist & table.meta['segments']
    return out


def get_times(table, etg):
    try:
        return table['time']
    except KeyError:
        if etg == 'pycbc_live':
            return table['end_time']
        # guess from mapped LIGO_LW table
        try:
            TableClass = get_etg_table(etg)
        except KeyError:
            pass
        else:
            tablename = TableClass.TableName(TableClass.tableName)
            if tablename.endswith('_burst'):
                return table['peak_time'] + table['peak_time_ns'] * 1e-9
            if tablename.endswith('_inspiral'):
                return table['end_time'] + table['end_time_ns'] * 1e-9
            if tablename.endswith('_ringdown'):
                return table['start_time'] + table['start_time_ns'] * 1e-9
        # use gwpy method (not guaranteed to work)
        return get_table_column(table, 'time').astype(float)


def get_etg_read_kwargs(config, etg, exclude=['columns']):
    """Read keyword arguments to pass to the trigger reader for a given etg
    """
    try:
        kwargs = dict(config.nditems(etg))
    except NoSectionError:
        return {}
    for key in kwargs.keys():
        if key in exclude:
            kwargs.pop(key)
            continue
        kwargs[key] = safe_eval(kwargs[key])
        if key.endswith('columns') and isinstance(kwargs[key], str):
            kwargs[key] = kwargs[key].replace(' ', '').split(',')
    return kwargs

def filter_triggers(table, filterstr):
    for token in filterstr.split(' '):
        try:
            for col in table.dtype.names:
                token = token.replace(col, "table.%s" % col)
            table = table[eval(token)]
        except:
            raise ValueError('Failed to parse trigger filter str: %s' % filterstr)
    return table

def read_cache(cache, segments, etg, nproc=1, filterstr=None, **kwargs):
    """Read a table of events from a cache

    This function is mainly meant for use from the `get_triggers method

    Parameters
    ----------
    cache : :class:`glue.lal.Cache`
        the formatted list of files to read
    segments : `~gwpy.segments.SegmentList`
        the list of segments to read
    etg : `str`
        the name of the trigger generator that created the files
    nproc : `int`, optional
        the number of parallel processes to use when reading
    **kwargs
        other keyword arguments are passed to the `EventTable.read` or
        `{tableclass}.read` methods

    Returns
    -------
    table : `~gwpy.table.EventTable`, `None`
        a table of events, or `None` if the cache has no overlap with
        the segments
    """
    if isinstance(cache, Cache):
        cache = cache.sieve(segmentlist=segments)
        cache = cache.checkfilesexist()[0]
        cache.sort(key=lambda x: x.segment[0])
        if etg == 'pycbc_live':  # remove empty HDF5 files
            cache = type(cache)(
                filter_pycbc_live_files(cache, ifo=kwargs['ifo']))
    # if no files, skip
    if len(cache) == 0:
        return
    # read triggers
    table = EventTable.read(cache, nproc=nproc, **kwargs)

    # Filter table
    if filterstr is not None:
        table = filter_triggers(table, filterstr)

    # append new events to existing table
    try:
        csegs = cache_segments(cache)
    except AttributeError:
        csegs = SegmentList()
    table.meta['segments'] = csegs
    return keep_in_segments(table, segments, etg)
