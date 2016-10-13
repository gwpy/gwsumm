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

import numpy
from numpy.lib import recfunctions

from glue.lal import Cache
from glue.ligolw.table import (StripTableName as strip_table_name,
                               CompareTableNames as compare_table_names)
from glue.ligolw.ligolw import PartialLIGOLWContentHandler

from gwpy.io.cache import cache_segments
from gwpy.table import (lsctables, GWRecArray)
from gwpy.table.utils import get_table_column
from gwpy.table.io.pycbc import filter_empty_files as filter_pycbc_live_files
from gwpy.table.io.hacr import get_hacr_triggers
from gwpy.segments import (DataQualityFlag, SegmentList)

try:
    import trigfind
except ImportError:
    from gwpy.table.io import trigfind

from . import globalv
from .utils import (re_cchar, vprint, count_free_cores, safe_eval)
from .config import (GWSummConfigParser, NoSectionError, NoOptionError)
from .channels import get_channel

TRIGFIND_FORMAT = {
    re.compile('omicron', re.I): 'sngl_burst',
    trigfind.pycbc_live: 'pycbc_live',
    trigfind.kleinewelle: 'sngl_burst',
    trigfind.dmt_omega: 'sngl_burst',
}

ETG_TABLE = lsctables.TableByName.copy()
ETG_TABLE.update({
    # single-IFO burst
    'omicron': lsctables.SnglBurstTable,
    'omega': lsctables.SnglBurstTable,
    'omegadq': lsctables.SnglBurstTable,
    'kleinewelle': lsctables.SnglBurstTable,
    'kw': lsctables.SnglBurstTable,
    'dmt_omega': lsctables.SnglBurstTable,
    'dmt_wsearch': lsctables.SnglBurstTable,
    # multi-IFO burst
    'cwb': lsctables.SnglBurstTable,
    # single-IFO inspiral
    'daily_ihope': lsctables.SnglInspiralTable,
    'daily_ahope': lsctables.SnglInspiralTable,
})
CONTENT_HANDLERS = {}


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
        return ETG_TABLE[etg.lower()]
    except KeyError as e:
        try:
            return ETG_TABLE[re_cchar.sub('_', etg.lower())]
        except KeyError:
            e.args = ('No LIGO_LW table registered to etg %r' % etg,)
            raise


def register_etg_table(etg, table, force=False):
    """Register a specific LIGO_LW table to an ETG

    Parameters
    ----------
    etg : `str`
        name of Event Trigger Generator to register
    table : `~gwpy.table.Table`, `str`
        `Table` class to register, or the ``tableName`` of the relevant class
    force : `bool`, optional, default: `False`
        overwrite an existing registration for the given ETG

    Raises
    ------
    KeyError
        if a `str` table cannot be resolved to a specific class
    """
    if isinstance(table, str):
        try:
            table = lsctables.TableByName[table]
        except KeyError as e:
            e.args = ('Cannot parse table name %r' % table,)
    if etg.lower() in ETG_TABLE and not force:
        raise KeyError('LIGO_LW table already registered to etg %r' % etg,)
    ETG_TABLE[etg.lower()] = table
    return table


def get_partial_contenthandler(table):
    """Build a `PartialLIGOLWContentHandler` for the given table

    Parameters
    ----------
    table : `type`
        the table class to be read

    Returns
    -------
    contenthandler : `type`
        a subclass of `~glue.ligolw.ligolw.PartialLIGOLWContentHandler` to
        read only the given `table`
    """
    def _filter_func(name, attrs):
        if name == table.tagName and attrs.has_key('Name'):
            return compare_table_names(attrs.get('Name'), table.tableName) == 0
        else:
            return False

    class _ContentHandler(PartialLIGOLWContentHandler):
        def __init__(self, document):
            super(_ContentHandler, self).__init__(document, _filter_func)

    return _ContentHandler


def get_triggers(channel, etg, segments, config=GWSummConfigParser(),
                 cache=None, columns=None, query=True, multiprocess=False,
                 ligolwtable=None, return_=True):
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

    # read segments from global memory
    try:
        havesegs = globalv.TRIGGERS[key].segments
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
        if 'format' not in kwargs and 'ahope' not in etg.lower():
            kwargs['format'] = etg.lower()

        # if single file
        if cache is not None and len(cache) == 1:
            trigs = read_cache(cache, new, etg, nproc=nproc,
                               tableclass=TableClass, **kwargs)
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
                    else:
                        for regex in TRIGFIND_FORMAT:
                            if regex.match(etg):
                                kwargs['format'] = TRIGFIND_FORMAT[regex]
                                if TRIGFIND_FORMAT[regex] in lsctables.TableByName:
                                    kwargs['get_as_columns'] = True
                                break
                    if etg.lower() == 'omega':
                        kwargs.setdefault('format', 'omega')
                    else:
                        kwargs.setdefault('format', 'ligolw')
                elif cache is not None:
                    segcache = cache
                # read table
                if etg.lower() == 'hacr':
                    trigs = get_hacr_triggers(channel, segment[0], segment[1],
                                              columns=columns)
                    trigs.segments = SegmentList([segment])
                else:
                    trigs = read_cache(segcache, SegmentList([segment]), etg,
                                       nproc=nproc, tableclass=TableClass,
                                       **kwargs)
                if trigs is not None:
                    add_triggers(trigs, key)
                    ntrigs += len(trigs)
                vprint(".")
        vprint(" | %d events read\n" % ntrigs)

    # if asked to read triggers, but didn't actually read any,
    # create an empty table so that subsequent calls don't raise KeyErrors
    if query and key not in globalv.TRIGGERS:
        if TableClass is not None:
            tab = lsctables.New(TableClass, columns=columns).to_recarray(
                get_as_columns=True)
        else:
            tab = GWRecArray((0,), dtype=[(c, float) for c in columns])
        tab.segments = SegmentList()
        add_triggers(tab, key)

    # work out time function
    if return_:
        return keep_in_segments(globalv.TRIGGERS[key], segments, etg)
    else:
        return


def add_triggers(table, key, segments=None):
    """Add a `GWRecArray` to the global memory cache
    """
    if segments is None:
        segments = table.segments
    try:
        old = globalv.TRIGGERS[key]
    except KeyError:
        globalv.TRIGGERS[key] = table
        globalv.TRIGGERS[key].segments = segments
    else:
        segs = old.segments
        globalv.TRIGGERS[key] = recfunctions.stack_arrays(
            (old, table), asrecarray=True, usemask=False,
            autoconvert=True).view(type(table))
        globalv.TRIGGERS[key].segments = segs + segments
    globalv.TRIGGERS[key].segments.coalesce()


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
    out.segments = segmentlist & table.segments
    return out


def get_times(table, etg):
    if etg == 'pycbc_live':
        return table['end_time']
    # guess from mapped LIGO_LW table
    try:
        TableClass = get_etg_table(etg)
    except KeyError:
        pass
    else:
        tablename = strip_table_name(TableClass.tableName)
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
    return kwargs


def read_cache(cache, segments, etg, nproc=1, tableclass=None, **kwargs):
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
    tableclass : `type`, optional
        the :class:`glue.ligolw.table.Table` sub-class for this ETG
    **kwargs
        other keyword arguments are passed to the `GWRecArray.read` or
        `{tableclass}.read` methods

    Returns
    -------
    table : `~gwpy.table.GWRecArray`, `None`
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
    if (kwargs.get('format', None) == 'ligolw' and tableclass is not None and
            'contenthandler' not in kwargs):
        if tableclass not in CONTENT_HANDLERS:  # only generate once
            CONTENT_HANDLERS[tableclass] = get_partial_contenthandler(
                tableclass)
        kwargs['contenthandler'] = CONTENT_HANDLERS[tableclass]
        lsctables.use_in(kwargs['contenthandler'])
    try:  # try directly reading a numpy.recarray
        table = GWRecArray.read(cache, nproc=nproc, **kwargs)
    except Exception as e:  # back up to LIGO_LW
        if tableclass is not None and 'No reader' in str(e):
            try:
                table = tableclass.read(cache, **kwargs)
            except Exception:
                raise e
            else:
                table = table.to_recarray(get_as_columns=True)
        else:
            raise
    # append new events to existing table
    try:
        csegs = cache_segments(cache)
    except AttributeError:
        csegs = SegmentList()
    table.segments = csegs
    return keep_in_segments(table, segments, etg)
