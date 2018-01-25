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

import warnings

from astropy.table import vstack as vstack_tables

from glue.lal import (Cache, CacheEntry)
from glue.ligolw import lsctables

from gwpy.io.cache import cache_segments
from gwpy.table import (EventTable, filters as table_filters)
from gwpy.table.filter import parse_column_filters
from gwpy.table.io.pycbc import filter_empty_files as filter_pycbc_live_files
from gwpy.segments import (DataQualityFlag, SegmentList)

import trigfind

from . import globalv
from .utils import (re_cchar, vprint, count_free_cores, safe_eval)
from .config import (GWSummConfigParser, NoSectionError)
from .channels import get_channel

# build list of default keyword arguments for reading ETGs
ETG_READ_KW = {
    'cwb': {
        'format': 'root',
        'treename': 'waveburst',
    },
    'daily_ihope': {
        'format': 'ligolw',
        'tablename': 'sngl_inspiral',
        'use_numpy_dtypes': True,
    },
    'daily_ahope': {
        'format': 'ligolw',
        'tablename': 'sngl_inspiral',
        'use_numpy_dtypes': True,
    },
    'dmt_omega': {
        'format': 'ligolw',
        'tablename': 'sngl_burst',
        'use_numpy_dtypes': True,
    },
    'dmt_wsearch': {
        'format': 'ligolw',
        'tablename': 'sngl_burst',
        'use_numpy_dtypes': True,
    },
    'kleinewelle': {
        'format': 'ligolw',
        'tablename': 'sngl_burst',
        'use_numpy_dtypes': True,
    },
    'kw': {
        'format': 'ligolw',
        'tablename': 'sngl_burst',
        'use_numpy_dtypes': True,
    },
    'omega': {
        'format': 'ascii',
    },
    'omegadq': {
        'format': 'ascii',
    },
    'omicron': {
        'format': 'ligolw',
        'tablename': 'sngl_burst',
        'use_numpy_dtypes': True,
    },
    'pycbc_live': {
        'format': 'hdf5.pycbc_live',
        'timecolumn': 'end_time',
        'extended_metadata': False,
    },
}

# set default for all LIGO_LW
for name in lsctables.TableByName:
    ETG_READ_KW[name] = {
        'format': 'ligolw',
        'tablename': name,
        'use_numpy_dtype': True,
    }


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
        kw_ = get_etg_read_kwargs(etg)
        form = kw_['format']
        tablename = kw_['tablename']
    except KeyError as e:
        e.args = ('No LIGO_LW table registered to etg %r' % etg,)
        raise
    if form == 'ligolw':
        return lsctables.TableByName[tablename]
    raise KeyError("No LIGO_LW table registered to etg %r" % etg)


def get_triggers(channel, etg, segments, config=GWSummConfigParser(),
                 cache=None, columns=None, format=None, query=True,
                 multiprocess=False, ligolwtable=None, filter=None,
                 timecolumn=None, return_=True):
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

    # get read keywords for this etg
    read_kw = get_etg_read_kwargs(etg, config=config, exclude=[])

    # extract columns (using function keyword if given)
    if columns:
          read_kw['columns'] = columns
    columns = read_kw.pop('columns', None)

    # parse filters
    if filter:
        read_kw.setdefault('selection', [])
        read_kw['selection'].extend(parse_column_filters(filter))

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

        # -- setup ----------

        # get find/read kwargs
        trigfindkwargs = dict(
             (k[9:], read_kw.pop(k)) for k in read_kw.keys() if
             k.startswith('trigfind-'))
        trigfindetg = trigfindkwargs.pop('etg', etg)

        # override with user options
        if format:
            read_kw['format'] = format
        elif not read_kw.get('format', None):
            read_kw['format'] = etg.lower()
        if timecolumn:
            read_kw['timecolumn'] = timecolumn
        elif columns is not None and 'time' in columns:
            read_kw['timecolumn'] = 'time'

        # replace columns keyword
        if read_kw['format'].startswith('ascii.'):
            read_kw['include_names'] = columns
        else:
            read_kw['columns'] = columns

        # customise kwargs for this ETG
        if etg.lower().replace('-', '_') in ['pycbc_live']:
            read_kw['ifo'] = get_channel(channel).ifo
        if etg.lower() in ['kw', 'kleinewelle']:
            read_kw['selection'].append('channel == "%s"' % channel)

        # filter on segments
        if 'timecolumn' in read_kw:
            read_kw['selection'].append((
                read_kw['timecolumn'], table_filters.in_segmentlist, new))

        # -- read -----------

        # if single file
        if cache is not None and len(cache) == 1:
            trigs = read_cache(cache, new, etg, nproc=nproc, **read_kw)
            if trigs is not None:
                add_triggers(trigs, key)
                ntrigs += len(trigs)
        # otherwise, loop over segments
        else:
            for segment in new:
                # find trigger files
                if cache is None and not etg.lower() == 'hacr':
                    try:
                        segcache = trigfind.find_trigger_files(
                            str(channel), trigfindetg, segment[0], segment[1],
                            **trigfindkwargs)
                    except ValueError as e:
                        warnings.warn("Caught %s: %s"
                                      % (type(e).__name__, str(e)))
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
                                       nproc=nproc, **read_kw)
                if trigs is not None:
                    add_triggers(trigs, key)
                    ntrigs += len(trigs)
                vprint(".")
        vprint(" | %d events read\n" % ntrigs)

    # if asked to read triggers, but didn't actually read any,
    # create an empty table so that subsequent calls don't raise KeyErrors
    if query and key not in globalv.TRIGGERS:
        # find LIGO_LW table for this ETG
        try:
            if columns is not None:  # don't need to map to LIGO_LW
                raise KeyError
            TableClass = get_etg_table(etg)
        except KeyError:  # build simple table
            tab = EventTable(names=columns)
        else:  # map to LIGO_LW table with full column listing
            tab = EventTable(lsctables.New(TableClass), get_as_columns=True)
        tab.meta['segments'] = SegmentList()
        add_triggers(tab, key)

    # work out time function
    if return_:
        try:
            return keep_in_segments(globalv.TRIGGERS[key], segments, etg)
        except KeyError:  # filtering didn't work, don't really care
            return globalv.TRIGGERS[key]
    else:
        return


def add_triggers(table, key, segments=None):
    """Add a `EventTable` to the global memory cache
    """
    if segments is not None:
        table.meta['segments'] = segments
    try:
        old = globalv.TRIGGERS[key]
    except KeyError:
        globalv.TRIGGERS[key] = table
    else:
        globalv.TRIGGERS[key] = vstack_tables((old, table))
        try:
            globalv.TRIGGERS[key].meta['segments'].coalesce()
        except KeyError:
            globalv.TRIGGERS[key].meta['segments'] = SegmentList()


def keep_in_segments(table, segmentlist, etg=None):
    """Return a view of the table containing only those rows in the segmentlist
    """
    times = get_times(table, etg)
    keep = table_filters.in_segmentlist(times, segmentlist)
    out = table[keep]
    out.meta['segments'] = segmentlist & table.meta['segments']
    return out


def get_times(table, etg):
    # allow user to have selected the time column
    if table.meta.get('timecolumn'):
        return table[table.meta['timecolumn']]
    # otherwise search for it
    try:
        return table['time']
    except KeyError:
        # shortcut pycbc
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
        raise


def get_time_column(table, etg):
    """Get the time column for this table, etg pair
    """
    if table.meta.get('timecolumn'):
        return table.meta['timecolumn']
    if 'time' in table.colnames:
        return 'time'
    if etg == 'pycbc_live':
        return 'end_time'

    # handle LIGO_LW tables (as best we can)
    tablename = ETG_READ_KW.get(etg, {}).get('tablename', None)
    if tablename in ('sngl_burst', 'multi_burst'):
        return 'peak'
    if tablename in ('sngl_inspiral', 'coinc_inspiral', 'multi_inspiral'):
        return 'end'
    if tablename in ('sngl_ringdown', 'multi_ringdown'):
        return 'start'
    raise ValueError("unknown time column for table")


def get_etg_read_kwargs(etg, config=None, exclude=['columns']):
    """Read keyword arguments to pass to the trigger reader for a given etg
    """
    # use global defaults
    kwargs = {
        'format': None,
        'selection': [],
    }

    # update with ETG defaults
    kwargs.update(ETG_READ_KW.get(re_cchar.sub('_', etg).lower(), {}))

    # get kwargs from config
    if config is not None:
        try:
            kwargs.update(config.nditems(etg))
        except NoSectionError:
            try:
                kwargs.update(config.nditems(etg.lower()))
            except NoSectionError:
                return {}

    # format kwargs
    for key in list(kwargs.keys()):
        # remove unwanted keys
        if key in exclude:
            kwargs.pop(key)
            continue
        # eval string to object (safely)
        kwargs[key] = safe_eval(kwargs[key])
        if key.endswith('columns') and isinstance(kwargs[key], str):
            kwargs[key] = kwargs[key].replace(' ', '').split(',')

    if 'selection' in kwargs:
        kwargs['selection'] = parse_column_filters(kwargs['selection'])

    return kwargs


def read_cache(cache, segments, etg, nproc=1, timecolumn=None, **kwargs):
    """Read a table of events from a cache

    This function is mainly meant for use from the `get_triggers` method

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
        cache = cache.pfnlist()  # some readers only like filenames
    if etg == 'pycbc_live':  # remove empty HDF5 files
        cache = filter_pycbc_live_files(cache, ifo=kwargs['ifo'])

    if len(cache) == 0:
        return

    # read triggers
    table = EventTable.read(cache, **kwargs)
    if timecolumn:
        table.meta['timecolumn'] = timecolumn

    # get back from cache entry
    if isinstance(cache, CacheEntry):
        cache = Cache([cache])

    # append new events to existing table
    try:
        csegs = cache_segments(cache) & segments
    except (AttributeError, TypeError):
        csegs = SegmentList()
    table.meta['segments'] = csegs

    if timecolumn:  # already filtered on-the-fly
        return table
    # filter now
    return keep_in_segments(table, segments, etg)
