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

import os.path
import glob
import re
import warnings
try:
    from configparser import (ConfigParser, NoSectionError, NoOptionError)
except ImportError:
    from ConfigParser import (ConfigParser, NoSectionError, NoOptionError)

from glue.lal import (Cache, CacheEntry)
from glue.ligolw.table import (StripTableName as strip_table_name,
                               CompareTableNames as compare_table_names)
from glue.ligolw.ligolw import PartialLIGOLWContentHandler

from gwpy.io.cache import cache_segments
from gwpy.table import lsctables

from gwpy.time import to_gps

from gwpy.table.utils import (get_table_column, get_row_value)
from gwpy.segments import (DataQualityFlag, SegmentList, Segment)

try:
    import trigfind
except ImportError:
    from gwpy.table.io import trigfind

from . import globalv
from .utils import (re_cchar, vprint)
from .channels import get_channel

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


def get_triggers(channel, etg, segments, config=ConfigParser(), cache=None,
                 query=True, multiprocess=False, tablename=None,
                 columns=None, contenthandler=None, return_=True):
    """Read a table of transient event triggers for a given channel.
    """
    key = '%s,%s' % (str(channel), etg.lower())
    if isinstance(segments, DataQualityFlag):
        segments = segments.active
    segments = SegmentList(segments)

    # get LIGO_LW table for this etg
    if tablename:
        TableClass = lsctables.TableByName[tablename]
        register_etg_table(etg, TableClass, force=True)
    elif key in globalv.TRIGGERS:
        TableClass = type(globalv.TRIGGERS[key])
    else:
        TableClass = get_etg_table(etg)

    # work out columns
    if columns is None:
        try:
            columns = config.get(etg, 'columns').split(',')
        except (NoSectionError, NoOptionError):
            if etg.lower() in ['cwb', 'cwb-ascii']:
                columns = None
            else:
                columns = TableClass.validcolumns.keys()
    if columns is not None:
        for col in ['process_id', 'search', 'channel']:
            if col not in columns:
                columns.append(col)

    # read segments from global memory
    try:
        havesegs = globalv.TRIGGERS[key].segments
    except KeyError:
        new = segments
        globalv.TRIGGERS.setdefault(
            key, lsctables.New(TableClass, columns=columns))
        globalv.TRIGGERS[key].segments = type(segments)()
    else:
        new = segments - havesegs

    # read new triggers
    query &= (abs(new) != 0)
    if query:
        # store read kwargs
        kwargs = {'columns': columns}

        # set content handler
        if contenthandler is None:
            contenthandler = get_partial_contenthandler(TableClass)
        lsctables.use_in(contenthandler)

        for segment in new:
            kwargs['filt'] = (
                lambda t: float(get_row_value(t, 'time')) in segment)
            # find trigger files
            if cache is None and etg.lower() == 'hacr':
                raise NotImplementedError("HACR parsing has not been "
                                          "implemented.")
            if cache is None and etg.lower() in ['kw', 'kleinewelle']:
                kwargs['filt'] = lambda t: (
                    float(get_row_value(t, 'time')) in segment and
                    t.channel == str(channel))
            if cache is None:
                try:
                    segcache = trigfind.find_trigger_urls(str(channel), etg,
                                                          segment[0],
                                                          segment[1])
                except ValueError as e:
                    warnings.warn("Caught %s: %s" % (type(e).__name__, str(e)))
                    continue
                if etg.lower() == 'omega':
                    kwargs['format'] = 'omega'
                else:
                    kwargs['format'] = 'ligolw'
            elif isinstance(cache, Cache):
                segcache = cache.sieve(segment=segment)
            else:
                segcache = cache
            if isinstance(segcache, Cache):
                segcache = segcache.checkfilesexist()[0]
            if 'format' not in kwargs and 'ahope' not in etg.lower():
                kwargs['format'] = etg.lower()
            if (issubclass(TableClass, lsctables.SnglBurstTable) and
                    etg.lower().startswith('cwb')):
                kwargs['ifo'] = get_channel(channel).ifo
            # read triggers and store
            if len(segcache) == 0:
                continue
            if kwargs.get('format', None) == 'ligolw':
                kwargs['contenthandler'] = contenthandler
            table = TableClass.read(segcache, **kwargs)
            globalv.TRIGGERS[key].extend(table)
            try:
                csegs = cache_segments(segcache)
            except AttributeError:
                csegs = SegmentList()
            try:
                globalv.TRIGGERS[key].segments.extend(csegs)
            except AttributeError:
                globalv.TRIGGERS[key].segments = csegs
            finally:
                globalv.TRIGGERS[key].segments.coalesce()
            vprint('\r')

    # work out time function
    if return_:
        times = get_table_column(globalv.TRIGGERS[key], 'time').astype(float)

        # return correct triggers
        out = lsctables.New(TableClass, columns=columns)
        out.channel = str(channel)
        out.etg = str(etg)
        out.extend(t for (i, t) in enumerate(globalv.TRIGGERS[key]) if
                   times[i] in segments)
        out.segments = segments & globalv.TRIGGERS[key].segments
        return out
    else:
        return
