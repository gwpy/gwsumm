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

try:
    from configparser import (ConfigParser, NoSectionError, NoOptionError)
except ImportError:
    from ConfigParser import (ConfigParser, NoSectionError, NoOptionError)

from glue.ligolw.table import StripTableName as strip_table_name

from gwpy.table import lsctables
from gwpy.table.io import trigfind

from gwpy.table.utils import (get_table_column, get_row_value)
from gwpy.segments import (DataQualityFlag, SegmentList)

from . import globalv
from .utils import (re_cchar, vprint)
from .data import find_cache_segments

ETG_TABLE = {
    # single-IFO burst
    'omicron': lsctables.SnglBurstTable,
    'omega': lsctables.SnglBurstTable,
    'omegadq': lsctables.SnglBurstTable,
    'kleinewelle': lsctables.SnglBurstTable,
    'kw': lsctables.SnglBurstTable,
    'dmtomega': lsctables.SnglBurstTable,
    'dmt_wsearch': lsctables.SnglBurstTable,
    # multi-IFO burst
    'cwb': lsctables.MultiBurstTable,
    # single-IFO inspiral
    'daily_ihope': lsctables.SnglInspiralTable,
    'daily_ahope': lsctables.SnglInspiralTable,
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
        return ETG_TABLE[etg.lower()]
    except KeyError as e:
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
            table = lsctables.TableByName[etg]
        except KeyError as e:
            e.args = ('Cannot parse table name %r' % table,)
    if etg.lower() in ETG_TABLE and not force:
        raise KeyError('LIGO_LW table already registered to etg %r' % etg,)
    ETG_TABLE[etg.lower()] = table
    return table


def get_triggers(channel, etg, segments, config=ConfigParser(), cache=None,
                 query=True, multiprocess=False, tablename=None,
                 return_=True):
    """Read a table of transient event triggers for a given channel.
    """
    key = '%s,%s' % (str(channel), etg.lower())
    if isinstance(segments, DataQualityFlag):
        segments = segments.active
    segments = SegmentList(segments)

    # get LIGO_LW table for this etg
    if tablename:
        TableClass = lsctables.TableByName[tablename]
    else:
        TableClass = get_etg_table(etg)

    # work out columns
    try:
        columns = config.get(etg, 'columns').split(',')
    except (NoSectionError, NoOptionError):
        columns = TableClass.validcolumns.keys()
    else:
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
        for segment in new:
            filter_ = lambda t: float(get_row_value(t, 'time')) in segment
            # find trigger files
            if cache is None and etg.lower() == 'hacr':
                raise NotImplementedError("HACR parsing has not been "
                                          "implemented.")
            else:
                if cache is None:
                    segcache = trigfind.find_trigger_urls(str(channel), etg,
                                                          segment[0],
                                                          segment[1])
                    form = 'ligolw'
                else:
                    segcache = cache.sieve(segment=segment)
                    form = etg.lower()
                # read triggers and store
                segcache = segcache.checkfilesexist()[0]
                table = TableClass.read(segcache, columns=columns,
                                        format=form, filt=filter_)
            globalv.TRIGGERS[key].extend(table)
            csegs = find_cache_segments(segcache)
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
        return out
    else:
        return
