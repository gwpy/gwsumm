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
try:
    from configparser import (ConfigParser, NoSectionError, NoOptionError)
except ImportError:
    from ConfigParser import (ConfigParser, NoSectionError, NoOptionError)

from glue.lal import (Cache, CacheEntry)
from glue.ligolw.table import (StripTableName as strip_table_name,
                               CompareTableNames as compare_table_names)
from glue.ligolw.ligolw import PartialLIGOLWContentHandler

from gwpy.table import lsctables
from gwpy.table.io import trigfind
from gwpy.time import to_gps

from gwpy.table.utils import (get_table_column, get_row_value)
from gwpy.segments import (DataQualityFlag, SegmentList, Segment)

from . import globalv
from .utils import (re_cchar, vprint)
from .channels import get_channel
from .data import find_cache_segments

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
                 contenthandler=None, return_=True):
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
    try:
        columns = config.get(etg, 'columns').split(',')
    except (NoSectionError, NoOptionError):
        if etg.lower() in ['cwb', 'cwb-ascii']:
            columns = None
        else:
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
            if cache is None and re.match('dmt(.*)omega', etg.lower()):
                segcache = find_dmt_omega(channel, segment[0], segment[1])
                kwargs['format'] = 'ligolw'
            elif cache is None and etg.lower() in ['kw', 'kleinewelle']:
                segcache = find_kw(channel, segment[0], segment[1])
                kwargs['format'] = 'ligolw'
                kwargs['filt'] = lambda t: (
                    float(get_row_value(t, 'time')) in segment and
                    t.channel == str(channel))
            elif cache is None:
                segcache = trigfind.find_trigger_urls(str(channel), etg,
                                                      segment[0],
                                                      segment[1])
                kwargs['format'] = 'ligolw'
            else:
                segcache = cache.sieve(segment=segment)
            if 'format' not in kwargs and 'ahope' not in etg.lower():
                kwargs['format'] = etg.lower()
            if (issubclass(TableClass, lsctables.SnglBurstTable) and
                    etg.lower().startswith('cwb')):
                kwargs['ifo'] = get_channel(channel).ifo
            # read triggers and store
            segcache = segcache.checkfilesexist()[0]
            if len(segcache) == 0:
                continue
            if kwargs.get('format', None) == 'ligolw':
                kwargs['contenthandler'] = contenthandler
            table = TableClass.read(segcache, **kwargs)
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
        out.segments = segments & globalv.TRIGGERS[key].segments
        return out
    else:
        return


def find_dmt_omega(channel, start, end, base=None):
    """Find DMT-Omega trigger XML files
    """
    span = Segment(to_gps(start), to_gps(end))
    channel = get_channel(channel)
    ifo = channel.ifo
    if base is None and channel.name.split(':', 1)[-1] == 'GDS-CALIB_STRAIN':
        base = '/gds-%s/dmt/triggers/%s-Omega_hoft' % (
            ifo.lower(), ifo[0].upper())
    elif base is None:
        raise NotImplementedError("This method doesn't know how to locate DMT "
                                  "Omega trigger files for %r" % str(channel))
    gps5 = int('%.5s' % start)
    end5 = int('%.5s' % end)
    out = Cache()
    append = out.append
    while gps5 <= end5:
        trigglob = os.path.join(
            base, str(gps5),
            '%s-%s_%s_%s_OmegaC-*-*.xml' % (
                ifo, channel.system, channel.subsystem, channel.signal))
        found = glob.glob(trigglob)
        for f in found:
            ce = CacheEntry.from_T050017(f)
            if ce.segment.intersects(span):
                append(ce)
        gps5 += 1
    out.sort(key=lambda e: e.path)
    vprint("    Found %d files for %s (DMT-Omega)\n"
           % (len(out), channel.ndsname))
    return out


def find_kw(channel, start, end, base=None):
    """Find KW trigger XML files
    """
    span = Segment(to_gps(start), to_gps(end))
    channel = get_channel(channel)
    ifo = channel.ifo
    if base is None and channel.name.split(':', 1)[-1] == 'GDS-CALIB_STRAIN':
        tag = '%s-KW_HOFT' % ifo[0].upper()
        base = '/gds-%s/dmt/triggers/%s' % (ifo.lower(), tag)
    elif base is None:
        tag = '%s-KW_TRIGGERS' % ifo[0].upper()
        base = '/gds-%s/dmt/triggers/%s' % (ifo.lower(), tag)
    gps5 = int('%.5s' % start)
    end5 = int('%.5s' % end)
    out = Cache()
    append = out.append
    while gps5 <= end5:
        trigglob = os.path.join(
            base, '%s-%d' % (tag, gps5), '%s-*-*.xml' % tag)
        found = glob.glob(trigglob)
        for f in found:
            ce = CacheEntry.from_T050017(f)
            if ce.segment.intersects(span):
                append(ce)
        gps5 += 1
    out.sort(key=lambda e: e.path)
    vprint("    Found %d files for %s (KW)\n"
           % (len(out), channel.ndsname))
    return out
