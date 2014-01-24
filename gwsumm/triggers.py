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

from laldetchar import triggers as trigio
from laldetchar.triggers import trigfind

from glue.ligolw.lsctables import TableByName

from gwpy.plotter.table import get_table_column
from gwpy.segments import DataQualityFlag

from . import globalv
from .utils import (re_cchar, vprint)
from .data import find_cache_segments


def get_triggers(channel, etg, segments, config=ConfigParser(), cache=None,
                 query=True, multiprocess=False):
    """Read a table of transient event triggers for a given channel.
    """
    key = '%s,%s' % (str(channel), etg)
    if isinstance(segments, DataQualityFlag):
        segments = segments.active
    # get default table type for this generator
    TableClass = TableByName[trigio.utils.which_table(etg)]

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
            key, trigio.utils.new_ligolw_table(etg, columns=columns))
        globalv.TRIGGERS[key].segments = type(segments)()
    else:
        new = segments - havesegs

    # read new triggers
    query &= (abs(new) != 0)
    if query:
        for segment in new:
            # find trigger files
            if cache is None and etg.lower() == 'hacr':
                raise NotImplementedError("HACR parsing has not been "
                                          "implemented.")
            else:
                if cache is None:
                    segcache = trigfind.find_trigger_urls(str(channel), etg,
                                                          segment[0],
                                                          segment[1])
                else:
                    segcache = cache.sieve(segment=segment)
                # read triggers and store
                table = trigio.from_files(segcache, etg, columns=columns,
                                          channel=str(channel),
                                          start=segment[0],
                                          end=segment[1],
                                          verbose=globalv.VERBOSE)
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
    times = get_table_column(globalv.TRIGGERS[key], 'time').astype(float)

    # return correct triggers
    out = trigio.utils.new_ligolw_table(etg, columns=columns)
    out.extend(t for (i, t) in enumerate(globalv.TRIGGERS[key]) if
               times[i] in segments)
    return out
