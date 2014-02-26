# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013)
#
# This file is part of GWSumm
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
# along with GWSumm.  If not, see <http://www.gnu.org/licenses/>

"""Custom `SummaryTab` for the output of the Daily iHope pipeline.
"""

import os

from lal import gpstime

from glue.lal import Cache
from glue.ligolw import (utils as llwutils)
from glue.ligolw.lsctables import (SnglInspiralTable, SummValueTable)

from gwpy.timeseries import (TimeSeries, TimeSeriesList)
from gwpy.plotter.table import (get_table_column, get_row_value)

from .. import (version, html, globalv)
from ..config import GWSummConfigParser
from ..data import find_cache_segments
from ..triggers import get_triggers
from ..utils import re_quote
from .registry import (get_tab, register_tab)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

base = get_tab('default')


class DailyIhopeTab(base):
    """Custom tab displaying a summary of Daily iHope results.
    """
    type = 'daily-ihope'

    @classmethod
    def from_ini(cls, config, section, plotdir=os.curdir, base=''):
        """Define a new `DailyIhopeTab` from a `ConfigParser`.
        """
        # parse generic configuration
        new = super(DailyIhopeTab, cls).from_ini(config, section,
                                                 plotdir=plotdir, base=base)
        new.channel = re_quote.sub('', config.get(section, 'channel'))

        # work out day directory and url
        utc = gpstime.gps_to_utc(new.span[0])
        basedir = os.path.normpath(config.get(section, 'base-directory'))
        daydir = os.path.join(basedir, utc.strftime('%Y%m'),
                              utc.strftime('%Y%m%d'))
        home_, postbase = daydir.split('/public_html/', 1)
        user = os.path.split(home_)[1]
        new.ihopepage = '/~%s/%s/' % (user, postbase)

        # get cache options
        cachefile = config.get(section, 'inspiral-cache')
        inspiralcache = os.path.join(daydir,  cachefile)
        with open(inspiralcache, 'r') as fobj:
            new.inspiralcache = Cache.fromfile(fobj).sieve(segment=new.span)
        cachefile = config.get(section, 'tmpltbank-cache')
        tmpltbankcache = os.path.join(daydir, cachefile)
        with open(tmpltbankcache, 'r') as fobj:
            new.tmpltbankcache = Cache.fromfile(fobj).sieve(segment=new.span)

        # get loudest options
        if config.has_option(section, 'loudest'):
            new.loudest = {'N': config.getint(section, 'loudest')}
            if config.has_option(section, 'loudest-columns'):
                new.loudest['columns'] = map(
                    lambda s: re_quote.sub('', s),
                    config.get(section, 'loudest-columns').split(','))
            else:
                new.loudest['columns'] = ['end', 'new_snr', 'snr', 'mchirp',
                                          'mass1', 'mass2', 'reduced_chisq']
            if config.has_option(section, 'loudest-labels'):
                new.loudest['labels'] = map(
                    lambda s: re_quote.sub('', s),
                    config.get(section, 'loudest-labels').split(','))
            else:
                new.loudest['labels'] = [' '.join(map(str.title, s.split('_')))
                                         for s in new.loudest['columns']]
            if config.has_option(section, 'loudest-rank'):
                new.loudest['rank'] = re_quote.sub(
                    '', config.get(section, 'loudest-rank'))
            else:
                new.loudest['rank'] = 'new_snr'
        else:
            new.loudest = None

        return new

    def get_tmpltbank_data(self):
        """Read inspiral horizon data from the TmpltBank cache
        """
        tmpltsegs = find_cache_segments(self.tmpltbankcache)
        ifo = self.channel.split(':')[0]
        rangechannel = '%s:horizon_distance' % ifo
        sizechannel = '%s:tmpltbank_size' % ifo
        globalv.DATA[rangechannel] = TimeSeriesList()
        globalv.DATA[sizechannel] = TimeSeriesList()
        for seg in tmpltsegs:
            segcache = self.tmpltbankcache.sieve(segment=seg)
            rangedata = []
            sizedata = []
            for ce in segcache:
                xmldoc = llwutils.load_filename(ce.path)
                svtable = SummValueTable.get_table(xmldoc)
                svtable.sort(key=lambda row: float(row.comment.split('_')[0]))
                rangedata.append(svtable[0].value * (1.4)**(5/6.))
                sizedata.append(len(SnglInspiralTable.get_table(xmldoc)))
            if rangedata:
                dt = float(abs(segcache[0].segment))
                epoch = segcache[0].segment[0] + dt/2.
                globalv.DATA[rangechannel].append(
                    TimeSeries(rangedata, sample_rate=1/dt, epoch=epoch,
                               name=rangechannel))
                globalv.DATA[rangechannel].coalesce()
                globalv.DATA[sizechannel].append(
                    TimeSeries(sizedata, sample_rate=1/dt, epoch=epoch,
                               name=sizechannel))
                globalv.DATA[sizechannel].coalesce()

    def process_state(self, state, nds='guess', multiprocess=False,
                      config=GWSummConfigParser()):
        self.get_tmpltbank_data()
        super(DailyIhopeTab, self).process_state(state, nds=nds,
                                                 multiprocess=multiprocess,
                                                 config=config,
                                                 datacache=Cache(),
                                                 trigcache=self.inspiralcache)

    def build_inner_html(self, state):
        """Write the '#main' HTML content for this `DailyIhopeTab`.
        """
        page = self.scaffold_plots(state)

        # link full results
        page.div(class_='btn-group')
        page.a('Click here for the full Daily Ihope results',
               href=self.ihopepage, rel='external', target='_blank',
               class_='btn btn-default btn-info btn-xl')
        page.div.close()

        if self.loudest:
            table = get_triggers(self.channel, 'Daily Ihope', state,
                                 query=False)
            rank = get_table_column(table, self.loudest['rank']).argsort()
            indexes = rank[-self.loudest['N']:][::-1]
            page.h1('Loudest events')
            page.p('The following table displays the %d loudest events as '
                   'recorded by Daily Ihope.' % self.loudest['N'])
            headers = self.loudest['labels']
            if 'time' in headers[0]:
                headers.insert(1, 'UTC time')
                date = True
            else:
                data = False
            data = []
            for idx in indexes:
                row = table[idx]
                data.append([])
                for column in self.loudest['columns']:
                    data[-1].append('%.3f' % float(get_row_value(row, column)))
                if date:
                    data[-1].insert(1, gpstime.tconvert(
                        row.get_end(), '%B %d %Y, %H:%M:%S.%f')[:-3])
            page.add(str(html.data_table(headers, data, table='data')))

        return page

register_tab(DailyIhopeTab)
