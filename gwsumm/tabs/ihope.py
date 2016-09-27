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
from warnings import warn

from astropy.io.registry import (register_reader, get_reader)

from glue.lal import Cache
from glue.ligolw import (utils as llwutils)
from glue.ligolw.lsctables import (SnglInspiralTable, SummValueTable)

from gwpy.segments import DataQualityFlag
from gwpy.time import from_gps
from gwpy.timeseries import (TimeSeries, TimeSeriesList)
from gwpy.plotter.table import (get_table_column, get_row_value)

from .. import (html, globalv)
from ..config import (GWSummConfigParser, NoOptionError, DEFAULTSECT)
from ..data import find_cache_segments
from ..triggers import (get_triggers, register_etg_table)
from ..utils import re_quote
from ..state import SummaryState
from ..mode import (Mode, get_mode)
from .registry import (get_tab, register_tab)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['DailyAhopeTab']

base = get_tab('default')


class DailyAhopeTab(base):
    """Custom tab displaying a summary of Daily iHope results.
    """
    type = 'daily-ahope'

    def __init__(self, *args, **kwargs):
        super(DailyAhopeTab, self).__init__(*args, **kwargs)
        register_etg_table(self.name.lower(), SnglInspiralTable)
        register_reader(self.name.lower(), SnglInspiralTable,
                        get_reader('ligolw', SnglInspiralTable))

    @classmethod
    def from_ini(cls, config, section, **kwargs):
        """Define a new `DailyAhopeTab` from a `ConfigParser`.
        """
        # parse states
        ifo = config.get(DEFAULTSECT, 'ifo')
        start = config.getint(DEFAULTSECT, 'gps-start-time')
        end = config.getint(DEFAULTSECT, 'gps-end-time')
        if config.has_option(section, 'states'):
            raise ValueError("DailyAhopeTab does not support configuration of "
                             "multiple states, please use the 'state' option "
                             "to name the Hveto state")
        try:
            state = re_quote.sub('', config.get(section, 'state'))
        except NoOptionError:
            state = 'Daily Ahope'
        if state in globalv.STATES:
            raise ValueError("State name for DailyAhopeTab must be unique, "
                             "please do not select '%s'" % state)
        globalv.STATES[state] = SummaryState(state, known=(start, end))
        globalv.STATES[state].definition = '%s:ahope' % ifo
        config.set(section, 'states', state)

        # parse generic configuration
        new = super(DailyAhopeTab, cls).from_ini(config, section, **kwargs)
        new.channel = re_quote.sub('', config.get(section, 'channel'))
        for p in new.plots + new.subplots:
            p.etg = new.name.lower()

        # work out day directory and url
        utc = from_gps(new.span[0])
        basedir = os.path.normpath(config.get(section, 'base-directory'))
        daydir = os.path.join(basedir, utc.strftime('%Y%m'),
                              utc.strftime('%Y%m%d'))
        home_, postbase = daydir.split('/public_html/', 1)
        user = os.path.split(home_)[1]
        new.ihopepage = '/~%s/%s/' % (user, postbase)

        # get cache options
        cachefile = config.get(section, 'inspiral-cache')
        new.inspiralcachefile = os.path.join(daydir,  cachefile)
        cachefile = config.get(section, 'tmpltbank-cache')
        new.tmpltbankcachefile = os.path.join(daydir, cachefile)
        segfile = config.get(section, 'segment-file')
        new.segmentfile = os.path.join(daydir, segfile)

        # get loudest options
        if config.has_option(section, 'loudest'):
            # set defaults
            new.loudest = {
                'N': config.getint(section, 'loudest'),
                'columns': ['end', 'new_snr', 'snr', 'mchirp',
                            'mass1', 'mass2', 'reduced_chisq'],
                'rank': 'new_snr',
                'dt': 8,
            }
            # override from config
            if config.has_option(section, 'loudest-columns'):
                new.loudest['columns'] = map(
                    lambda s: re_quote.sub('', s),
                    config.get(section, 'loudest-columns').split(','))
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
            if config.has_option(section, 'loudest-dt'):
                new.loudest['dt'] = config.getfloat(section, 'loudest-dt')
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
                try:
                    globalv.DATA[rangechannel].coalesce()
                except ValueError:
                    pass
                globalv.DATA[sizechannel].append(
                    TimeSeries(sizedata, sample_rate=1/dt, epoch=epoch,
                               name=sizechannel))
                try:
                    globalv.DATA[sizechannel].coalesce()
                except ValueError:
                    pass

    def process(self, *args, **kwargs):
        # read the segment files
        if os.path.isfile(self.segmentfile):
            segs = DataQualityFlag.read(self.segmentfile, coalesce=False)
            self.states[0].known = segs.known
            self.states[0].active = segs.active
            self.states[0].ready = True
        else:
            warn('Segment file %s not found.' % self.segmentfile)
            return
        if len(self.states[0].active) == 0:
            warn('No segments analysed by daily ahope.')
            return
        # read the cache files
        if os.path.isfile(self.inspiralcachefile):
            with open(self.inspiralcachefile, 'r') as fobj:
                try:
                    self.inspiralcache = Cache.fromfile(fobj).sieve(
                                             segment=self.span)
                except ValueError as e:
                    if "could not convert \'\\n\' to CacheEntry" in str(e):
                        self.inspiralcache = Cache()
                    else:
                        raise
        else:
            warn("Cache file %s not found." % self.inspiralcachefile)
            return
        if os.path.isfile(self.tmpltbankcachefile):
            with open(self.tmpltbankcachefile, 'r') as fobj:
                try:
                    self.tmpltbankcache = Cache.fromfile(fobj).sieve(
                                              segment=self.span)
                except ValueError:
                    if "could not convert \'\\n\' to CacheEntry" in str(e):
                        self.tmpltbankcache = Cache()
                    else:
                        raise
        else:
            warn("Cache file %s not found." % self.tmpltbankcachefile)
            self.tmpltbankcache = Cache()

        # only process if the cachfile was found
        super(DailyAhopeTab, self).process(*args, **kwargs)

    def process_state(self, state, nds=None, multiprocess=False,
                      config=GWSummConfigParser(),
                      segdb_error='raise', trigcache=None, datacache=None):
        if trigcache is None:
            trigcache = self.inspiralcache
        if datacache is None:
            datacache = Cache()
        super(DailyAhopeTab, self).process_state(
            state, nds=nds, multiprocess=multiprocess, config=config,
            datacache=datacache, trigcache=trigcache)

    def write_state_html(self, state):
        """Write the '#main' HTML content for this `DailyAhopeTab`.
        """
        daydir = os.path.split(self.segmentfile)[0]
        # did it run
        if not os.path.isdir(daydir):
            page = html.markup.page()
            page.div(class_='alert alert-warning')
            page.p("No analysis was performed for this period, "
                   "please try again later.")
            page.p("If you believe these data should have been found, please "
                   "contact %s."
                   % html.markup.oneliner.a('the CBC DQ group',
                                            class_='alert-link',
                                            href='mailto:cbc+dq@ligo.org'))
            page.div.close()
        elif (not os.path.isfile(self.segmentfile) or
              len(self.states[0].active) != 0 and
              not os.path.isfile(self.inspiralcachefile)):
            page = html.markup.page()
            page.div(class_='alert alert-danger')
            page.p("This analysis seems to have failed.")
            page.p("If you believe these data should have been found, please "
                   "contact %s."
                   % html.markup.oneliner.a('the CBC DQ group',
                                            class_='alert-link',
                                            href='mailto:cbc+dq@ligo.org'))
            page.div.close()
        elif len(self.states[0].active) == 0:
            page = html.markup.page()
            page.div(class_='alert alert-info')
            page.p("This analysis found no segments over which to run.")
            page.p("If you believe this to be an error, please contact %s."
                   % html.markup.oneliner.a('the CBC DQ group',
                                            class_='alert-link',
                                            href='mailto:cbc+dq@ligo.org'))
            page.div.close()
        else:
            # otherwise, carry on...
            page = self.scaffold_plots(state=state)

            # link full results
            page.hr(class_='row-divider')
            page.div(class_='btn-group')
            page.a('Click here for the full Daily Ahope results',
                   href=self.ihopepage, rel='external', target='_blank',
                   class_='btn btn-default btn-info btn-xl')
            page.div.close()
            page.hr(class_='row-divider')

            if self.loudest:
                table = get_triggers(self.channel, self.plots[0].etg, state,
                                     query=False)
                rank = get_table_column(
                    table, self.loudest['rank']).argsort()[::-1]
                loudest = []
                i = 0
                while len(loudest) < self.loudest['N'] and i < rank.size:
                    t = table[rank[i]]
                    if i == 0 or all([abs(float(t.get_end()) -
                                          float(t2.get_end())) >=
                                      self.loudest['dt'] for t2 in loudest]):
                        loudest.append(t)
                    i += 1
                page.h1('Loudest events')
                page.p('The following table displays the %d loudest events as '
                       'recorded by Daily Ahope (with at least %s-second '
                       'separation).' % (self.loudest['N'], self.loudest['dt']))
                headers = self.loudest['labels']
                if 'time' in headers[0]:
                    headers.insert(1, 'UTC time')
                    date = True
                else:
                    data = False
                data = []
                for row in loudest:
                    data.append([])
                    for column in self.loudest['columns']:
                        data[-1].append('%.3f' % float(get_row_value(row,
                                                                     column)))
                    if date:
                        data[-1].insert(1, from_gps(row.get_end()).strftime(
                                               '%B %d %Y, %H:%M:%S.%f')[:-3])
                page.add(str(html.data_table(headers, data, table='data')))

            if self.subplots:
                page.hr(class_='row-divider')
                page.h1('Sub-plots')
                layout = get_mode() == Mode.week and [7] or [4]
                plist = [p for p in self.subplots if p.state in [state, None]]
                page.add(str(self.scaffold_plots(plots=plist, state=state,
                                                 layout=layout)))

                # link full results
                page.hr(class_='row-divider')
                page.div(class_='btn-group')
                page.a('Click here for the full Daily Ahope results',
                       href=self.ihopepage, rel='external', target='_blank',
                       class_='btn btn-default btn-info btn-xl')
                page.div.close()
                page.hr(class_='row-divider')

        # write to file
        idx = self.states.index(state)
        with open(self.frames[idx], 'w') as fobj:
            fobj.write(str(page))
        return self.frames[idx]

register_tab(DailyAhopeTab)
register_tab(DailyAhopeTab, name='archived-daily-ihope')

register_reader('daily ahope', SnglInspiralTable,
                get_reader('ligolw', SnglInspiralTable))
