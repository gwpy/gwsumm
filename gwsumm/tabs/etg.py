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

"""Custom `SummaryTab` for the output of an ETG.
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

from .. import (version, html, globalv)
from ..config import (GWSummConfigParser, NoOptionError, DEFAULTSECT)
from ..data import get_channel
from ..triggers import (get_etg_table, get_triggers, register_etg_table)
from ..utils import re_quote
from ..state import SummaryState
from ..mode import (get_mode, MODE_ENUM)
from .registry import (get_tab, register_tab)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

base = get_tab('default')


class EventTriggerTab(base):
    """Custom tab displaying a summary of event trigger generation results.
    """
    type = 'archived-triggers'

    def __init__(self, name, start, end, channel=None, etg=None, table=None,
                 cachefile=None, url=None, **kwargs):
        super(EventTriggerTab, self).__init__(name, start, end, **kwargs)
        self.channel = channel and get_channel(channel) or None
        self.cachefile = cachefile
        self.url = url
        # parse ETG and LIGO_LW table class
        if etg is None:
            etg = self.name
        self.etg = etg
        if table is None or isinstance(table, str):
            tablename = isinstance(table, str) and table or self.etg
            try:
                table = get_etg_table(tablename)
            except KeyError as e:
                e.args = ("Cannot automatically determine LIGO_LW table for "
                          "etg %r, please specify in configuration file or "
                          "when creating EventTriggerTab" % tablename,)
                raise
        # register custom readers for this type
        try:
            register_etg_table(self.etg.lower(), table)
        except KeyError:
            pass
        try:
            register_reader(self.etg.lower(), table,
                            get_reader('ligolw', table))
        except Exception as e:
            if 'already defined' in str(e):
                pass
            else:
                raise

    @classmethod
    def from_ini(cls, config, section, **kwargs):
        """Define a new `EventTriggerTab` from a `ConfigParser`.
        """
        for key in ['channel', 'etg', 'url', 'cachefile', 'table']:
            try:
                kwargs.setdefault(key, config.get(section, key))
            except NoOptionError:
                pass
        new = super(EventTriggerTab, cls).from_ini(config, section, **kwargs)

        # set ETG for plots
        for p in new.plots + new.subplots:
            p.etg = new.etg.lower()

        # get loudest options
        if config.has_option(section, 'loudest'):
            # set defaults
            new.loudest = {
                'N': config.getint(section, 'loudest'),
                'rank': ['new_snr'],
                'dt': 8,
            }
            if 'hope' in new.etg.lower():
                new.loudest['columns'] = [
                    'end', 'new_snr', 'snr', 'mchirp', 'mass1', 'mass2',
                    'reduced_chisq',
                ]
            else:
                new.loudest['columns'] = [
                    'peak', 'duration', 'snr', 'peak_frequency', 'bandwidth',
                    'Q',
                ]
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
                new.loudest['rank'] = map(
                    lambda c: re_quote.sub('', c),
                    config.get(section, 'loudest-rank').split(','))
            if config.has_option(section, 'loudest-dt'):
                new.loudest['dt'] = config.getfloat(section, 'loudest-dt')
        else:
            new.loudest = None

        return new

    def process(self, *args, **kwargs):
        # read the cache files
        if self.cachefile is None:
            self.cache = None
        elif os.path.isfile(self.cachefile):
            with open(self.cachefile, 'r') as fobj:
                try:
                    self.cache = Cache.fromfile(fobj).sieve(
                                             segment=self.span)
                except ValueError as e:
                    if "could not convert \'\\n\' to CacheEntry" in str(e):
                        self.cache = Cache()
                    else:
                        raise
        else:
            warn("Cache file %s not found." % self.cachefile)
            return

        # only process if the cachfile was found
        if kwargs.get('trigcache', None) is None:
            kwargs['trigcache'] = self.cache
        super(EventTriggerTab, self).process(*args, **kwargs)

    def write_state_html(self, state):
        """Write the '#main' HTML content for this `EventTriggerTab`.
        """
        # did it run
        if self.cachefile and not os.path.isfile(self.cachefile):
            page = html.markup.page()
            page.div(class_='alert alert-danger')
            page.p("This analysis seems to have failed.")
            page.p("If you believe these data should have been found, please "
                   "contact %s."
                   % html.markup.oneliner.a('the DetChar group',
                                            class_='alert-link',
                                            href='mailto:detchar@ligo.org'))
            page.div.close()
        elif len(self.states[0].active) == 0:
            page = html.markup.page()
            page.div(class_='alert alert-info')
            page.p("This analysis found no segments over which to run.")
            page.p("If you believe this to be an error, please contact %s."
                   % html.markup.oneliner.a('the DetChar group',
                                            class_='alert-link',
                                            href='mailto:detchar@ligo.org'))
            page.div.close()
        else:
            # otherwise, carry on...
            page = self.scaffold_plots(state=state)

            # link full results
            if self.url:
                page.hr(class_='row-divider')
                page.div(class_='btn-group')
                page.a('Click here for the full %s results' % self.name,
                       href=self.url, rel='external', target='_blank',
                       class_='btn btn-default btn-info btn-xl')
                page.div.close()
                page.hr(class_='row-divider')

            if self.loudest:
                page.h1('Loudest events')
                page.p('The following table(s) displays the %d loudest events '
                       'as recorded by %s (with at least %s-second '
                       'separation).'
                       % (self.loudest['N'], self.etg, self.loudest['dt']))
                # get triggers
                table = get_triggers(self.channel, self.plots[0].etg, state,
                                     query=False)
                # set table headers
                headers = list(self.loudest['labels'])
                if 'time' in headers[0]:
                    headers.insert(1, 'UTC time')
                    date = True
                else:
                    date = False
                # loop over rank columns
                for rank in self.loudest['rank']:
                    try:
                        rankstr = self.loudest['labels'][
                            self.loudest['columns'].index(rank)]
                    except ValueError:
                        rankstr = repr(rank)
                    page.h3('Loudest events by %s' % rankstr)
                    rank = get_table_column(table, rank).argsort()[::-1]
                    loudest = []
                    i = 0
                    while len(loudest) < self.loudest['N'] and i < rank.size:
                        t = table[rank[i]]
                        if i == 0 or all([
                                abs(float(get_row_value(t, 'time')) -
                                    float(get_row_value(t2, 'time'))) >=
                                self.loudest['dt'] for t2 in loudest]):
                            loudest.append(t)
                        i += 1
                    data = []
                    for row in loudest:
                        data.append([])
                        for column in self.loudest['columns']:
                            data[-1].append(
                                '%.3f' % float(get_row_value(row, column)))
                        if date:
                            data[-1].insert(
                                1,
                                from_gps(get_row_value(
                                    row, self.loudest['columns'][0])).strftime(
                                    '%B %d %Y, %H:%M:%S.%f')[:-3])
                    page.add(str(html.data_table(headers, data, table='data')))

            if self.subplots:
                page.hr(class_='row-divider')
                page.h1('Sub-plots')
                layout = get_mode() == MODE_ENUM['WEEK'] and [7] or [4]
                plist = [p for p in self.subplots if p.state in [state, None]]
                page.add(str(self.scaffold_plots(plots=plist, state=state,
                                                 layout=layout)))

                # link full results
                if self.url:
                    page.hr(class_='row-divider')
                    page.div(class_='btn-group')
                    page.a('Click here for the full %s results' % self.name,
                           href=self.url, rel='external', target='_blank',
                           class_='btn btn-default btn-info btn-xl')
                    page.div.close()
                    page.hr(class_='row-divider')

        # write to file
        idx = self.states.index(state)
        with open(self.frames[idx], 'w') as fobj:
            fobj.write(str(page))
        return self.frames[idx]

register_tab(EventTriggerTab)
