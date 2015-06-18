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

from gwpy.time import from_gps
from gwpy.plotter.table import (get_table_column, get_row_value)

from .. import (version, html)
from ..config import NoOptionError
from ..data import get_channel
from ..triggers import (get_etg_table, get_triggers, register_etg_table)
from ..utils import re_quote
from ..mode import (get_mode, MODE_ENUM)
from .registry import (get_tab, register_tab)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


class EventTriggerTab(get_tab('default')):
    """Custom `DataTab` displaying a summary of event trigger generation

    Parameters
    ----------
    name : `str`
        name of this tab (required)
    start : `LIGOTimeGPS`, `str`
        start time of this `DataTab`, anything that can be parsed by
        `~gwpy.time.to_gps` is fine
    end : `LIGOTimeGPS`, `str`
        end time of this `DataTab`, format as for `start`
    channel : `str`
        name of the channel of interest
    etg : `str`, optional
        name of this event trigger generator (ETG), defaults to the name
        of this `EventTriggerTab`
    states : `list` of `states <gwsumm.state.SummaryState>`
        the `list` of states (`~gwsumm.state.SummaryState`) over which
        this `DataTab` should be processed. More states can be added
        later (but before running :meth:`~DataTab.process`) via
        :meth:`~DataTab.add_state`.
    table : `type`, `str`, optional
        LIGO_LW `~glue.ligolw.table.Table` class to use for this ETG,
        e.g. use `~glue.ligolw.lsctables.SnglBurstTable` for Omicron, or
        `~glue.ligolw.lsctables.SnglInspiralTable` for CBC
    cache : `~glue.lal.Cache`, `str`, optional
        `Cache` object, or path to a LAL-format cache file on disk,
        from which to read the event triggers. If no cache is given,
        the ~gwpy.table.io.trigfind` module will be used to automatically
        locate the trigger files.
    url : `str`, optional
        URL for linking to more details results for this tab.
    **kwargs
        all other keyword arguments accepted by the `DataTab`

    See Also
    --------
    gwsumm.tabs.DataTab
        for details on the other keyword arguments (``**kwargs``)
        accepted by the constructor for this `EventTriggerTab`.
    """
    type = 'archived-triggers'

    def __init__(self, name, start, end, channel=None, etg=None, table=None,
                 cache=None, url=None, **kwargs):
        """Create a new `EventTriggerTab`
        """
        super(EventTriggerTab, self).__init__(name, start, end, **kwargs)
        self.channel = channel and get_channel(channel) or None
        self.cache = cache
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

        Parameters
        ----------
        cp : :class:`~gwsumm.config.GWConfigParser`
            customised configuration parser containing given section
        section : `str`
            name of section to parse
        *args, **kwargs
            other positional and keyword arguments to pass to the class
            constructor (`__init__`)

        Notes
        -----
        In addition to those attributes parsed by the parent
        `DataTab.from_ini` method, this method will parse the following
        attributes from a `ConfigParser`::

        .. autosummary::

           ~EventTriggerTab.channel
           ~EventTriggerTab.etg
           ~EventTriggerTab.url
           ~EventTriggerTab.cache
           ~EventTriggerTab.table

        Additionally, the loudest events table is configured by giving the
        following options

        -----------------  ---------------------------------------------
        `loudest`          the number of events to include in the table
        `loudest-rank`     the statistic to use in sorting the table
        `loudest-dt`       the minimum time separation for unique events
        `loudest-columns`  the columns to print in the table
        `loudest-labels`   the labels for each of the given columns
        -----------------  ---------------------------------------------

        See Also
        --------
        DataTab.from_ini
            the parent parsing method that handles parsing plot defintions,
            amongst other things
        """
        for key in ['channel', 'etg', 'url', 'cache', 'table']:
            try:
                kwargs.setdefault(
                    key, re_quote.sub('', config.get(section, key)))
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
        if isinstance(self.cache, str) and os.path.isfile(self.cache):
            with open(self.cache, 'r') as fobj:
                try:
                    self.cache = Cache.fromfile(fobj).sieve(
                                             segment=self.span)
                except ValueError as e:
                    if "could not convert \'\\n\' to CacheEntry" in str(e):
                        self.cache = Cache()
                    else:
                        raise
        elif isinstance(self.cache, str):
            warn("Cache file %s not found." % self.cache)
            return
        elif self.cache is not None and not isinstance(self.cache, Cache):
            raise ValueError("Cannot parse EventTriggerTab.cache of type %r"
                             % type(self.cache))

        # only process if the cachfile was found
        if kwargs.get('trigcache', None) is None:
            kwargs['trigcache'] = self.cache
        super(EventTriggerTab, self).process(*args, **kwargs)

    def write_state_html(self, state):
        """Write the '#main' HTML content for this `EventTriggerTab`.
        """
        # did it run
        if not self.cache:
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
