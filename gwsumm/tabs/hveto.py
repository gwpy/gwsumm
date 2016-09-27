# coding=utf-8
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

"""Custom `SummaryTab` for the output of the HierarchicalVeto algorithm.
"""

import os
import re
from glob import glob

from numpy import loadtxt

from astropy.io.registry import register_reader

from glue.lal import (Cache, CacheEntry, LIGOTimeGPS)

from gwpy.table import lsctables
from gwpy.segments import (SegmentList, DataQualityFlag)

from .registry import (get_tab, register_tab)

from .. import (html, globalv)
from ..mode import Mode
from ..config import (GWSummConfigParser, NoSectionError, NoOptionError,
                      DEFAULTSECT)
from ..channels import get_channel
from ..segments import get_segments
from ..state import SummaryState
from ..triggers import (register_etg_table, get_triggers, add_triggers)
from ..plot import (get_plot, register_plot)
from ..utils import re_quote

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['HvetoTab']

base = get_tab('default')
SummaryPlot = get_plot(None)
SegmentPlot = get_plot('segments')

HVETO_COLUMNS = ['peak_time', 'peak_time_ns', 'peak_frequency', 'snr']


class HvetoTab(base):
    """Custom tab displaying a summary of Hveto results.
    """
    type = 'hveto'
    summaryrows = ['Winning channel', 'Time Window [s]', 'SNR Thresh.',
                   'Significance', 'N. trigs',
                   'Use %', 'Efficiency [%]', 'Deadtime [%]',
                   'Cum. Efficiency [%]', 'Cum. Deadtime [%]']

    def __init__(self, *args, **kwargs):
        if kwargs['mode'] != Mode.day:
            raise RuntimeError("HvetoTab is only available in %s mode."
                               % Mode.day.name)
        super(HvetoTab, self).__init__(*args, **kwargs)

    @classmethod
    def from_ini(cls, config, section, **kwargs):
        """Define a new `HvetoTab` from a `ConfigParser`.
        """
        # set state information
        start = config.getint(DEFAULTSECT, 'gps-start-time')
        end = config.getint(DEFAULTSECT, 'gps-end-time')
        ifo = config.get('DEFAULT', 'ifo')

        if config.has_option(section, 'states'):
            raise ValueError("HvetoTab does not support configuration of "
                             "multiple states, please use the 'state' option "
                             "to name the Hveto state")
        try:
            state = re_quote.sub('', config.get(section, 'state'))
        except NoOptionError:
            state = 'Hveto'
        if state in globalv.STATES:
            raise ValueError("State name for HvetoTab must be unique, "
                             "please do not select '%s'" % state)
        globalv.STATES[state] = SummaryState(state, known=(start, end))
        globalv.STATES[state].definition = '%s:hveto' % ifo
        config.set(section, 'states', state)

        # create hveto ETG section
        try:
            config.get('hveto', 'columns')
        except NoSectionError:
            config.add_section('hveto')
            config.set('hveto', 'columns', ','.join(HVETO_COLUMNS))
        except NoOptionError:
            pass

        # parse generic configuration
        new = super(HvetoTab, cls).from_ini(config, section, **kwargs)

        # work out day directory and url
        basedir = os.path.normpath(config.get(section, 'base-directory'))
        daydir = os.path.realpath(
            os.path.join(basedir, config.get(section, 'directory-tag')))
        new.directory = daydir
        home_, postbase = daydir.split('/public_html/', 1)
        user = os.path.split(home_)[1]
        new.url = '/~%s/%s/' % (user, postbase)
        return new

    def process(self, config=GWSummConfigParser(), **kwargs):

        # set params
        self.rounds = None

        if not os.path.isdir(self.directory):
            self.rounds = None
            return

        # get some basic info
        ifo = config.get('DEFAULT', 'ifo')

        # read the configuration
        d = os.path.realpath(self.directory).rstrip('/')
        self.conf = dict()
        confs = glob(os.path.join(d, '%s-HVETO_CONF-*-*.txt' % ifo))
        if len(confs) != 1:
            self.rounds = 'FAIL'
            return
        conffile = confs[0]
        try:
            with open(conffile) as f:
                self.conf = dict()
                lines = f.readlines()[3:]
                for line in lines:
                    try:
                        key, val = line.split(': ', 1)
                        self.conf[key.strip()] = eval(val)
                    except (ValueError, SyntaxError, NameError):
                        pass
        except IOError:
            self.rounds = 'FAIL'
            return
        else:
            etg = self.conf.pop('AUXtype', None)
            if 'DEfnm' in self.conf:
                name = re_quote.sub('', self.conf['DEfnm'])
                self.primary = '%s:%s' % (ifo, name)
                if 'DEtype' in self.conf:
                    hetg = re_quote.sub('', self.conf['DEtype'])
                    if re.search('_%s\Z' % hetg, self.primary, re.I):
                        self.primary = self.primary[:-len(hetg)-1]
            else:
                self.primary = None

        # find the segments
        try:
            ce = CacheEntry.from_T050017(conffile)
        except (ValueError):
            start = int(self.span[0])
            duration = int(abs(self.span))
            span = self.span
        else:
            start = int(ce.segment[0])
            duration = int(abs(ce.segment))
            span = ce.segment
        try:
            statefile = self.conf['dqfnm']
        except KeyError:
            statefile = '%s-HVETO_DQ_SEGS-%d-%d.txt' % (ifo, start, duration)
        if not os.path.isfile(os.path.join(self.directory, statefile)):
            self.rounds = 'NOSEGMENTS'
            return

        # find the results table
        resultsfile = os.path.join(self.directory, 'summary_stats.txt')
        if not os.path.isfile(resultsfile):
            self.rounds = 'FAIL'
            return

        # determine the Hveto state
        cache = Cache([CacheEntry.from_T050017(
                           os.path.join(self.directory, statefile))])
        segments = SegmentList.read(cache)
        globalv.SEGMENTS[self.states[0].definition] = DataQualityFlag(
            self.states[0].definition, known=[span], active=segments)
        self.finalize_states(config=config, query=False)

        # read results file
        self.rounds = []
        with open(resultsfile, 'r') as f:
            for line in f.readlines():
                self.rounds.append(dict(zip(self.summaryrows,
                                            line.split(' ')[1:])))
                # fix channel name
                c = '%s:%s' % (ifo, self.rounds[-1]['Winning channel'])
                if etg and re.search('_%s\Z' % etg, c, re.I):
                     c = c.rsplit('_', 1)[0]
                self.rounds[-1]['Winning channel'] = c

        # read starting triggers
        rawfile = ('%s-HVETO_RAW_TRIGS_ROUND_0-%d-%d.txt'
                   % (ifo, start, duration))
        cache = Cache([CacheEntry.from_T050017(
                           os.path.join(self.directory, rawfile))])
        get_triggers('%s:hveto_start' % ifo, 'hveto', [self.span],
                     config=config, cache=cache, return_=False)
        get_triggers('%s:hveto_vetoed_all' % ifo, 'hveto', [self.span],
                     config=config, cache=Cache(), return_=False)

        for r in range(1, len(self.rounds) + 1):
            # read round veto triggers
            rawfile = ('%s-HVETO_VETOED_TRIGS_ROUND_%d-%d-%d.txt'
                       % (ifo, r, start, duration))
            cache = Cache([CacheEntry.from_T050017(
                               os.path.join(self.directory, rawfile))])
            trigs = get_triggers('%s:hveto_vetoed_round %d' % (ifo, r), 'hveto',
                         [self.span], config=config, cache=cache)
            add_triggers(trigs, '%s:hveto_vetoed_all,hveto' % ifo,
                         segments=SegmentList([self.span]))
            # read round veto segments
            segfile = ('%s-HVETO_VETO_SEGS_ROUND_%d-%d-%d.txt'
                       % (ifo, r, start, duration))
            cache = Cache([CacheEntry.from_T050017(
                               os.path.join(self.directory, segfile))])
            get_segments('%s:hveto_veto_segs_round_%d' % (ifo, r), [self.span],
                         config=config, cache=cache, return_=False)

        for plot in self.plots:
            if isinstance(plot, HvetoSegmentSummaryPlot):
                plot.find_flags()

        kwargs['trigcache'] = Cache()
        kwargs['segmentcache'] = Cache()
        super(HvetoTab, self).process(config=config, **kwargs)

        # find some plots
        for plot in ['OVERAL_HISTOGRAM', 'OVERAL_EFF_DT'][::-1]:
             filename = (
                 '%s-HVETO_%s-%d-%d.png' % (ifo, plot, start, duration))
             plotfile = os.path.join(self.directory, filename)
             if os.path.isfile(plotfile):
                 p = SummaryPlot(os.path.join(self.url, filename), new=False)
                 p.state = self.states[0]
                 self.plots.insert(0, p)

        # delete data from archive
        del globalv.SEGMENTS[self.states[0].definition]
        for row in range(1, len(self.rounds) + 1):
            del globalv.SEGMENTS['%s:hveto_veto_segs_round_%s' % (ifo, row)]

    def write_state_html(self, state):
        """Write the '#main' HTML content for this `HvetoTab`.
        """
        page = html.markup.page()

        # run failed
        url = html.markup.oneliner.a("This analysis", href=self.url,
                                     class_='alert-link')
        if self.rounds is None:
            page.div(class_='alert alert-info')
            page.p("%s has not been performed." % url)
            page.p("If you believe this represents a problem, please "
                   "contact %s."
                   % html.markup.oneliner.a('the DetChar group',
                                            class_='alert-link',
                                            href='mailto:detchar@ligo.org'))
            page.div.close()
        elif self.rounds == 'FAIL':
            page.div(class_='alert alert-warning')
            page.p("%s seems to have failed." % url)
            page.p("If you believe these data should have been found, please "
                   "contact %s."
                   % html.markup.oneliner.a('the DetChar group',
                                            class_='alert-link',
                                            href='mailto:detchar@ligo.org'))
            page.div.close()
        # no segments
        elif self.rounds == 'NOSEGMENTS':
            page.div(class_='alert alert-info')
            page.p("%s found no segments." % url)
            page.p("If you believe some segments should have been found, "
                   "please contact %s."
                   % html.markup.oneliner.a('the DetChar group',
                                            class_='alert-link',
                                            href='mailto:detchar@ligo.org'))
            page.div.close()
        # otherwise...
        else:
            # print results table
            page.div(class_='scaffold well')
            if self.primary:
                channel = get_channel(self.primary)
                page.p()
                page.strong('Primary channel: ')
                if channel.url:
                    page.add(html.markup.oneliner.a(self.primary,
                                                    href=channel.url,
                                                    target='_blank'))
                else:
                    page.add(self.primary)
                page.p.close()
            headers = list(self.summaryrows)
            data = []
            for i, round in enumerate(self.rounds):
                data.append([str(i + 1)] + [str(round[key]) for key in headers])
                channel = get_channel(data[-1][1])
                # format CIS url and type
                if re.search('\.[a-z]+\Z', channel.name):
                    name, ctype = channel.name.rsplit('.', 1)
                    c2 = get_channel(name)
                    ctype = ctype in ['rms'] and ctype.upper() or ctype.title()
                else:
                    c2 = channel
                    ctype = 'Raw'
                if c2.url:
                    data[-1][1] = html.markup.oneliner.a(str(channel),
                                                  href=c2.url,
                                                  target='_blank')
                else:
                    data[-1][1] = str(channel)
            page.add(str(html.data_table(['Round'] + headers, data)))
            page.div.close()

            # add plots
            page.hr(class_='row-divider')
            page.add(str(self.scaffold_plots(state=state)))
            page.hr(class_='row-divider')

            # link full results
            page.div(class_='btn-group')
            page.a('Click here for the full Hveto results',
                   href=self.url, rel='external', target='_blank',
                   class_='btn btn-default btn-info btn-xl')
            page.div.close()

        # write to file
        idx = self.states.index(state)
        with open(self.frames[idx], 'w') as fobj:
            fobj.write(str(page))
        return self.frames[idx]

register_tab(HvetoTab)
register_etg_table('hveto', 'sngl_burst')


def read_hveto_triggers(f, columns=HVETO_COLUMNS, filt=None, nproc=1):
    """Read a `SnglBurstTable` of triggers from an Hveto txt file.
    """
    # allow multiprocessing
    if nproc != 1:
        from gwpy.table.io.cache import read_cache
        return read_cache(f, lsctables.SnglBurstTable.tableName,
                          columns=columns, nproc=nproc, format='hveto')

    # format list of files
    if isinstance(f, CacheEntry):
        files = [f.path]
    elif isinstance(f, (str, unicode)) and f.endswith(('.cache', '.lcf')):
        files = open_cache(f).pfnlist()
    elif isinstance(f, (str, unicode)):
        files = f.split(',')
    elif isinstance(f, Cache):
        files = f.pfnlist()
    else:
        files = list(f)

    # generate output
    out = lsctables.New(lsctables.SnglBurstTable, columns=columns)
    append = out.append

    # iterate over files
    for f in files:
        trigs = loadtxt(f, dtype=float)
        for t, f, snr in trigs:
            b = lsctables.SnglBurst()
            b.set_peak(LIGOTimeGPS(float(t)))
            b.peak_frequency = f
            b.snr = snr
            if filt is None or filt(b):
                append(b)
    return out

register_reader('hveto', lsctables.SnglBurstTable, read_hveto_triggers)


class HvetoSegmentSummaryPlot(SegmentPlot):
    """Custom SegmentSummaryPlot to handle unkown numbers of hveto rounds
    """
    type = 'hveto-segments'
    defaults = SegmentPlot.defaults.copy()
    defaults.update({
        'on_is_bad': True,
        'valid': None,
        'insetlabels': False,
    })

    def find_flags(self):
        # work out flags on-the-fly
        if not self.flags:
            tag = 'hveto_veto_segs_round_1'
            flag = None
            for key in globalv.SEGMENTS:
                if key.endswith(tag):
                    flag = key
                    break
            i = 1
            while True:
                if flag in globalv.SEGMENTS:
                    self.add_flag(flag)
                else:
                    break
                i += 1
                flag = '%s_%s' % (flag.rsplit('_', 1)[0], str(i))
            self.pargs.setdefault('labels', ['Round %d' % (i+1) for
                                             i in range(len(self.flags))])

    def init_plot(self, **kwargs):
        plot, axes = super(HvetoSegmentSummaryPlot, self).init_plot(**kwargs)
        for ax in axes:
            ax.tick_params(axis='y', labelsize=18)
        return plot, axes

register_plot(HvetoSegmentSummaryPlot)
