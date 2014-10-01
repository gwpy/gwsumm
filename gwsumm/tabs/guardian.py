# coding=utf-8
# Copyright (C) Duncan Macleod (2014)
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

"""Definition of the `GuardianTab`
"""

from __future__ import print_function

import re

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from dateutil import tz

import numpy

from astropy.time import Time

from gwpy.segments import DataQualityDict

from ..config import GWSummConfigParser
from ..state import ALLSTATE
from .registry import (get_tab, register_tab)
from .. import (globalv, version, html)
from ..data import get_timeseries
from ..segments import get_segments
from ..plot.registry import (get_plot, register_plot)
from ..utils import vprint

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

Tab = get_tab('default')
UTC = tz.gettz('UTC')

class GuardianTab(Tab):
    """Summarises the data recorded by an Advanced LIGO Guardian node.

    Each guardian node controls and monitors state transitions for a
    specific subsystem of the Advanced LIGO interferometer. The
    `GuardianTab` summarises those data with a state segments plot, a
    transitions summary table, and a detailed list of transitions and
    segments for each listed state.
    """
    type = 'guardian'

    @classmethod
    def from_ini(cls, config, section, plotdir='plots', base=''):
        """Define a new `GuardianTab`.
        """
        new = super(Tab, cls).from_ini(
                  config, section, base=base)
        if len(new.states) > 1 or new.states[0].name != ALLSTATE:
            raise ValueError("GuardianTab does not accept state selection")
        new.plots = []
        new.plotdir = plotdir
        new.ifo = config.get(section, 'ifo')

        # record node and states
        new.node = config.get(section, 'node')
        new.grdstates = OrderedDict()
        for key, name in config.nditems(section):
            try:
                key = int(key)
            except ValueError:
                continue
            else:
                new.grdstates[int(key)] = name

        # build plots
        new.segmenttag = '%s:%s %%s' % (new.ifo, new.node)
        labels = new.grdstates.values()[::-1]
        flags = [new.segmenttag % name for name in labels]
        new.plots.append(get_plot('segments')(
            flags, new.span[0], new.span[1], labels=labels, outdir=plotdir,
            valid=None,
            title='%s Guardian %s state' % (
                new.ifo, new.node.replace('_', r'\_')),
            tag='GRD_%s_SEGMENTS' % re.sub('[-\s]', '_', new.node)))

        return new

    def process(self, nds='guess', multiprocess=True,
                config=GWSummConfigParser(), datacache=None,
                **kwargs):
        """Process data for the given state.
        """
        ifo = self.ifo

        for p in self.plots:
            if p.outputfile in globalv.WRITTEN_PLOTS:
                p.new = False

        # --------------------------------------------------------------------
        # work out which channels are needed

        prefix = '%s:GRD-%s_%%s' % (self.ifo, self.node)

        state = sorted(self.states, key=lambda s: abs(s.active))[0]
        grddata = get_timeseries(prefix % 'STATE_N', state, config=config,
                                 nds=nds, multiprocess=multiprocess,
                                 cache=datacache)
        vprint("    All time-series data loaded\n")

        # --------------------------------------------------------------------
        # find segments and transitions

        self.transitions = dict((v, []) for v in self.grdstates)

        for data in grddata:
            segs = DataQualityDict()
            for v, name in self.grdstates.iteritems():
                tag = self.segmenttag % name
                instate = data == v
                segs[tag] = instate.to_dqflag(name=name)
                for trans in (
                        numpy.diff(instate.astype(int)) == 1).nonzero()[0]:
                    t = data.times[trans+1]
                    from_ = data[trans].value
                    self.transitions[v].append((t, from_))

            globalv.SEGMENTS += segs

        super(GuardianTab, self).process(
            config=config, nds=nds, multiprocess=multiprocess,
            datacache=datacache, **kwargs)

    def write_state_html(self, state):
        """Write the HTML for the given state of this `GuardianTab`
        """
        page = self.scaffold_plots()
        page.div(class_='row')
        page.div(class_='col-md-12')

        # draw table
        page.h2('State transitions')
        page.p("The table below lists for each state (row), the number of "
               "transitions into that state from each other state (columns).")
        page.p("Only those states named in the configuration are shown, but "
               "the 'Total' includes transitions from any and all states.")
        page.table(class_='transitions data')
        page.tr(class_='header')
        for th in ['State'] + self.grdstates.values() + ['Total']:
            page.th(th)
        page.tr.close()
        for i, bit in enumerate(self.grdstates):
            page.tr()
            name = self.grdstates[bit]
            page.th(name)
            for j, bit2 in enumerate(self.grdstates):
                if i == j:
                    page.td('-', class_='IOP')
                    continue
                count = len([t for t in self.transitions[bit] if
                             t[1] == bit2])
                page.td(str(count))
            page.th(str(len(self.transitions[bit])))
            page.tr.close()
        page.table.close()
        page.div.close()
        page.div.close()

        # summarise each state
        page.div(class_='row')
        page.div(class_='col-md-12')
        page.h2('State details')
        page.div(class_='panel-group', id='accordion')
        for i, bit in enumerate(self.grdstates):
            name = self.grdstates[bit]
            page.div(class_='panel panel-default', id=str(bit))
            # heading
            page.a(href='#collapse%d' % bit,
                   **{'data-toggle': 'collapse', 'data-parent': '#accordion'})
            page.div(class_='panel-heading')
            page.h4(name, class_='panel-title')
            page.div.close()
            page.a.close()

            # body
            page.div(id='collapse%d' % bit, class_='panel-collapse collapse')
            page.div(class_='panel-body')

            # print transitions
            page.p('This state was entered %d times as follows:'
                   % len(self.transitions[bit]))
            headers = ['GPS time', 'UTC time', 'Local time', 'Transition from']
            data = []
            if self.ifo in ['H1', 'C1', 'P1']:
                localzone = tz.gettz('America/Los_Angeles')
            elif self.ifo in ['L1']:
                localzone = tz.gettz('America/Chicago')
            else:
                localzone = tz.gettz('Europe/Berlin')
            for t, from_ in self.transitions[bit]:
                t2 = Time(t, format='gps', scale='utc')
                tlocal = Time(
                    t2.datetime.replace(tzinfo=UTC).astimezone(localzone),
                    format='datetime', scale='utc')
                data.append((t, t2.iso, tlocal.iso, '%s [%d]' % (
                             self.grdstates.get(from_, 'Unknown'), from_)))
            page.add(str(html.data_table(headers, data, table='guardian data')))

            # print segments
            page.p('This state was active during the following segments:')
            flag = get_segments(self.segmenttag % name, state.active,
                                query=False).copy()
            page.add(str(self.print_segments(flag)))
            page.div.close()
            page.div.close()
            page.div.close()
        page.div.close()
        page.div.close()
        page.div.close()

        return super(Tab, self).write_state_html(state, plots=False,
                                                     pre=page)
register_tab(GuardianTab)

