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

import re
from collections import OrderedDict
from configparser import NoOptionError

from dateutil import tz

import numpy

from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize

from astropy.time import Time

from glue.lal import Cache

from gwpy.segments import DataQualityDict

from gwdetchar.io import html

from .. import globalv
from ..config import GWSummConfigParser
from ..data import get_timeseries_dict
from ..plot.registry import get_plot
from ..segments import get_segments
from ..state import ALLSTATE
from ..utils import vprint
from .registry import (get_tab, register_tab)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['GuardianTab']

DataTab = get_tab('default')
UTC = tz.gettz('UTC')
REQUESTSTUB = '+request'
NOMINALSTUB = '+nominal'
MODE_COLORS = ['grey', 'magenta', 'red', 'saddlebrown']

re_guardian_index = re.compile(r'\[(?P<idx>.*)\] (?P<label>.*)')


class GuardianTab(DataTab):
    """Summarises the data recorded by an Advanced LIGO Guardian node.

    Each guardian node controls and monitors state transitions for a
    specific subsystem of the Advanced LIGO interferometer. The
    `GuardianTab` summarises those data with a state segments plot, a
    transitions summary table, and a detailed list of transitions and
    segments for each listed state.
    """
    type = 'guardian'

    @classmethod
    def from_ini(cls, config, section, plotdir='plots', **kwargs):
        """Define a new `GuardianTab`.
        """
        new = super(DataTab, cls).from_ini(config, section, **kwargs)
        new.plots = []
        new.plotdir = plotdir
        new.ifo = config.get(section, 'IFO')

        # record node and states
        new.node = config.get(section, 'node')
        new.grdstates = OrderedDict()
        plot = []
        for key, name in config.nditems(section):
            try:
                key = int(key)
            except ValueError:
                continue
            else:
                if name[0] == '*':
                    plot.append(False)
                else:
                    plot.append(True)
                new.grdstates[int(key)] = name.strip('*')
        try:
            new.transstates = map(
                int, config.get(section, 'transitions').split(','))
        except NoOptionError:
            new.transstates = list(new.grdstates)

        # -- build plots ------------------------

        new.set_layout([1, 2])
        grdidxs = dict((state, idx) for idx, state in
                       new.grdstates.items())
        new.segmenttag = '%s:%s %%s' % (new.ifo, new.node)
        pstates = [l for i, l in enumerate(list(new.grdstates.values())[::-1])
                   if plot[-i-1]]
        flags = [new.segmenttag % name for name in pstates]
        labels = ['[%d] %s' % (grdidxs[state], state) for state in pstates]

        # define colours
        cmap = get_cmap('brg')(Normalize(-1, 1)(
            numpy.linspace(-1, 1, len(pstates))))[::-1]
        th = len(flags) > 8 and (new.span[1] - new.span[0])/200. or 0

        for state in new.states:
            # get common plot tag prefix
            tagprefix = 'GRD_%s' % re.sub(r'[-\s]', '_', new.node.upper())
            if state.name.lower() != ALLSTATE:  # include state name if not All
                tagprefix = '%s_%s' % (state.name, tagprefix)
            # segment plot
            new.plots.append(get_plot('guardian')(
                flags, new.span[0], new.span[1], state=state,
                labels=labels, outdir=plotdir,
                known={'hatch': 'x', 'alpha': 0.1, 'facecolor': 'none'},
                tag='%s_SEGMENTS' % tagprefix,
                title='%s Guardian %s state' % (
                    new.ifo, new.node.replace('_', r'\_')), zorder=2))

            # pie
            new.plots.append(get_plot('segment-pie')(
                flags, new.span[0], new.span[1], state=state,
                labels=pstates, colors=cmap,
                tag='%s_SEGMENT_PIE' % tagprefix,
                startangle=180, counterclock=False, wedge_linewidth=0.01,
                outdir=plotdir, title='%s Guardian %s state' % (
                    new.ifo, new.node.replace('_', r'\_')),
                legend_fontsize=16, legend_sorted=True, legend_threshold=th,
            ))

            # bar
            new.plots.append(get_plot('segment-bar')(
                flags, new.span[0], new.span[1], state=state,
                labels=pstates, sorted=True,
                tag='%s_SEGMENT_BAR' % tagprefix,
                outdir=plotdir, title='%s Guardian %s state' % (
                    new.ifo, new.node.replace('_', r'\_')),
            ))
        return new

    def process(self, nds=None, nproc=1,
                config=GWSummConfigParser(), datacache=None,
                segmentcache=Cache(), datafind_error='raise', **kwargs):
        """Process data for the given state.
        """
        # finalize state information
        self.finalize_states(
            config=config, segdb_error=kwargs.get('segdb_error', 'raise'),
            datafind_error=datafind_error)
        vprint("States finalised [%d total]\n" % len(self.states))
        vprint("    Default state: %r\n" % str(self.defaultstate))

        # remove plots that have already been generated
        for p in self.plots:
            if p.outputfile in globalv.WRITTEN_PLOTS:
                p.new = False

        # --------------------------------------------------------------------
        # work out which channels are needed

        prefix = '%s:GRD-%s_%%s' % (self.ifo, self.node)
        state = sorted(self.states, key=lambda s: abs(s.active))[-1]

        prefices = ['STATE_N', 'REQUEST_N', 'NOMINAL_N', 'OK', 'MODE', 'OP']
        alldata = list(get_timeseries_dict(
            [prefix % x for x in prefices],
            state, config=config, nds=nds, nproc=nproc,
            cache=datacache, datafind_error=datafind_error,
            dtype='int32').values())
        vprint("    All time-series data loaded\n")

        # --------------------------------------------------------------------
        # find segments and transitions

        self.transitions = dict((v, []) for v in self.grdstates)

        for sdata, rdata, ndata, okdata in zip(*alldata[:4]):
            ssegs = DataQualityDict()
            rsegs = DataQualityDict()
            nsegs = DataQualityDict()
            oksegs = (okdata == 1).to_dqflag(name='Node OK')
            for v, name in self.grdstates.items():
                # get segments for state
                tag = self.segmenttag % name
                instate = sdata == v
                ssegs[tag] = instate.to_dqflag(name=name)
                diff_ = numpy.diff(instate.value.astype(int))
                transin = (diff_ == 1).nonzero()[0] + 1
                transout = (diff_ == -1).nonzero()[0] + 1
                for i, j in zip(transin, transout):
                    t = sdata.times[i].value
                    from_ = sdata[i-1].value
                    to_ = sdata[j].value
                    self.transitions[v].append((t, from_, to_))
                # get segments for request
                tag = self.segmenttag % name + REQUESTSTUB
                instate = rdata == v
                rsegs[tag] = instate.to_dqflag(name=name)
                # get segments for nominal
                tag = self.segmenttag % name + NOMINALSTUB
                nom = ndata == v
                nsegs[tag] = nom.to_dqflag(name=name)

            globalv.SEGMENTS += ssegs
            globalv.SEGMENTS += rsegs
            globalv.SEGMENTS += nsegs
            globalv.SEGMENTS += {self.segmenttag % 'OK': oksegs}

        super(GuardianTab, self).process(
            config=config, nds=nds, nproc=nproc,
            datacache=datacache, segmentcache=segmentcache, **kwargs)

    def write_state_html(self, state):
        """Write the HTML for the given state of this `GuardianTab`
        """
        page = self.scaffold_plots(state=state)
        page.div(class_='alert alert-info alert-dismissible')
        # add close button
        page.button(type="button", class_="close", **{'data-dismiss': "alert"})
        page.span('&times;', **{'aria-hidden': "true"})
        page.span('Close', class_="sr-only")
        page.button.close()
        # show message
        page.p('For all of the following data, "Unknown" simply labels any '
               'state in this node that was not chosen for display, and does '
               'not mean that the state was unrecognised by the Guardian '
               'system. Transitions from an "Unkown" state are not listed in '
               'the below table, but are included in the totals.')
        page.div.close()

        page.div(class_='row')
        page.div(class_='col-md-12')

        # draw table
        page.h2('%s state transitions' % self.node)
        id_ = '{}-state-transitions'.format(self.ifo.lower())
        page.table(class_='table table-condensed table-hover '
                          'table-responsive transitions', id_=id_)
        page.caption("Transitions into each state (row) from each other "
                     "state (column). Only those states named in the "
                     "configuration are shown, but the 'Total' includes "
                     "transitions from any and all states. '-' indicates no "
                     "transitions from that state.")
        page.thead()
        page.tr()
        for th in ['State'] + list(self.grdstates.values()) + ['Total']:
            page.th(th)
        page.tr.close()
        page.thead.close()
        page.tbody()
        for i, bit in enumerate(self.transstates):
            page.tr()
            name = self.grdstates[bit].strip('*')
            page.th(name)
            for j, bit2 in enumerate(self.grdstates):
                if i == j:
                    page.td('-', class_='ignore')
                    continue
                count = len([t for t in self.transitions[bit] if
                             t[1] == bit2])
                if count:
                    page.td(str(count))
                else:
                    page.td('-')
            page.th(str(len(self.transitions[bit])))
            page.tr.close()
        page.tbody.close()
        page.table.close()
        page.button(
            'Export to CSV', class_='btn btn-default btn-table',
            onclick="exportTableToCSV('{name}.csv', '{name}')".format(
                name=id_))
        page.div.close()
        page.div.close()

        # summarise each state
        page.div(class_='row')
        page.div(class_='col-md-12')
        page.h2('State details')
        page.div(class_='panel-group', id='accordion')
        for i, bit in enumerate(self.grdstates):
            name = self.grdstates[bit].strip('*')
            page.div(class_='panel well panel-primary', id=str(bit))
            # heading
            page.div(class_='panel-heading')
            page.a(href='#collapse-%d' % bit,
                   **{'data-toggle': 'collapse', 'data-parent': '#accordion'})
            page.h4('%s [%d]' % (name, bit), class_='panel-title')
            page.a.close()
            page.div.close()

            # body
            page.div(id='collapse-%d' % bit, class_='panel-collapse collapse')
            page.div(class_='panel-body')

            # print transitions
            headers = ['Transition GPS time', 'UTC time', 'Local time',
                       'Transition from', 'Exited to']
            data = []
            if self.ifo in ['H1', 'C1', 'P1']:
                localzone = tz.gettz('America/Los_Angeles')
            elif self.ifo in ['L1']:
                localzone = tz.gettz('America/Chicago')
            else:
                localzone = tz.gettz('Europe/Berlin')
            for t, from_, to_ in self.transitions[bit]:
                t2 = Time(t, format='gps', scale='utc')
                tlocal = Time(
                    t2.datetime.replace(tzinfo=UTC).astimezone(localzone),
                    format='datetime', scale='utc')
                data.append((
                    t, t2.iso, tlocal.iso,
                    '%s [%d]' % (self.grdstates.get(from_, 'Unknown'), from_),
                    '%s [%d]' % (self.grdstates.get(to_, 'Unknown'), to_)))
            page.add(str(html.table(
                headers, data,
                id='%s-guardian-%s' % (self.ifo.lower(), str(bit)),
                caption="Transitions for %s %r state" % (self.node, name))))

            # print segments
            flag = get_segments(self.segmenttag % name, state.active,
                                query=False).copy()
            livetime = abs(flag.active)
            try:
                duty = abs(flag.active) / float(abs(flag.known)) * 100.
            except ZeroDivisionError:
                duty = 0
            page.p('This state was active for %.2f seconds (%.2f%%) during '
                   'the following segments:' % (livetime, duty))
            page.add(str(self.print_segments(flag)))
            page.div.close()  # panel-body
            page.div.close()  # panel-collapse
            page.div.close()  # panel
        page.div.close()
        page.div.close()
        page.div.close()

        return super(DataTab, self).write_state_html(state, plots=False,
                                                     pre=page)


register_tab(GuardianTab)
