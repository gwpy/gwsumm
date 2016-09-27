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

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from dateutil import tz

import numpy

from matplotlib.cm import get_cmap
from matplotlib.colors import Normalize

from astropy.time import Time

from glue.lal import Cache

from gwpy.segments import (DataQualityDict, SegmentList, Segment)
from gwpy.plotter import SegmentPlot

from ..config import (GWSummConfigParser, NoOptionError)
from ..state import ALLSTATE
from .registry import (get_tab, register_tab)
from .. import (globalv, html)
from ..data import (get_timeseries, get_timeseries_dict)
from ..segments import get_segments
from ..plot.registry import (get_plot, register_plot)
from ..utils import (vprint, re_quote)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['GuardianTab']

DataTab = get_tab('default')
UTC = tz.gettz('UTC')
REQUESTSTUB = '+request'
NOMINALSTUB = '+nominal'
MODE_COLORS = ['grey', 'magenta', 'red', 'saddlebrown']

re_guardian_index = re.compile('\[(?P<idx>.*)\] (?P<label>.*)')


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
        if len(new.states) > 1 or new.states[0].name != ALLSTATE:
            raise ValueError("GuardianTab does not accept state selection")
        new.plots = []
        new.plotdir = plotdir
        new.ifo = config.get(section, 'ifo')

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
            new.transstates = new.grdstates.keys()

        # -- build plots ------------------------

        new.set_layout([1, 2])
        grdidxs = dict((state, idx) for idx, state in
                       new.grdstates.iteritems())
        new.segmenttag = '%s:%s %%s' % (new.ifo, new.node)
        pstates = [l for i, l in enumerate(new.grdstates.values()[::-1])
                  if plot[-i-1]]
        flags = [new.segmenttag % name for name in pstates]
        labels = ['[%d] %s' % (grdidxs[state], state) for state in pstates]

        # segment plot
        new.plots.append(get_plot('guardian')(
            flags, new.span[0], new.span[1], labels=labels, outdir=plotdir,
            known={'hatch': 'x', 'alpha': 0.1, 'facecolor': 'none'},
            tag='GRD_%s_SEGMENTS' % re.sub('[-\s]', '_', new.node),
            title='%s Guardian %s state' % (
                new.ifo, new.node.replace('_', r'\_')), zorder=2))

        # pie
        cmap = get_cmap('brg')(Normalize(-1, 1)(
            numpy.linspace(-1, 1, len(pstates))))[::-1]
        th = len(flags) > 8 and (new.span[1] - new.span[0])/200. or 0
        new.plots.append(get_plot('segment-pie')(
            flags, new.span[0], new.span[1], labels=pstates, colors=cmap,
            tag='GRD_%s_SEGMENT_PIE' % re.sub('[-\s]', '_', new.node),
            startangle=180, counterclock=False, wedge_linewidth=0.01,
            outdir=plotdir, title='%s Guardian %s state' % (
                new.ifo, new.node.replace('_', r'\_')),
            legend_fontsize=16, legend_sorted=True, legend_threshold=th,
        ))

        # bar
        new.plots.append(get_plot('segment-bar')(
            flags, new.span[0], new.span[1], labels=pstates, sorted=True,
            tag='GRD_%s_SEGMENT_BAR' % re.sub('[-\s]', '_', new.node),
            outdir=plotdir, title='%s Guardian %s state' % (
                new.ifo, new.node.replace('_', r'\_')),
        ))
        return new

    def process(self, nds=None, multiprocess=True,
                config=GWSummConfigParser(), datacache=None,
                segmentcache=Cache(), datafind_error='raise', **kwargs):
        """Process data for the given state.
        """
        for p in self.plots:
            if p.outputfile in globalv.WRITTEN_PLOTS:
                p.new = False

        # --------------------------------------------------------------------
        # work out which channels are needed

        prefix = '%s:GRD-%s_%%s' % (self.ifo, self.node)
        state = sorted(self.states, key=lambda s: abs(s.active))[0]

        try:
            version = get_timeseries(
                prefix % 'VERSION', state, config=config, nds=nds,
                multiprocess=multiprocess, cache=datacache,
                datafind_error=datafind_error,
            ).join(gap='ignore').min().value
        except ValueError:
            version = 1201

        prefices = ['STATE_N', 'REQUEST_N', 'NOMINAL_N', 'OK', 'MODE']
        if version >= 1200:
            prefices.append('OP')
        alldata = get_timeseries_dict(
            [prefix % x for x in prefices],
            state, config=config, nds=nds, multiprocess=multiprocess,
            cache=datacache, datafind_error=datafind_error,
            dtype='int16').values()
        vprint("    All time-series data loaded\n")

        # --------------------------------------------------------------------
        # find segments and transitions

        self.transitions = dict((v, []) for v in self.grdstates)

        for sdata, rdata, ndata, okdata in zip(*alldata[:4]):
            ssegs = DataQualityDict()
            rsegs = DataQualityDict()
            nsegs = DataQualityDict()
            oksegs = (okdata == 1).to_dqflag(name='Node OK')
            for v, name in self.grdstates.iteritems():
                # get segments for state
                tag = self.segmenttag % name
                instate = sdata == v
                ssegs[tag] = instate.to_dqflag(name=name)
                transin = (numpy.diff(
                    instate.astype(int)) == 1).nonzero()[0] + 1
                transout = (numpy.diff(
                    instate.astype(int)) == -1).nonzero()[0] + 1
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
            config=config, nds=nds, multiprocess=multiprocess,
            datacache=datacache, segmentcache=segmentcache, **kwargs)

    def write_state_html(self, state):
        """Write the HTML for the given state of this `GuardianTab`
        """
        page = self.scaffold_plots()
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
        page.h2('State transitions')
        page.p("The table below lists for each state (row), the number of "
               "transitions into that state from each other state (columns).")
        page.p("Only those states named in the configuration are shown, but "
               "the 'Total' includes transitions from any and all states. "
               "'-' indicates no transitions from that state")
        page.table(class_='table table-condensed table-hover '
                          'table-responsive transitions')
        page.thead()
        page.tr()
        for th in ['State'] + self.grdstates.values() + ['Total']:
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
            page.a(href='#collapse%d' % bit,
                   **{'data-toggle': 'collapse', 'data-parent': '#accordion'})
            page.h4(name, class_='panel-title')
            page.a.close()
            page.div.close()

            # body
            page.div(id='collapse%d' % bit, class_='panel-collapse collapse')
            page.div(class_='panel-body')

            # print transitions
            page.p('This state was active %d times as follows:'
                   % len(self.transitions[bit]))
            headers = ['GPS time', 'UTC time', 'Local time',
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
            page.add(str(html.data_table(headers, data, table='guardian data')))

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
            page.div.close()
            page.div.close()
            page.div.close()
        page.div.close()
        page.div.close()
        page.div.close()

        return super(DataTab, self).write_state_html(state, plots=False,
                                                     pre=page)
register_tab(GuardianTab)


class GuardianStatePlot(get_plot('segments')):
    type = 'guardian'
    defaults = get_plot('segments').defaults.copy()
    defaults.update({
        'color': None,
        'edgecolor': 'black',
        'linewidth': 0.5,
        'requestcolor': (0., .4, 1.),
        'nominalcolor': (1.0, 0.7, 0.0),
        'legend_loc': 'upper left',
        'legend_bbox_to_anchor': (1.01, 1),
        'legend_borderaxespad': 0.,
        'legend_fontsize': 12,
        'legend_title': 'Node state',
        'ytick.labelsize': 10,
    })

    def __init__(self, *args, **kwargs):
        super(GuardianStatePlot, self).__init__(*args, **kwargs)
        self.preview_labels = True

    @property
    def node(self):
        return self.flags[0].split(' ', 1)[0][3:]

    @property
    def ifo(self):
        return list(self.ifos)[0]

    def draw(self):
        (plot, axes) = self.init_plot(plot=SegmentPlot)
        ax = axes[0]

        # get labels
        flags = map(lambda f: str(f).replace('_', r'\_'), self.flags)
        labels = self.pargs.pop('labels', self.pargs.pop('label', flags))
        ax.set_insetlabels(self.pargs.pop('insetlabels', True))
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        labels = map(lambda s: re_quote.sub('', str(s).strip('\n ')), labels)

        # parse plotting arguments
        legendargs = self.parse_legend_kwargs()
        plotargs = self.parse_plot_kwargs()[0]
        plotargs.pop('label')
        activecolor = plotargs.pop('facecolor')
        nominalcolor = self.pargs.pop('nominalcolor')
        requestcolor = self.pargs.pop('requestcolor')
        actargs = plotargs.copy()
        plotargs.update({
            'facecolor': nominalcolor,
            'edgecolor': 'none',
            'known': {'alpha': 0.1, 'facecolor': 'lightgray'},
        })
        reqargs = plotargs.copy()
        reqargs.update({
            'facecolor': requestcolor,
            'known': None,
            })
        actargs.update({
            'facecolor': activecolor,
            'known': None,
        })

        if self.state and not self.all_data:
            valid = self.state.active
        else:
            valid = SegmentList([self.span])

        # plot segments
        for y, (flag, label) in enumerate(zip(self.flags, labels)[::-1]):
            inreq = str(flag) + REQUESTSTUB
            nominal = str(flag) + NOMINALSTUB
            segs = get_segments([flag, inreq, nominal], validity=valid,
                                query=False)
            # format label
            if self.fileformat != 'svg':
                try:
                    idx, label = re_guardian_index.match(label).groups()
                except AttributeError:
                    pass
                else:
                    x = float(self.span[0]) - (float(abs(self.span)) * 0.005)
                    ax.text(x, y, '[%s]' % idx, ha='right', va='center',
                            fontsize=12)
            # plot segments
            ax.plot(segs[nominal], label=label, y=y, height=1., **plotargs)
            ax.plot(segs[inreq], label=label, y=y, collection='ignore',
                    **reqargs)
            ax.plot(segs[flag], label=label, y=y, collection='ignore',
                    height=.6, **actargs)

        # make custom legend
        seg = Segment(self.start - 10, self.start - 9)
        v = plotargs.pop('known', None)
        if v:
            v.pop('collection', None)
            v = ax.build_segment(seg, 0, **v)
        n = ax.build_segment(seg, 0, facecolor=nominalcolor,
                             edgecolor=plotargs['edgecolor'])
        a = ax.build_segment(seg, 0, facecolor=requestcolor,
                             edgecolor=plotargs['edgecolor'])
        b = ax.build_segment(seg, 0, facecolor=activecolor,
                             edgecolor=actargs['edgecolor'])
        if v:
            ax.legend([v, n, a, b], ['Alive', 'Nominal', 'Request', 'Active'],
                      **legendargs)
        else:
            ax.legend([n, a, b], ['Nominal', 'Request', 'Active'], **legendargs)

        # customise plot
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if 'ylim' not in self.pargs:
            ax.set_ylim(-0.5, len(self.flags) - 0.5)

        # fake colorbar
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)

        # add node MODE along the bottom
        sax = None
        legentry = OrderedDict()
        grdmode = get_timeseries(
           '%s:GRD-%s_MODE' % (self.ifo, self.node),
           valid, query=False).join(gap='pad', pad=-1)
        try:
           grdop = get_timeseries(
              '%s:GRD-%s_OP' % (self.ifo, self.node),
              valid, query=False).join(gap='pad', pad=-1)
        except KeyError:
            modes = [(grdmode, (None, 'PAUSE', 'EXEC', 'MANAGED'))]
            colors = ('yellow', (0., .4, 1.), (.5, .0, .75))
        else:
            modes = [(grdmode, ('EXEC', 'MANAGAED', 'MANUAL')),
                     (grdop, (None, 'PAUSE', None))]
            colors = ((0., .4, 1.), (.5, .0, .75), 'hotpink', 'yellow')
        cidx = 0
        for i, (data, mstate) in enumerate(modes):
            for j, m in enumerate(mstate):
                if m is None:
                    continue
                try:
                    x = (data == j).to_dqflag()
                except KeyError:
                    pass
                else:
                    if sax is None:
                        sax = plot.add_state_segments(
                            x.active,
                            plotargs={'y': 0, 'facecolor': colors[cidx],
                                      'edgecolor': 'none',
                                      'collection': 'ignore'})
                    else:
                        sax.plot_segmentlist(
                            x.active, facecolor=colors[cidx], y=0,
                            edgecolor='none', collection='ignore')
                legentry[m.title()] = ax.build_segment(
                    seg, 0, facecolor=colors[cidx], edgecolor='none')
                cidx += 1

        # add OK segments along the bottom
        try:
            ok = get_segments(flag.split(' ', 1)[0] + ' OK', validity=valid,
                              query=False)
        except NameError:
            pass
        else:
            # just plot OK
            if sax is not None:
                sax.plot_segmentlist(ok.active, facecolor=activecolor, y=0,
                                     edgecolor=actargs['edgecolor'],
                                     collection='ignore', label='Node `OK\'')
                # add gap in legend
                legentry['$--$'] = ax.build_segment(
                    seg, 0, facecolor='w', fill=False, edgecolor='none',
                    linewidth=0)
                # add OK to legend
                legentry['Node `OK\''] = ax.build_segment(
                    seg, 0, facecolor=activecolor,
                    edgecolor=actargs['edgecolor'])
            else:
                sax = plot.add_state_segments(ok, ax, plotargs={
                    'label': 'Node `OK\'', 'y': 0, 'facecolor': activecolor,
                    'known': {'facecolor': 'red', 'edgecolor': 'black'}})

        # make custom legend
        if sax is not None:
            if legentry:
                sax.legend(legentry.values(), legentry.keys(), loc='lower left',
                           bbox_to_anchor=(1.01, 0.), borderaxespad=0,
                           fontsize=12, title='Node mode')
            sax.tick_params(axis='y', which='major', labelsize=12)
            sax.set_epoch(float(self.pargs.get('epoch', self.start)))
            ax.set_epoch(sax.get_epoch())

        return self.finalize()

register_plot(GuardianStatePlot)
