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

"""Definition of the `AccountingTab`
"""

import re
from collections import OrderedDict
from configparser import NoOptionError

from glue.lal import Cache

from gwpy.segments import (DataQualityDict, SegmentList)

from gwdetchar.io import html

from ..config import GWSummConfigParser
from .registry import (get_tab, register_tab)
from .. import globalv
from ..data import get_timeseries
from ..segments import get_segments
from ..plot.registry import get_plot
from ..utils import vprint

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['AccountingTab']

ParentTab = get_tab('default')


class AccountingTab(ParentTab):
    """Summarise the data recorded by the operating mode channels
    """
    type = 'accounting'

    @classmethod
    def from_ini(cls, config, section, plotdir='plots', **kwargs):
        new = super(ParentTab, cls).from_ini(config, section, **kwargs)

        # add information
        new.plots = []
        new.channel = config.get(section, 'channel')
        new.ifo = config.get(section, 'IFO')
        new.modes = OrderedDict((int(idx), name) for (idx, name) in
                                config.nditems(section) if idx.isdigit())

        # -----------
        # build plots

        new.set_layout([2, 2])
        new.segmenttag = '%s:%%s' % (new.channel)
        tag = new.channel.split(':', 1)[-1].replace('-', '_')

        groups = type(new.modes)(
            (idx, name.strip('*')) for (idx, name) in new.modes.items()
            if int(idx) % 10 == 0)
        pstates = OrderedDict(
            (idx, name) for (idx, name) in new.modes.items() if
            not name.startswith('*'))

        # parse colors
        colors = dict(
            (int(key[6:]), eval(color)) for (key, color) in
            config.nditems(section) if re.match(r'color-\d+', key))
        for flag in pstates:
            group = int(flag // 10 * 10)
            if flag not in colors:
                colors[flag] = colors.get(group, None)

        # plot segments
        for flags, ptag in zip([groups, pstates], ['overview', 'details']):
            segcolors = [colors.get(flag, None) for flag in flags]
            if not any(segcolors):
                segcolors = None
            for state in new.states:
                new.plots.append(get_plot('segments')(
                    [new.segmenttag % idx for idx in flags],
                    new.span[0], new.span[1], labels=list(flags.values()),
                    state=state, outdir=plotdir,
                    pid='%s_SEGMENTS_%s' % (tag, ptag.upper()),
                    active=segcolors,
                    known={'alpha': 0.1, 'facecolor': 'lightgray'},
                    title='%s operating mode %s' % (new.ifo, ptag)))

        # plot pie charts
        try:
            explode = map(float, config.get(section, 'pie-explode').split(','))
        except NoOptionError:
            explode = None
        ptag = 'overview'
        piecolors = [colors.get(flag, None) for flag in groups]
        if not any(piecolors):
            piecolors = None

        for state in new.states:
            labels = list(groups.values())
            new.plots.append(get_plot('segment-pie')(
                [new.segmenttag % idx for idx in groups],
                new.span[0], new.span[1], state=state, labels=labels,
                outdir=plotdir, pid='%s_PIE_%s' % (tag, ptag.upper()),
                colors=piecolors, explode=explode,
                title='%s operating mode %s' % (new.ifo, ptag)))

            shortlabels = ['\n'.join('{0:.9}.'.format(w) if len(w) > 9 else w
                                     for w in l.split()) for
                           l in labels]
            new.plots.append(get_plot('duty')(
                [new.segmenttag % idx for idx in groups],
                new.span[0], new.span[1], state=state, labels=shortlabels,
                outdir=plotdir, pid='%s_SEGMENT_BAR_%s' % (tag, ptag.upper()),
                colors=piecolors, stacked=True, ylim=[0, 100],
                ylabel=r'Percentage [\%] of available time',
                legend_loc='upper left', legend_bbox_to_anchor=(1., 1),
                legend_fontsize=12, legend_borderaxespad=0,
                legend_frameon=False,
                legend_handlelength=1.0, legend_handletextpad=.5,
                title='%s operating mode %s' % (new.ifo, ptag)))

        return new

    def process(self, nds=None, nproc=1,
                config=GWSummConfigParser(), datacache=None,
                datafind_error='raise', **kwargs):
        """Process time accounting data
        """
        for p in self.plots:
            if p.outputfile in globalv.WRITTEN_PLOTS:
                p.new = False

        # get archived GPS time
        tag = self.segmenttag % list(self.modes)[0]
        try:
            lastgps = globalv.SEGMENTS[tag].known[-1][-1]
        except (IndexError, KeyError):
            lastgps = self.span[0]

        # get data
        new = SegmentList([type(self.span)(lastgps, self.span[1])])
        data = get_timeseries(self.channel, new,
                              config=config, nds=nds, dtype='int16',
                              datafind_error=datafind_error,
                              nproc=nproc, cache=datacache)
        vprint("    All time-series data loaded\n")

        # find segments
        for ts in data:
            modesegments = DataQualityDict()
            for idx, name in self.modes.items():
                # get segments for state
                tag = self.segmenttag % idx
                instate = ts == idx * ts.unit
                modesegments[tag] = instate.to_dqflag(name=name.strip('*'))
                # append segments for group
                group = int(idx // 10. * 10)
                gtag = self.segmenttag % group
                try:
                    modesegments[gtag] += modesegments[tag]
                except KeyError:
                    modesegments[gtag] = modesegments[tag]
            globalv.SEGMENTS += modesegments

        kwargs['segdb_error'] = 'ignore'
        super(AccountingTab, self).process(
            config=config, nds=nds, nproc=nproc,
            segmentcache=Cache(), datacache=datacache,
            datafind_error=datafind_error, **kwargs)

    def write_state_html(self, state):
        """Write the HTML for the given state of this `GuardianTab`
        """
        page = self.scaffold_plots(state=state)

        page.div(class_='row')
        page.div(class_='col-md-12')

        # get segment data
        groups = type(self.modes)((idx, name) for idx, name in
                                  self.modes.items() if idx % 10 == 0)
        modes = type(self.modes)(
            (idx, name) for idx, name in self.modes.items() if
            not name.startswith('*'))
        headers = ['Index', 'Name', 'Active seconds',
                   'Hours', '%']
        prefix = self.ifo.lower()
        for flags, title, id in zip(
                [groups, modes],
                ['Top-level mode information', 'Detailed mode information'],
                ['%s-top-level-mode' % prefix, '%s-detailed-mode' % prefix]):
            page.h1(title)
            data = []
            pc = float(abs(state.active) / 100.)
            tots = 0
            toth = 0
            totp = 0
            for idx, name in flags.items():
                flag = get_segments(self.segmenttag % idx,
                                    state.active, query=False).copy()
                actives = abs(flag.active)
                activeh = actives / 3600.
                try:
                    activep = actives / pc
                except ZeroDivisionError:
                    activep = 0
                tots += actives
                toth += activeh
                totp += activep
                data.append([str(idx), name.strip('*'), '%.1f' % actives,
                             '%.1f' % activeh, '%.1f' % activep])
            data.append(map(
                lambda x: '<strong>%s</strong>' % x,
                ['', 'Total:', '%.1f' % tots, '%.1f' % toth, '%.1f' % totp]))
            page.add(str(html.table(
                headers, data, id=id,
                caption="%s observatory mode statistics as recorded in "
                        "<samp>%s</samp>" % (title.split()[0], self.channel))))
        return super(ParentTab, self).write_state_html(state, plots=False,
                                                       pre=page)


register_tab(AccountingTab)
