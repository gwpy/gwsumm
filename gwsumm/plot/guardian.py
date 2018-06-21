# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2018)
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

"""Plots of Guardian data
"""

from collections import OrderedDict

from six import string_types

from gwpy.segments import (Segment, SegmentList)
from gwpy.plot.segments import SegmentRectangle

from ..data import get_timeseries
from ..segments import get_segments
from ..utils import re_quote
from .registry import (get_plot, register_plot)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


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
        from ..tabs.guardian import (REQUESTSTUB, NOMINALSTUB,
                                     re_guardian_index)

        plot = self.init_plot()
        ax = plot.gca()

        # get labels
        flags = map(lambda f: str(f).replace('_', r'\_'), self.flags)
        labels = self.pargs.pop('labels', self.pargs.pop('label', flags))
        ax.set_insetlabels(self.pargs.pop('insetlabels', True))
        if isinstance(labels, string_types):
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
        n = SegmentRectangle(seg, 0, facecolor=nominalcolor,
                             edgecolor=plotargs['edgecolor'])
        a = SegmentRectangle(seg, 0, facecolor=requestcolor,
                             edgecolor=plotargs['edgecolor'])
        b = SegmentRectangle(seg, 0, facecolor=activecolor,
                             edgecolor=actargs['edgecolor'])
        if v:
            v.pop('collection', None)
            v = SegmentRectangle(seg, 0, **v)
            ax.legend([v, n, a, b], ['Alive', 'Nominal', 'Request', 'Active'],
                      **legendargs)
        else:
            ax.legend([n, a, b], ['Nominal', 'Request', 'Active'],
                      **legendargs)

        # customise plot
        for key, val in self.pargs.items():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if 'ylim' not in self.pargs:
            ax.set_ylim(-0.5, len(self.flags) - 0.5)

        # add node MODE along the bottom
        sax = None
        legentry = OrderedDict()
        grdmode = get_timeseries(
            '%s:GRD-%s_MODE' % (self.ifo, self.node),
            valid, query=False).join(gap='pad', pad=-1)
        grdop = get_timeseries(
            '%s:GRD-%s_OP' % (self.ifo, self.node),
            valid, query=False).join(gap='pad', pad=-1)
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
                    seg_kw = {'y': 0, 'facecolor': colors[cidx],
                              'edgecolor': 'none', 'collection': 'ignore'}
                    if sax is None:
                        sax = plot.add_segments_bar(x.active, **seg_kw)
                    else:
                        sax.plot_segmentlist(x.active, **seg_kw)
                legentry[m.title()] = SegmentRectangle(
                    seg, 0, facecolor=colors[cidx], edgecolor='none')
                cidx += 1

        # add OK segments along the bottom
        try:
            ok = get_segments(flag.split(' ', 1)[0] + ' OK', validity=valid,
                              query=False)
        except NameError:
            pass
        else:
            seg_kw = {'y': 0, 'facecolor': activecolor, 'label': 'Node `OK\''}
            # just plot OK
            if sax is not None:
                sax.plot_segmentlist(ok.active, edgecolor=actargs['edgecolor'],
                                     collection='ignore', **seg_kw)
                # add gap in legend
                legentry['$--$'] = SegmentRectangle(
                    seg, 0, facecolor='w', fill=False, edgecolor='none',
                    linewidth=0)
                # add OK to legend
                legentry['Node `OK\''] = SegmentRectangle(
                    seg, 0, facecolor=activecolor,
                    edgecolor=actargs['edgecolor'])
            else:
                sax = plot.add_segments_bar(ok, ax, known='red', **seg_kw)

        # make custom legend
        if sax is not None:
            if legentry:
                sax.legend(legentry.values(), legentry.keys(),
                           loc='lower left', bbox_to_anchor=(1.01, 0.),
                           borderaxespad=0, fontsize=12, title='Node mode')
            sax.tick_params(axis='y', which='major', labelsize=12)
            sax.set_epoch(float(self.pargs.get('epoch', self.start)))
            ax.set_epoch(sax.get_epoch())

        return self.finalize()

register_plot(GuardianStatePlot)
