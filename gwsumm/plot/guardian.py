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
from gwpy.plot.colors import tint
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
        'linewidth': 0.5,
        'requestcolor': '#0066ff',
        'nominalcolor': '#ffb200',
        'legend_loc': 'upper left',
        'legend_bbox_to_anchor': (1., 1.),
        'legend_borderaxespad': 0.,
        'legend_frameon': False,
        'legend_fontsize': 12,
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
        height = plotargs.pop('height', .8)
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
            'height': height,
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
        for y, (flag, label) in enumerate(list(zip(self.flags, labels))[::-1]):
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
                             edgecolor=nominalcolor)
        a = SegmentRectangle(seg, 0, facecolor=requestcolor,
                             edgecolor=requestcolor)
        b = SegmentRectangle(seg, 0, facecolor=activecolor,
                             edgecolor=actargs['edgecolor'])
        handles = [n, a, b]
        labels = ['Nominal', 'Request', 'Active']
        if v:
            v.pop('collection', None)
            v['edgecolor'] = v.get('edgecolor') or tint(v['facecolor'], .5)
            handles.insert(0, SegmentRectangle(seg, 0, **v))
            labels.insert(0, 'Alive')
        ax.legend(handles, labels, title='Node state', **legendargs)

        # customise plot
        for key, val in self.pargs.items():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if 'ylim' not in self.pargs:
            ax.set_ylim(-0.5, len(self.flags) - 0.5)

        epoch = ax.get_epoch()

        # add node MODE along the bottom
        sax = None
        seg_kw = {'y': 0, 'edgecolor': 'none', 'label': 'Mode'}
        legentry = OrderedDict()
        colors = iter(('#0066ff', '#8000bf', 'hotpink', 'yellow'))
        for ctag, mstate in [
            # bit listing for channels (Nones are ignored)
            ('MODE', ('EXEC', 'MANAGED', 'MANUAL')),
            ('OP', (None, 'PAUSE', None)),
        ]:
            # get data
            data = get_timeseries(
                '{0.ifo}:GRD-{0.node}_{1}'.format(self, ctag), valid,
                query=False).join(gap='pad', pad=-1)
            for i, m in filter(lambda x: x[1] is not None, enumerate(mstate)):
                x = (data == i).to_dqflag()
                fc = next(colors)
                if sax is None:
                    sax = plot.add_segments_bar(
                        x.active, facecolor=fc, **seg_kw)
                else:
                    sax.plot_segmentlist(
                        x.active, facecolor=fc, **seg_kw)
                legentry[m.title()] = SegmentRectangle(seg, 0, facecolor=fc,
                                                       edgecolor=fc)
                seg_kw.update({'collection': 'ignore', 'label': None})

        # add OK segments along the bottom
        ok = get_segments(flag.split(' ', 1)[0] + ' OK', validity=valid,
                          query=False)
        seg_kw.pop('edgecolor', None)
        sax.plot_segmentlist(ok.active, facecolor=activecolor,
                             edgecolor=actargs['edgecolor'],
                             **seg_kw)
        legentry['`OK\''] = SegmentRectangle(seg, 0, facecolor=activecolor,
                                             edgecolor=actargs['edgecolor'])

        # make custom legend
        legendargs.update({
            'bbox_to_anchor': (1., 0.),
            'loc': 'lower left',
        })
        sax.legend(list(legentry.values()), list(legentry),
                   title='Node mode', **legendargs)
        sax.tick_params(axis='y', which='major', labelsize=12)
        for ax_ in (ax, sax):
            ax_.set_epoch(epoch)

        return self.finalize()


register_plot(GuardianStatePlot)
