# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013)
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

"""Definitions for the standard plots
"""

from __future__ import division

import bisect
from itertools import (cycle, combinations)
from numbers import Number
from collections import OrderedDict
from configparser import NoOptionError

from six import string_types

import numpy

from dateutil.relativedelta import relativedelta

from matplotlib import rcParams
from matplotlib.artist import setp
from matplotlib.cbook import iterable
from matplotlib.colors import (rgb2hex, is_color_like)
from matplotlib.patches import Rectangle

from glue import iterutils

from gwpy.plot.colors import (GW_OBSERVATORY_COLORS, tint)
from gwpy.plot.segments import SegmentRectangle
from gwpy.segments import (Segment, SegmentList, DataQualityFlag)
from gwpy.time import (from_gps, to_gps)

from .. import globalv
from ..mode import (Mode, get_mode)
from ..utils import (re_quote, get_odc_bitmask, re_flagdiv, safe_eval)
from ..channels import (get_channel, re_channel)
from ..data import get_timeseries
from ..segments import (get_segments, format_padding)
from ..state import ALLSTATE
from .core import (BarPlot, PiePlot, format_label)
from .registry import (get_plot, register_plot)
from .mixins import SegmentLabelSvgMixin
from .utils import (hash, usetex_tex)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

TimeSeriesDataPlot = get_plot('timeseries')
GREEN = '#33cc33'


def tint_hex(*args, **kwargs):
    return rgb2hex(tint(*args, **kwargs))


def common_limits(datasets, default_min=0, default_max=0):
    """Find the global maxima and minima of a list of datasets.

    Parameters
    ----------
    datasets : `iterable`
        list (or any other iterable) of data arrays to analyse.

    default_min : `float`, optional
        fall-back minimum value if datasets are all empty.

    default_max : `float`, optional
        fall-back maximum value if datasets are all empty.

    Returns
    -------
    (min, max) : `float`
        2-tuple of common minimum and maximum over all datasets.
    """
    if isinstance(datasets, numpy.ndarray) or not iterable(datasets[0]):
        datasets = [datasets]
    max_stat = max(list(iterutils.flatten(datasets)) + [-numpy.inf])
    min_stat = min(list(iterutils.flatten(datasets)) + [numpy.inf])
    if numpy.isinf(-max_stat):
        max_stat = default_max
    if numpy.isinf(min_stat):
        min_stat = default_min
    return min_stat, max_stat


class SegmentDataPlot(SegmentLabelSvgMixin, TimeSeriesDataPlot):
    """Segment plot of one or more `DataQualityFlags <DataQualityFlag>`.
    """
    type = 'segments'
    data = 'segments'
    defaults = TimeSeriesDataPlot.defaults.copy()
    defaults.update({
        'mask': None,
        'color': None,
        'on-is-bad': False,
        'insetlabels': 'inset',
        'legend-bbox_to_anchor': (1., 1.),
        'legend-loc': 'upper left',
        'legend-borderaxespad': 0,
        'legend-fontsize': 12,
        'legend-frameon': False,
        'legend-handletextpad': .5,
    })
    DRAW_PARAMS = TimeSeriesDataPlot.DRAW_PARAMS + [
        'known', 'height', 'y', 'facecolor', 'edgecolor',
    ]

    def __init__(self, flags, start, end, state=None, outdir='.', **kwargs):
        padding = kwargs.pop('padding', None)
        super(SegmentDataPlot, self).__init__([], start, end, state=state,
                                              outdir=outdir, **kwargs)
        self._allflags = []
        self.flags = flags
        self.preview_labels = False
        self.padding = padding

    def get_channel_groups(self, *args, **kwargs):
        return [(f, [f]) for f in self.flags]

    @property
    def flags(self):
        return [f.name for f in self._flags]

    @flags.setter
    def flags(self, flist):
        if isinstance(flist, string_types):
            flist = [f.strip('\n ') for f in flist.split(',')]
        self._flags = []
        for f in flist:
            self.add_flag(f)

    def add_flag(self, f):
        # append flag to main list
        if isinstance(f, DataQualityFlag):
            self._flags.append(f)
        else:
            self._flags.append(DataQualityFlag(f))
        # append raw flags to 'allflags' property
        flags = re_flagdiv.split(str(f))[::2]
        for f in flags:
            if not f:
                continue
            self._allflags.append(DataQualityFlag(f))

    @property
    def allflags(self):
        return [f.name for f in self._allflags]

    @property
    def padding(self):
        return OrderedDict((f.name, f.padding) for f in self._allflags)

    @padding.setter
    def padding(self, pad):
        for f, p in format_padding(self._allflags, pad).items():
            if isinstance(p, (float, int)):
                f.padding = (p, p)
            else:
                f.padding = p

    @property
    def ifos(self):
        """Interferometer set for this `SegmentDataPlot`
        """
        return set([f.strip('!&-_')[:2] for f in self.allflags])

    @property
    def pid(self):
        """File pid for this `DataPlot`.
        """
        try:
            return self._pid
        except AttributeError:
            self._pid = hash("".join(map(str, self.flags)))
            return self.pid

    @pid.setter
    def pid(self, id_):
        self._pid = str(id_)

    @classmethod
    def from_ini(cls, config, section, start, end, flags=None, state=ALLSTATE,
                 **kwargs):
        # get padding
        try:
            kwargs.setdefault(
                'padding', config.get(section, 'padding'))
        except NoOptionError:
            pass
        if 'padding' in kwargs:
            kwargs['padding'] = list(eval(kwargs['padding']))
        # build figure
        new = super(SegmentDataPlot, cls).from_ini(config, section, start,
                                                   end, state=state, **kwargs)
        # get flags
        if flags is None:
            flags = dict(config.items(section)).pop('flags', [])
        if isinstance(flags, string_types):
            flags = [f.strip('\n ') for f in flags.split(',')]
        new.flags = flags
        return new

    def init_plot(self, projection='segments', **kwargs):
        return super(SegmentDataPlot, self).init_plot(
            projection=projection, **kwargs)

    def get_segment_color(self):
        """Parse the configured ``pargs`` and determine the colors for
        active and valid segments.
        """
        active = safe_eval(
            self.pargs.pop('active', self.pargs.pop('facecolor', None)))
        known = safe_eval(self.pargs.pop('known', 0))
        # neither known nor active defined
        if active is None and known == 0:
            self.pargs['facecolor'] = '#33cc33'
            self.pargs['known'] = 'red'
        # only active is defined
        elif known == 0:
            if isinstance(active, dict):
                self.pargs.update(active)
                active = active.get('facecolor')
            else:
                self.pargs['facecolor'] = active
            if (isinstance(active, string_types) and
                    active.lower() in ('red', '#ff0000')):
                self.pargs['known'] = 'dodgerblue'
            else:
                self.pargs['known'] = 'red'
        # only known is defined
        elif active is None:
            self.pargs['known'] = known
            if known in ['#33cc33', 'green', 'g']:
                self.pargs['facecolor'] = 'dodgerblue'
            else:
                self.pargs['facecolor'] = '#33cc33'
        # both are given
        else:
            if isinstance(active, dict):
                self.pargs.update(active)
            else:
                self.pargs['facecolor'] = active
            self.pargs['known'] = known

        # format defaults
        self.pargs.setdefault('height', .8)
        if isinstance(self.pargs['known'], string_types):
            self.pargs['known'] = {'facecolor': self.pargs['known']}

        for dtup in (self.pargs, self.pargs['known']):
            # allow user to give tuple of dicts for 'known'
            if not isinstance(dtup, (list, tuple)):
                dtup = [dtup]
            [SegmentDataPlot._set_default_edgecolor(d) for d in dtup if
             d is not None]

        # set default height for known
        if (isinstance(self.pargs['known'], dict) and
                isinstance(self.pargs['height'], Number)):
            self.pargs['known'].setdefault('height', self.pargs['height'] * .5)

        return self.pargs

    @staticmethod
    def _set_default_edgecolor(pargs):
        """Set the default edgecolor based on the given facecolor
        """
        fc = pargs.get('facecolor')
        ec = pargs.get('edgecolor')
        # if list of colors, map list of edgecolors
        if (not ec and isinstance(fc, (list, tuple)) and
                not is_color_like(fc)):
            pargs['edgecolor'] = [tint_hex(x, factor=.5) for x in fc]
        # otherwise map single color
        elif fc and not ec:
            pargs['edgecolor'] = tint_hex(fc, factor=.5)

    def parse_plot_kwargs(self, *args, **kwargs):
        self.get_segment_color()
        return super(SegmentDataPlot, self).parse_plot_kwargs(*args, **kwargs)

    def draw(self):
        # get labelsize
        labelsize = self.pargs.pop('labelsize', 12)
        if self.pargs.get('insetlabels', True) is False:
            rcParams['ytick.labelsize'] = labelsize

        # create figure
        plot = self.init_plot()
        ax = plot.gca()

        # extract plotting arguments
        legendargs = self.parse_legend_kwargs()
        plotargs = self.parse_plot_kwargs()
        legcolors = plotargs[0].copy()

        # plot segments
        for i, (flag, pargs) in enumerate(
                list(zip(self.flags, plotargs))[::-1]):
            label = re_quote.sub('', pargs.pop('label', str(flag)))
            if (self.fileformat == 'svg' and not str(flag) in label and
                    ax.get_insetlabels()):
                label = '%s [%s]' % (label, str(flag))
            elif self.fileformat == 'svg' and not str(flag) in label:
                label = '[%s] %s' % (label, str(flag))
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            segs = get_segments(flag, validity=valid, query=False,
                                padding=self.padding).coalesce()
            if self.pargs.get('on-is-bad', False):
                segs = ~segs
            pargs.setdefault('known', None)
            pargs.setdefault('y', i)
            ax.plot(segs, label=label, **pargs)

        # make custom legend
        if legcolors.get('known', None):
            self.add_legend(ax, legcolors, **legendargs)

        # customise plot
        if ax.get_autoscaley_on():
            self.pargs['ylim'] = (-.5, len(self.flags) - 0.5)
        self.apply_parameters(ax, **self.pargs)

        # finalise
        self.add_state_segments(ax)
        self.add_future_shade()

        return self.finalize()

    def add_legend(self, ax, colors, **kwargs):
        aface = colors['facecolor']
        aedge = colors['edgecolor']
        kface = colors['known']
        kedge = None
        if isinstance(kface, dict):
            kedge = kface['edgecolor']
            kface = kface['facecolor']

        # draw dummy segments for known, active and edges, and create legend
        seg = Segment(0, 1)
        k = SegmentRectangle(seg, 0, facecolor=kface, edgecolor=kedge)
        a = SegmentRectangle(seg, 0, facecolor=aface, edgecolor=aedge)
        if aedge not in (None, 'none', aface):
            t = SegmentRectangle(seg, 0, facecolor=aedge, edgecolor=aedge)
            return ax.legend([k, a, t], ['Known', 'Active', 'Transition'],
                             **kwargs)
        return ax.legend([k, a], ['Known', 'Active'], **kwargs)


register_plot(SegmentDataPlot)


class StateVectorDataPlot(TimeSeriesDataPlot):
    """DataPlot of some `StateVector` data.

    While technically a sub-class of the `TimeSeriesDataPlot`, for
    data access and processing reasons, the output shadows that of the
    `SegmentDataPlot` more closely.
    """
    type = 'statevector'
    data = 'statevector'
    defaults = SegmentDataPlot.defaults.copy()
    DRAW_PARAMS = list(SegmentDataPlot.DRAW_PARAMS)

    # copy from SegmentDataPlot
    flag = property(fget=SegmentDataPlot.flags.__get__,
                    fset=SegmentDataPlot.flags.__set__,
                    fdel=SegmentDataPlot.flags.__delete__,
                    doc="""List of flags generated for this
                        `StateVectorDataPlot`.""")
    get_segment_color = SegmentDataPlot.__dict__['get_segment_color']

    def __init__(self, *args, **kwargs):
        super(StateVectorDataPlot, self).__init__(*args, **kwargs)
        self.flags = []

    @property
    def pid(self):
        try:
            return self._pid
        except AttributeError:
            basis = "".join(map(str, self.channels))
            if self.pargs.get('bits', None):
                basis += str(self.pargs['bits'])
            self._pid = hash(basis)
            return self.pid

    def _parse_labels(self, defaults=[]):
        """Pop the labels for plotting from the `pargs` for this Plot

        This method overrides from the `TimeSeriesDataPlot` in order
        to set the bit names from the various channels as the defaults
        in stead of the channel names
        """
        chans = list(zip(*self.get_channel_groups()))[0]
        labels = list(self.pargs.pop('labels', defaults))
        if isinstance(labels, string_types):
            labels = labels.split(',')
        for i, l in enumerate(labels):
            if isinstance(l, (list, tuple)):
                labels[i] = list(labels[i])
                for j, l2 in enumerate(l):
                    labels[i][j] = format_label(str(l2))
            elif isinstance(l, str):
                labels[i] = format_label(str(l))
        while len(labels) < len(chans):
            labels.append(None)
        return labels

    def parse_plot_kwargs(self, *args, **kwargs):
        self.get_segment_color()
        return super(StateVectorDataPlot, self).parse_plot_kwargs(
            *args, **kwargs)

    def init_plot(self, *args, **kwargs):
        kwargs.setdefault('projection', 'segments')
        return super(StateVectorDataPlot, self).init_plot(*args, **kwargs)

    def draw(self):
        # make font size smaller
        labelsize = self.rcParams.get('ytick.labelsize', 12)
        if self.pargs.get('insetlabels', True) is False:
            rcParams['ytick.labelsize'] = labelsize

        plot = self.init_plot()
        ax = plot.gca()

        # get bit setting
        bits = self.pargs.pop('bits', None)
        if bits and len(self.channels) > 1:
            raise ValueError("Specifying 'bits' doesn't work for a "
                             "state-vector plot including multiple channels")

        # extract plotting arguments
        extraargs = self.parse_plot_kwargs()

        # plot segments
        nflags = 0
        for channel, pargs in zip(self.channels[::-1], extraargs[::-1]):
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            channel = get_channel(channel)
            if bits:
                bits_ = [x if i in bits else None for
                         (i, x) in enumerate(channel.bits)]
            else:
                try:
                    bits_ = channel.bits
                except AttributeError:
                    m = list(re_channel.findall(str(channel)))
                    if len(m) == 1 and hasattr(get_channel(m[0]), 'bits'):
                        bits_ = get_channel(m[0]).bits
                    else:
                        raise
            data = get_timeseries(str(channel), valid, query=False,
                                  statevector=True)
            flags = None
            for stateseries in data:
                if not stateseries.size:
                    stateseries.epoch = self.start
                    stateseries.dx = 0
                    if channel.sample_rate is not None:
                        stateseries.sample_rate = channel.sample_rate
                stateseries.bits = bits_
                if 'int' not in str(stateseries.dtype):
                    stateseries = stateseries.astype('uint32')
                newflags = list(stateseries.to_dqflags().values())
                if self.pargs.get('on-is-bad', False):
                    for i, flag in enumerate(newflags):
                        newflags[i] = ~newflags[i]
                if flags is None:
                    flags = newflags
                else:
                    for i, flag in enumerate(newflags):
                        flags[i] += flag
            if flags is None:
                flags = [DataQualityFlag(b) for b in channel.bits if
                         b not in [None, '']]
            nflags += len([m for m in bits_ if m is not None])
            labels = pargs.pop('label', [None]*len(flags))
            if isinstance(labels, str):
                labels = [labels]
            while len(labels) < len(flags):
                labels.append(None)
            for flag, label in list(zip(flags, labels))[::-1]:
                kwargs = pargs.copy()
                if label is not None:
                    kwargs['label'] = label
                ax.plot(flag, **kwargs)

        # customise plot
        if 'ylim' not in self.pargs:
            self.pargs['ylim'] = (-.5, nflags-.5)
        self.apply_parameters(ax, **self.pargs)

        self.add_state_segments(ax)
        self.add_future_shade()

        return self.finalize()


register_plot(StateVectorDataPlot)


class DutyDataPlot(SegmentDataPlot):
    """`DataPlot` of the duty-factor for a `SegmentList`
    """
    type = 'duty'
    data = 'segments'
    defaults = TimeSeriesDataPlot.defaults.copy()
    defaults.update({
        'alpha': 0.8,
        'sep': False,
        'side_by_side': False,
        'normalized': None,
        'cumulative': False,
        'stacked': False,
        'ylabel': r'Duty factor [\%]',
        'ylim': (0, 100),
    })

    def __init__(self, flags, start, end, state=None, outdir='.',
                 bins=None, **kwargs):
        kwargs.setdefault('fileformat', 'png')
        super(DutyDataPlot, self).__init__(flags, start, end, state=state,
                                           outdir=outdir, **kwargs)
        self.bins = bins

    @property
    def pid(self):
        try:
            return self._pid
        except AttributeError:
            super(DutyDataPlot, self).pid
            if self.pargs.get('cumulative', False):
                self._pid += '_CUMULATIVE'
            return self.pid

    @pid.setter
    def pid(self, p):
        self._pid = p

    def parse_plot_kwargs(self, *args, **kwargs):
        return super(SegmentDataPlot, self).parse_plot_kwargs(*args, **kwargs)

    def get_bins(self):
        """Work out the correct histogram binning for this `DutyDataPlot`
        """
        # if not given anything, work it out from the mode
        if self.bins is None:
            m = get_mode()
            duration = float(abs(self.span))
            # for year mode, use a month
            if m == Mode.year or duration >= 86400 * 300:
                dt = relativedelta(months=1)
            # for more than 8 weeks, use weeks
            elif duration >= 86400 * 7 * 8:
                dt = relativedelta(weeks=1)
            # for week and month mode, use daily
            elif m in [Mode.week, Mode.month] or duration >= 86400 * 7:
                dt = relativedelta(days=1)
            # for day mode, make hourly duty factor
            elif m == Mode.day:
                dt = relativedelta(hours=1)
            # otherwise provide 10 bins
            else:
                dt = relativedelta(seconds=float(abs(self.span))/10.)
        # if given a float, assume this is the bin size
        elif isinstance(self.bins, (float, int)):
            dt = relativedelta(seconds=self.bins)
        # if we don't have a list, we must have worked out dt
        if not isinstance(self.bins, (list, tuple, numpy.ndarray)):
            self.bins = []
            s = from_gps(self.start)
            e = from_gps(self.end)
            while s < e:
                t = int(to_gps(s + dt) - to_gps(s))
                self.bins.append(t)
                s += dt
        self.bins = numpy.asarray(self.bins)
        return self.bins

    def calculate_duty_factor(self, segments, bins=None, cumulative=False,
                              normalized=None):
        if normalized is None and cumulative:
            normalized = False
        elif normalized is None:
            normalized = 'percent'
        if normalized == 'percent':
            normalized = 100.
        else:
            normalized = float(normalized)
        if not bins:
            bins = self.get_bins()
        if isinstance(segments, DataQualityFlag):
            segments = segments.known & segments.active
        duty = numpy.zeros(len(bins))
        mean = numpy.zeros(len(bins))
        for i in range(len(bins)):
            bin = SegmentList([Segment(self.start + float(sum(bins[:i])),
                                       self.start + float(sum(bins[:i+1])))])
            d = float(abs(segments & bin))
            if normalized:
                d *= normalized / bins[i]
            duty[i] = d
            mean[i] = duty[:i+1].mean()
        if cumulative:
            duty = duty.cumsum()
        return duty, mean

    def draw(self, outputfile=None):
        sep = self.pargs.pop('sep', False)
        if sep:
            if self.pargs.get('side_by_side'):
                raise ValueError('DutyDataPlot parameters \'sep\' and '
                                 '\'side_by_side\' should not be used '
                                 'together')
            geometry = (len(self.flags), 1)
        else:
            geometry = (1, 1)

        plot = self.init_plot(geometry=geometry, projection='rectilinear',
                              sharex=True)
        axes = plot.axes

        # extract plotting arguments
        style = self.pargs.pop('style', 'bar')
        stacked = self.pargs.pop('stacked', False)
        sidebyside = self.pargs.pop('side_by_side', False)
        normalized = self.pargs.pop('normalized', True)
        cumulative = self.pargs.pop('cumulative', False)
        if normalized is None and not cumulative:
            normalized = 'percent'
        rollingmean = self.pargs.pop('rolling_mean',
                                     not stacked and not cumulative)
        plotargs = self.parse_plot_kwargs()
        legendargs = self.parse_legend_kwargs()
        if sep:
            legendargs.setdefault('loc', 'upper left')
            legendargs.setdefault('bbox_to_anchor', (1.01, 1))
            legendargs.setdefault('borderaxespad', 0)

        # work out times and plot mean for legend
        self.get_bins()
        times = float(self.start) + numpy.concatenate(
                                 ([0], self.bins[:-1].cumsum()))
        now = bisect.bisect_left(times, globalv.NOW)
        if rollingmean:
            axes[0].plot(times[:1], [-1], 'k--', label='Rolling mean')

        # get bar parameters
        try:
            bottom = axes[0].get_ylim()[0]
        except KeyError:
            bottom = 0
        bottom = numpy.zeros(times.size) + bottom

        # plot segments
        if self.state and not self.all_data:
            valid = self.state.active
        else:
            valid = SegmentList([self.span])
        for i, (ax, flag, pargs, propc) in enumerate(
                zip(cycle(axes), self.flags, plotargs,
                    cycle(rcParams['axes.prop_cycle']))):
            # get segments
            segs = get_segments(flag, validity=valid, query=False,
                                padding=self.padding)
            duty, mean = self.calculate_duty_factor(
                segs, normalized=normalized, cumulative=cumulative)
            # plot duty cycle
            if sep and pargs.get('label') == flag.replace('_', r'\_'):
                pargs.pop('label', None)
            elif 'label' in pargs and normalized == 'percent' and not stacked:
                if legendargs.get('loc', None) in ['upper left', 2]:
                    pargs['label'] = pargs['label'] + '\n[%.1f\\%%]' % mean[-1]
                else:
                    pargs['label'] = pargs['label'] + r' [%.1f\%%]' % mean[-1]
            color = pargs.pop('color', propc['color'])
            # plot in relevant style
            if style == 'line':
                lineargs = pargs.copy()
                lineargs.setdefault('drawstyle', 'steps-post')
                ax.plot(times[:now], duty[:now], color=color, **lineargs)
            elif style not in ['bar', 'fill']:
                raise ValueError("Cannot display %s with style=%r"
                                 % (type(self).__name__, style))
            else:
                # work out positions
                if sidebyside:
                    pad = .1
                    x = 1 - pad * 2
                    w = pargs.pop('width', 1.) * x / len(self.flags)
                    offset = pad + x/len(self.flags) * (i + 1/2.)
                elif stacked:
                    offset = .5
                    w = pargs.pop('width', .9)
                else:
                    offset = .5
                    w = pargs.pop('width', 1.)
                width = w * self.bins[:now]
                if stacked:
                    height = duty
                    pargs.setdefault('edgecolor', color)
                else:
                    height = duty - bottom
                if style == 'fill':
                    width = self.bins[:now]
                    ec = pargs.pop('edgecolor', 'black')
                    pargs['edgecolor'] = 'none'
                    lw = pargs.pop('linewidth', 1)
                    pargs['linewidth'] = 0
                ax.bar((times + self.bins * offset)[:now], height[:now],
                       bottom=bottom[:now], align='center',
                       width=width, color=color, **pargs)
                if style == 'fill':
                    ax.plot(times[:now+1], duty[:now+1],
                            drawstyle='steps-post', color=ec, linewidth=lw)

            # plot mean
            if rollingmean:
                t = [self.start] + list(times + self.bins/2.) + [self.end]
                mean = [mean[0]] + list(mean) + [mean[-1]]
                ax.plot(t, mean, color=sep and 'k' or color, linestyle='--')

            # record duty for stacked chart
            if stacked:
                bottom += height

        # customise plot
        for ax in axes:
            self.apply_parameters(ax, **self.pargs)
            if 'hours' in self.pargs.get('ylabel', ax.get_ylabel()):
                ax.get_yaxis().get_major_locator().set_params(
                    steps=[1, 2, 4, 8])
        if sep:
            # set text
            ylabel = axes[0].yaxis.get_label()
            y = axes[-1].get_position().y0 + (
                axes[0].get_position().y1 - axes[-1].get_position().y0)/2.
            t = plot.text(0.04, y, ylabel.get_text(), rotation=90, ha='center',
                          va='center')
            t.set_fontproperties(ylabel.get_font_properties())
            for i, ax in enumerate(axes):
                ax.set_ylabel('')
                if i:
                    ax.set_title('')
                if i < len(axes) - 1:
                    ax.set_xlabel('')
                    setp(ax.get_xticklabels(), visible=False)

        # add custom legend for mean
        if rollingmean:
            axsize = axes[0].get_position().size
            yoff = 0.01 * axsize[0] / axsize[1]
            lkwargs = legendargs.copy()
            lkwargs.update({
                'loc': 'lower right',
                'bbox_to_anchor': (1.0, 1. + yoff),
                'fontsize': 12,
                'borderaxespad': 0,
            })
            leg = axes[0].legend(['Rolling mean'], **lkwargs)
            if leg.get_frame().get_edgecolor() != 'none':
                leg.get_frame().set_edgecolor(rcParams['grid.color'])
            axes[0].add_artist(leg)
            axes[0].lines[0].set_label('_')

        # add legend
        for ax in axes:
            ax.legend(**legendargs)

        self.add_state_segments(axes[-1])
        self.add_future_shade()

        return self.finalize(outputfile=outputfile)


register_plot(DutyDataPlot)


class ODCDataPlot(SegmentLabelSvgMixin, StateVectorDataPlot):
    """Custom `StateVectorDataPlot` for ODCs with bitmasks
    """
    type = 'odc'
    data = 'odc'
    defaults = StateVectorDataPlot.defaults.copy()
    defaults.update({
        'no_summary_bit': False,
        'in_mask_color': (.0, .4, 1.),
        'masked_off_color': 'red',
        'unmasked_off_color': (1.0, 0.7, 0.0),
        'legend-loc': 'upper left',
        'legend-bbox_to_anchor': (1.01, 1),
        'legend-borderaxespad': 0.,
        'legend-fontsize': 10,
    })

    def __init__(self, *args, **kwargs):
        bitmaskc = kwargs.pop('bitmask_channel', None)
        super(ODCDataPlot, self).__init__(*args, **kwargs)
        if bitmaskc:
            self.bitmask = bitmaskc.split(',')
        else:
            self.bitmask = list(map(get_odc_bitmask, self.channels))

    def get_bitmask_channels(self):
        return type(self.channels)(list(map(get_channel, self.bitmask)))

    @property
    def pid(self):
        try:
            return self._pid
        except AttributeError:
            chans = "".join(map(str, self.channels))
            masks = "".join(map(str, self.get_bitmask_channels()))
            basis = chans + masks
            if self.pargs.get('bits', None):
                basis += str(self.pargs["bits"])
            self._pid = hash(basis)
            return self.pid

    def draw(self):
        # make font size smaller
        labelsize = self.pargs.pop('labelsize', 12)
        rcParams['ytick.labelsize'] = labelsize

        # make figure
        plot = self.init_plot()
        ax = plot.gca()
        ax.grid(False, which='both', axis='y')

        # extract plotting arguments
        nosummary = self.pargs.pop('no_summary_bit', False)
        activecolor = self.pargs.pop('active', GREEN)
        edgecolor = self.pargs.pop('edgecolor', 'black')
        maskoncolor = self.pargs.pop('masked_off_color', 'red')
        maskoffcolor = self.pargs.pop('unmasked_off_color', (1.0, 0.7, 0.0))
        inmaskcolor = self.pargs.pop('in_mask_color', (.0, .4, 1.))
        plotargs = {'facecolor': activecolor,
                    'edgecolor': edgecolor,
                    'height': .8}
        legendargs = self.parse_legend_kwargs()

        # plot segments
        nflags = 0
        for i, (channel, bitmaskchan) in enumerate(
                zip(self.channels, self.get_bitmask_channels())):
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            # read ODC and bitmask vector
            data = get_timeseries(str(channel), valid, query=False,
                                  statevector=True)
            bitmask = get_timeseries(bitmaskchan, valid, query=False,
                                     statevector=True)
            # plot bitmask
            flags = {}
            # plot bits
            for type_, svlist in zip(['bitmask', 'data'], [bitmask, data]):
                flags[type_] = None
                for stateseries in svlist:
                    if not stateseries.size:
                        stateseries.epoch = self.start
                        stateseries.dx = 0
                        if channel.sample_rate is not None:
                            stateseries.sample_rate = channel.sample_rate
                    stateseries.bits = channel.bits
                    if 'int' not in str(stateseries.dtype):
                        stateseries = stateseries.astype('uint32')
                    newflags = stateseries.to_dqflags()
                    if flags[type_] is None:
                        flags[type_] = newflags
                    else:
                        for i, flag in newflags.items():
                            flags[type_][i] += flag
            i = 0
            for i, bit in enumerate(channel.bits):
                if bit is None or bit == '':
                    continue
                try:
                    mask = flags['bitmask'][bit].active
                except TypeError:
                    continue
                segs = flags['data'][bit]
                label = '[%s] %s' % (i, segs.name)
                # plot summary bit
                if segs.name == channel.bits[0] and not nosummary:
                    summargs = plotargs.copy()
                    summargs['height'] *= 3
                    ax.plot(segs, y=-nflags - 1, label=label,
                            known=maskoncolor, **summargs)
                    nflags += 2
                # plot masks and separate masked/not masked
                else:
                    maskon = segs.copy()
                    maskon.known &= mask
                    maskon.active &= mask
                    maskoff = segs.copy()
                    maskoff.known -= mask
                    maskoff.active -= mask
                    # plot mask
                    ax.plot(mask, y=-nflags, facecolor=inmaskcolor,
                            edgecolor='none', height=1., label=None,
                            collection=False, zorder=-1001)
                    # plot mask
                    if maskoff:
                        ax.plot(maskoff, y=-nflags, label=label,
                                known=maskoffcolor, **plotargs)
                        label = None
                    if maskon:
                        ax.plot(maskon, y=-nflags, label=label,
                                known=maskoncolor, **plotargs)

                    label = '[%s] %s' % (i, segs.name)
                nflags += 1

        # make custom legend
        epoch = ax.get_epoch()
        xlim = ax.get_xlim()
        seg = Segment(self.start - 10, self.start - 9)
        m = SegmentRectangle(seg, y=0, facecolor=inmaskcolor, edgecolor='none')
        v = SegmentRectangle(seg, y=0, facecolor=maskoncolor,
                             edgecolor=edgecolor)
        x = SegmentRectangle(seg, y=0, facecolor=maskoffcolor,
                             edgecolor=edgecolor)
        a = SegmentRectangle(seg, y=0, facecolor=activecolor,
                             edgecolor=edgecolor)
        if edgecolor not in [None, 'none']:
            t = SegmentRectangle(seg, y=0, facecolor=edgecolor)
            ax.legend([m, v, x, a, t],
                      ['In bitmask', 'Bit masked\nand OFF',
                       'Bit unmasked\nand OFF',  'Bit ON',
                       'Transition'], **legendargs)
        else:
            ax.legend([m, v, x, a],
                      ['In bitmask', 'Bit masked\nand OFF',
                       'Bit unmasked\nand OFF', 'Bit ON'],
                      **legendargs)
        ax.set_epoch(epoch)
        ax.set_xlim(*xlim)

        # customise plot
        if ax.get_autoscaley_on():  # no user-set ylim
            self.pargs['ylim'] = (-nflags+.5, .5)
        self.apply_parameters(ax, **self.pargs)

        # add bit mask axes and finalise
        self.add_state_segments(ax)
        self.add_future_shade()
        out = self.finalize()
        return out


register_plot(ODCDataPlot)


class SegmentPiePlot(PiePlot, SegmentDataPlot):
    type = 'segment-pie'
    _single_call = True
    defaults = {
        'legend-loc': 'center left',
        'legend-bbox_to_anchor': (.8, .5),
        'legend-fontsize': 14,
        'legend-frameon': False,
        'wedge-width': .55,
        'wedge-edgecolor': 'white',
    }

    parse_plot_kwargs = TimeSeriesDataPlot.parse_plot_kwargs

    def parse_wedge_kwargs(self, defaults=dict()):
        wedgeargs = defaults.copy()
        for key in list(self.pargs):
            if key.startswith('wedge-') or key.startswith('wedge_'):
                wedgeargs[key[6:]] = self.pargs.pop(key)
        return wedgeargs

    def draw(self, outputfile=None):
        plot = self.init_plot(projection='rectilinear')
        ax = plot.gca()

        # extract plotting arguments
        future = self.pargs.pop('include_future', False)
        legendargs = self.parse_legend_kwargs()
        wedgeargs = self.parse_wedge_kwargs()
        plotargs = self.parse_plot_kwargs()

        # use state to generate suptitle with GPS span
        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        usetex_tex(str(self.state))))
        else:
            self.pargs.setdefault(
                'suptitle', '[%s-%s]' % (self.span[0], self.span[1]))

        # get segments
        data = []
        for flag in self.flags:
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            segs = get_segments(flag, validity=valid, query=False,
                                padding=self.padding).coalesce()
            data.append(float(abs(segs.active)))
        if future:
            total = sum(data)
            alltime = abs(self.span)
            data.append(alltime-total)
            if 'labels' in plotargs:
                plotargs['labels'] = list(plotargs['labels']) + [' ']
            if 'colors' in plotargs:
                plotargs['colors'] = list(plotargs['colors']) + ['white']

        # make pie
        labels = plotargs.pop('label')
        patches = ax.pie(data, **plotargs)[0]
        ax.axis('equal')

        # set wedge params
        for wedge in patches:
            for key, val in wedgeargs.items():
                getattr(wedge, 'set_%s' % key)(val)

        # make legend
        legendargs['title'] = ax.get_title()
        ax.set_title('')
        legth = legendargs.pop('threshold', 0)
        legsort = legendargs.pop('sorted', False)
        tot = float(sum(data))
        pclabels = []
        for d, label in zip(data, labels):
            if not label or label == ' ':
                pclabels.append(label)
            else:
                try:
                    pc = d/tot * 100
                except ZeroDivisionError:
                    pc = 0.0
                pclabels.append(usetex_tex(
                    '%s [%1.1f%%]' % (label, pc)).replace(r'\\', '\\'))

        # add time to top
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            extra = Rectangle((0, 0), 1, 1, fc='w', fill=False, ec='none',
                              linewidth=0)
        # sort entries
        if legsort:
            patches, pclabels, data = map(list, zip(*sorted(
                 list(zip(patches, pclabels, data)),
                 key=lambda x: x[2],
                 reverse=True)))
        # and restrict to the given threshold
        if legth:
            try:
                patches, pclabels, data = map(list, zip(*[
                    x for x in zip(patches, pclabels, data) if x[2] >= legth]))
            except ValueError:
                pass

        if suptitle:
            leg = ax.legend([extra]+patches, [suptitle]+pclabels, **legendargs)
            t = leg.get_texts()[0]
            t.set_fontproperties(t.get_fontproperties().copy())
            t.set_size(min(12, t.get_size()))
        else:
            leg = ax.legend(patches, pclabels, **legendargs)
        legt = leg.get_title()
        legt.set_fontsize(max(22, legendargs.get('fontsize', 22)+4))
        legt.set_ha('left')

        # customise plot
        self.apply_parameters(ax, **self.pargs)

        # copy title and move axes
        if ax.get_title():
            title = plot.suptitle(ax.get_title())
            title.update_from(ax.title)
            title.set_y(title._y + 0.05)
            ax.set_title('')
        axpos = ax.get_position()
        offset = -.2
        ax.set_position([axpos.x0+offset, .1, axpos.width, .8])

        # add bit mask axes and finalise
        self.pargs['xlim'] = None
        return self.finalize(outputfile=outputfile, pad_inches=0)


register_plot(SegmentPiePlot)


class NetworkDutyPiePlot(SegmentPiePlot):
    """Special case of the `SegmentPiePlot` for network duty factors
    """
    type = 'network-duty-pie'
    NETWORK_NAME = {
        0: 'no',
        1: 'single',
        2: 'double',
        3: 'triple',
        4: 'quadruple',
        5: 'quintuple',
        6: 'sextuple',
    }
    NETWORK_COLOR = GW_OBSERVATORY_COLORS.copy()
    NETWORK_COLOR.update({
        'no': 'black',
        'single': (1.0, 0.7, 0.0),
        'double': (0.0, 0.4, 1.0),
        'triple': 'pink',
        'quadruple': (1.0, 0.4, 0.0),
    })
    defaults = SegmentPiePlot.defaults.copy()
    defaults.update({
        'legend-fontsize': 24,
    })

    def draw(self):
        # get segments
        if self.state and not self.all_data:
            valid = self.state.active
        else:
            valid = SegmentList([self.span])
        # construct compound flags for each network size
        flags = dict((f[:2], f) for f in self.flags)
        network = ''.join(sorted(set(flags)))
        self.pargs.setdefault('title', '%s network duty factor' % network)
        networkflags = []
        colors = []
        labels = []
        exclude = DataQualityFlag()
        for i in list(range(len(flags)+1))[::-1]:
            name = self.NETWORK_NAME[i]
            flag = '%s:%s' % (network, name)
            networksegs = DataQualityFlag(flag, known=valid)
            for ifoset in combinations(flags, i):
                if not ifoset:
                    compound = '!%s' % '!'.join(list(flags.values()))
                else:
                    compound = '&'.join(flags[ifo] for ifo in ifoset)
                segs = get_segments(compound, validity=valid, query=False,
                                    padding=self.padding).coalesce()
                networksegs += segs
            globalv.SEGMENTS[flag] = networksegs.copy()
            globalv.SEGMENTS[flag].active -= exclude.active
            exclude = networksegs
            networkflags.append(flag)
            labels.append('%s interferometer' % name.title())
            colors.append(self.NETWORK_COLOR.get(name))

        self.pargs.setdefault('colors', colors)
        self.pargs.setdefault('labels', labels)

        # reset flags and generate plot
        flags_ = self.flags
        outputfile = self.outputfile
        self.flags = networkflags
        out = super(NetworkDutyPiePlot, self).draw(outputfile=outputfile)
        self.flags = flags_
        return out


register_plot(NetworkDutyPiePlot)


class SegmentBarPlot(BarPlot, SegmentDataPlot):
    type = 'segment-bar'
    _single_call = True
    defaults = {
        'scale': 'percent',
        'color': GREEN,
        'edgecolor': 'green',
        'alpha': .6,
        'align': 'edge',  # FIXME, updated for mpl 2.0, can simplify code
    }
    SCALE_UNIT = {
        None: 'seconds',
        1: 'seconds',
        'percent': r'\%',
        60: 'minutes',
        3600: 'hours',
    }

    def draw(self, outputfile=None):
        plot = self.init_plot(projection='rectilinear')
        ax = plot.gca()

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        usetex_tex(str(self.state))))
        else:
            self.pargs.setdefault(
                'suptitle', '[%s-%s]' % (self.span[0], self.span[1]))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.993, va='top')

        scale = self.pargs.pop('scale', 'percent')
        if scale == 'percent':
            self.pargs.setdefault('ylim', (0, 100))
        elif isinstance(scale, (int, float)):
            self.pargs.setdefault('ylim', (0, abs(self.span) / scale))
        try:
            self.pargs.setdefault('ylabel', 'Livetime [%s]'
                                  % self.SCALE_UNIT[scale])
        except KeyError:
            self.pargs.setdefault('ylabel', 'Livetime')

        # extract plotting arguments
        sort = self.pargs.pop('sorted', False)
        plotargs = self.parse_plot_kwargs()

        # get segments
        data = []
        labels = plotargs.pop('label', self.flags)
        for flag in self.flags:
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            segs = get_segments(flag, validity=valid, query=False,
                                padding=self.padding).coalesce()
            livetime = float(abs(segs.active))
            if scale == 'percent':
                try:
                    data.append(100 * livetime / float(abs(segs.known)))
                except ZeroDivisionError:
                    data.append(0)
            elif isinstance(scale, (float, int)):
                data.append(livetime / scale)

        if sort:
            data, labels = list(zip(*sorted(
                list(zip(data, labels)), key=lambda x: x[0], reverse=True)))

        # make bar chart
        width = plotargs.pop('width', .8)
        x = numpy.arange(len(data)) - width/2.
        ax.bar(x, data, width=width, **plotargs)

        # set labels
        ax.set_xticks(range(len(data)))
        ax.set_xticklabels(labels, rotation=30,
                           rotation_mode='anchor', ha='right', fontsize=13)
        ax.tick_params(axis='x', pad=2)
        ax.xaxis.labelpad = 2
        ax.xaxis.grid(False)
        self.pargs.setdefault('xlim', (-.5, len(data)-.5))

        # customise plot
        self.apply_parameters(ax, **self.pargs)

        # add bit mask axes and finalise
        self.pargs['xlim'] = None
        return self.finalize(outputfile=outputfile, transparent="True",
                             pad_inches=0)


register_plot(SegmentBarPlot)


class SegmentHistogramPlot(get_plot('histogram'), SegmentDataPlot):
    """Histogram of segment duration
    """
    type = 'segment-histogram'
    data = 'segments'
    defaults = {'ylabel': 'Number of segments',
                'log': False,
                'histtype': 'stepfilled',
                'bottom': 0,
                'rwidth': 1}

    parse_plot_kwargs = TimeSeriesDataPlot.parse_plot_kwargs

    def draw(self, outputfile=None):
        # make axes
        plot = self.init_plot(projection='rectilinear')
        axes = plot.axes

        # use state to generate suptitle with GPS span
        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        usetex_tex(str(self.state))))
        else:
            self.pargs.setdefault(
                'suptitle', '[%s-%s]' % (self.span[0], self.span[1]))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.993, va='top')

        # extract plotting arguments
        histargs = self.parse_plot_kwargs()

        # get segments
        data = []
        for flag in self.flags:
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            segs = get_segments(flag, validity=valid, query=False,
                                padding=self.padding).coalesce()
            data.append([float(abs(x)) for x in segs.active])

        # get range
        if 'range' not in histargs[0]:
            _lim = common_limits(data)
            for d in histargs:
                d['range'] = _lim

        # plot
        for ax, arr, pargs in zip(cycle(axes), data, histargs):
            if len(arr) == 0:
                kwargs = dict(
                    (k, pargs[k]) for k in ['label', 'color'] if pargs.get(k))
                ax.plot([], **kwargs)
            else:
                if pargs.get('normed', False) in ['N', 'num', 'number']:
                    pargs['normed'] = False
                    pargs.setdefault('weights', [1/len(arr)] * len(arr))
                ax.hist(arr, **pargs)

        # customise plot
        legendargs = self.parse_legend_kwargs()
        for i, ax in enumerate(axes):
            for key, val in self.pargs.items():
                if key == 'title' and i > 0:
                    continue
                if key == 'xlabel' and i < (len(axes) - 1):
                    continue
                if key == 'ylabel' and (
                        (len(axes) % 2 and i != len(axes) // 2) or
                        (len(axes) % 2 == 0 and i > 0)):
                    continue
                try:
                    getattr(ax, 'set_%s' % key)(val)
                except AttributeError:
                    setattr(ax, key, val)
            if len(self.flags) > 1:
                ax.legend(**legendargs)

        if len(axes) > 1 and axes[0].get_ylabel():
            # set text
            ylabel = axes[0].yaxis.get_label()
            y = axes[-1].get_position().y0 + (
                axes[0].get_position().y1 - axes[-1].get_position().y0)/2.
            t = plot.text(0.04, y, ylabel.get_text(), rotation=90, ha='center',
                          va='center')
            t.set_fontproperties(ylabel.get_font_properties())
            for i, ax in enumerate(axes):
                ax.set_ylabel('')
                if i:
                    ax.set_title('')
                if i < len(axes) - 1:
                    ax.set_xlabel('')
                    setp(ax.get_xticklabels(), visible=False)

        # set common ylim
        if 'ylim' not in self.pargs:
            y0 = min([ax.get_ylim()[0] for ax in axes])
            y1 = max([ax.get_ylim()[1] for ax in axes])
            for ax in axes:
                ax.set_ylim(y0, y1)

        # add bit mask axes and finalise
        return self.finalize(outputfile=outputfile, transparent="True",
                             pad_inches=0)


register_plot(SegmentHistogramPlot)
