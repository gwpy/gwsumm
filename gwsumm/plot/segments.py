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

import hashlib
import bisect
from itertools import (cycle, combinations)

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

import numpy

from dateutil.relativedelta import relativedelta

from matplotlib.patches import Rectangle

from gwpy.plotter import *
from gwpy.plotter.tex import label_to_latex
from gwpy.segments import (Segment, SegmentList, DataQualityFlag)
from gwpy.time import (from_gps, to_gps)

from .. import globalv
from ..mode import (Mode, get_mode)
from ..config import NoOptionError
from ..utils import (re_quote, get_odc_bitmask, re_flagdiv, safe_eval)
from ..data import (get_channel, get_timeseries)
from ..segments import (get_segments, format_padding)
from ..state import ALLSTATE
from .core import (BarPlot, PiePlot)
from .registry import (get_plot, register_plot)
from .mixins import *

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

TimeSeriesDataPlot = get_plot('timeseries')
GREEN = (0.2, 0.8, 0.2)


class SegmentDataPlot(SegmentLabelSvgMixin, TimeSeriesDataPlot):
    """Segment plot of one or more `DataQualityFlags <DataQualityFlag>`.
    """
    type = 'segments'
    data = 'segments'
    defaults = {'mask': None,
                'color': None,
                'on-is-bad': False,
                'insetlabels': 'inset',
                'legend-bbox_to_anchor': (1.01, 1.),
                'legend-loc': 'upper left',
                'legend-borderaxespad': 0,
                'legend-fontsize': 12}
    DRAW_PARAMS = TimeSeriesDataPlot.DRAW_PARAMS + ['known']

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
        if isinstance(flist, str):
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
        for f, p in format_padding(self._allflags, pad).iteritems():
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
            self._pid = hashlib.md5(
                "".join(map(str, self.flags))).hexdigest()[:6]
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
        if isinstance(flags, str):
            flags = [f.strip('\n ') for f in flags.split(',')]
        new.flags = flags
        return new

    def get_segment_color(self):
        """Parse the configured ``pargs`` and determine the colors for
        active and valid segments.
        """
        active = safe_eval(
            self.pargs.pop('active', self.pargs.pop('facecolor', None)))
        known = safe_eval(self.pargs.pop('known', 0))
        # neither known nor active defined
        if active is None and known == 0:
            if bool(self.pargs.pop('on-is-bad', False)):
                self.pargs['facecolor'] = 'red'
                self.pargs.setdefault('edgecolor', 'darkred')
                self.pargs['known'] = {'facecolor': GREEN}
            else:
                self.pargs['facecolor'] = GREEN
                self.pargs.setdefault('edgecolor', 'darkred')
                self.pargs['known'] = 'red'
        # only active is defined
        elif known == 0:
            if isinstance(active, dict):
                self.pargs.update(active)
                active = active.get('facecolor')
            else:
                self.pargs['facecolor'] = active
            if isinstance(active, str) and active.lower() == 'red':
                self.pargs['known'] = 'dodgerblue'
            else:
                self.pargs['known'] = 'red'
        # only known is defined
        elif active is None:
            self.pargs['known'] = known
            if known in [GREEN, str(GREEN), 'green', 'g']:
                self.pargs['facecolor'] = 'dodgerblue'
                self.pargs.setdefault('edgecolor', 'blue')
            else:
                self.pargs['facecolor'] = GREEN
                self.pargs.setdefault('edgecolor', 'green')
        # both are given
        else:
            if isinstance(active, dict):
                self.pargs.update(active)
            else:
                self.pargs['facecolor'] = active
            self.pargs['known'] = known
        self.pargs.setdefault('edgecolor', 'black')
        return self.pargs

    def parse_plot_kwargs(self, *args, **kwargs):
        self.get_segment_color()
        return super(SegmentDataPlot, self).parse_plot_kwargs(*args, **kwargs)

    def draw(self):
        # get labelsize
        _labelsize = rcParams['ytick.labelsize']
        labelsize = self.pargs.pop('labelsize', 12)
        if self.pargs.get('insetlabels', True) is False:
            rcParams['ytick.labelsize'] = labelsize

        # create figure
        (plot, axes) = self.init_plot(plot=SegmentPlot)
        ax = axes[0]

        # extract plotting arguments
        legendargs = self.parse_legend_kwargs()
        mask = self.pargs.pop('mask')
        plotargs = self.parse_plot_kwargs()
        legcolors = plotargs[0].copy()

        # plot segments
        for i, (flag, pargs) in enumerate(
                zip(self.flags, plotargs)[::-1]):
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
            pargs.setdefault('known', None)
            ax.plot(segs, y=i, label=label, **pargs)

        # make custom legend
        known = legcolors.pop('known', None)
        if known:
            active = legcolors.pop('facecolor')
            edgecolor = legcolors.pop('edgecolor')
            epoch = ax.get_epoch()
            xlim = ax.get_xlim()
            seg = SegmentList([Segment(self.start - 10, self.start - 9)])
            if isinstance(known, dict):
                known = known['facecolor']
            v = ax.plot(seg, facecolor=known, collection=False)[0][0]
            a = ax.plot(seg, facecolor=active, edgecolor=edgecolor,
                        collection=False)[0][0]
            if edgecolor not in [None, 'none']:
                t = ax.plot(seg, facecolor=edgecolor, collection=False)[0][0]
                ax.legend([v, a, t], ['Known', 'Active', 'Transition'],
                          **legendargs)
            else:
                ax.legend([v, a], ['Known', 'Active'], **legendargs)
            ax.set_epoch(epoch)
            ax.set_xlim(*xlim)

        # customise plot
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if 'ylim' not in self.pargs:
            ax.set_ylim(-0.5, len(self.flags) - 0.5)

        # add bit mask axes and finalise
        if mask is None and not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        elif mask is not None:
            plot.add_bitmask(mask, topdown=True)
        if self.state and self.state.name != ALLSTATE:
            self.add_state_segments(ax)

        rcParams['ytick.labelsize'] = _labelsize
        return self.finalize()

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
        except:
            chans = "".join(map(str, self.channels))
            self._pid = hashlib.md5(chans).hexdigest()[:6]
            if self.pargs.get('bits', None):
                self._pid = hashlib.md5(
                    self._pid + str(self.pargs['bits'])).hexdigest()[:6]
            return self.pid

    def _parse_labels(self, defaults=[]):
        """Pop the labels for plotting from the `pargs` for this Plot

        This method overrides from the `TimeSeriesDataPlot` in order
        to set the bit names from the various channels as the defaults
        in stead of the channel names
        """
        chans = zip(*self.get_channel_groups())[0]
        labels = list(self.pargs.pop('labels', defaults))
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        for i, l in enumerate(labels):
            if isinstance(l, (list, tuple)):
                labels[i] = list(labels[i])
                for j, l2 in enumerate(l):
                    labels[i][j] = rUNDERSCORE.sub(r'\_', str(l2).strip('\n '))
            elif isinstance(l, str):
                labels[i] = rUNDERSCORE.sub(r'\_', str(l).strip('\n '))
        while len(labels) < len(chans):
            labels.append(None)
        return labels

    def parse_plot_kwargs(self, *args, **kwargs):
        self.get_segment_color()
        return super(StateVectorDataPlot, self).parse_plot_kwargs(
            *args, **kwargs)

    def draw(self):
        # make font size smaller
        labelsize = self.rcParams.get('ytick.labelsize', 12)
        if self.pargs.get('insetlabels', True) is False:
            rcParams['ytick.labelsize'] = labelsize

        (plot, axes) = self.init_plot(plot=SegmentPlot)
        ax = axes[0]

        # get bit setting
        bits = self.pargs.pop('bits', None)
        if bits and len(self.channels) > 1:
            raise ValueError("Specifying 'bits' doesn't work for a "
                             "state-vector plot including multiple channels")

        # extract plotting arguments
        mask = self.pargs.pop('mask')
        ax.set_insetlabels(self.pargs.pop('insetlabels', True))
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
                bits_ = channel.bits
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
                if not 'int' in str(stateseries.dtype):
                    stateseries = stateseries.astype('uint32')
                newflags = stateseries.to_dqflags().values()
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
            for flag, label in zip(flags, labels)[::-1]:
                kwargs = pargs.copy()
                if label is not None:
                    kwargs['label'] = label
                ax.plot(flag, **kwargs)

        # customise plot
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if 'ylim' not in self.pargs:
            ax.set_ylim(-0.5, nflags - 0.5)

        # add bit mask axes and finalise
        if mask is None and not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        elif mask is not None:
            plot.add_bitmask(mask, topdown=True)
        if self.state and self.state.name != ALLSTATE:
            self.add_state_segments(ax)

        return self.finalize()

register_plot(StateVectorDataPlot)


class DutyDataPlot(SegmentDataPlot):
    """`DataPlot` of the duty-factor for a `SegmentList`
    """
    type = 'duty'
    data = 'segments'
    defaults = {'alpha': 0.8,
                'sep': False,
                'side_by_side': False,
                'normalized': None,
                'cumulative': False,
                'stacked': False,
                'ylabel': r'Duty factor [\%]'}

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
        except:
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
                                 '\'side_by_side\' should not be used together')
            geometry = (len(self.flags), 1)
        else:
            geometry = (1, 1)

        (plot, axes) = self.init_plot(plot=TimeSeriesPlot, geometry=geometry)

        # extract plotting arguments
        style = self.pargs.pop('style', 'bar')
        stacked = self.pargs.pop('stacked', False)
        sidebyside = self.pargs.pop('side_by_side', False)
        normalized = self.pargs.pop('normalized', True)
        cumulative = self.pargs.pop('cumulative', False)
        if normalized is None and not cumulative:
            normalized = 'percent'
        plotargs = self.parse_plot_kwargs()
        legendargs = self.parse_legend_kwargs()
        if sep:
            legendargs.setdefault('loc', 'upper left')
            legendargs.setdefault('bbox_to_anchor', (1.01, 1))
            legendargs.setdefault('borderaxespad', 0)
        rollingmean = self.pargs.pop('rolling_mean',
                                     not stacked and not cumulative)

        # work out times and plot mean for legend
        self.get_bins()
        times = float(self.start) + numpy.concatenate(
                                 ([0], self.bins[:-1].cumsum()))
        now = bisect.bisect_left(times, globalv.NOW)
        if rollingmean:
            axes[0].plot(times[:1], [-1], 'k--', label='Rolling mean')

        # get bar parameters
        try:
            bottom = self.pargs['ylim'][0]
        except KeyError:
            bottom = 0
        bottom = numpy.zeros(times.size) + bottom

        # plot segments
        if self.state and not self.all_data:
            valid = self.state.active
        else:
            valid = SegmentList([self.span])
        for i, (ax, flag, pargs, color) in enumerate(
                zip(cycle(axes), self.flags, plotargs,
                    cycle(rcParams['axes.color_cycle']))):
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
            color = pargs.pop('color', color)
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
                    print(w, offset)
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
                    ax.plot(times[:now], duty[:now], drawstyle='steps-post',
                            color=ec, linewidth=lw)

            # plot mean
            if rollingmean:
                t = [self.start] + list(times + self.bins/2.) + [self.end]
                mean = [mean[0]] + list(mean) + [mean[-1]]
                ax.plot(t, mean, color=sep and 'k' or color, linestyle='--')

            # record duty for stacked chart
            if stacked:
                bottom += height

        # customise plot
        for key, val in self.pargs.iteritems():
            for ax in axes:
                try:
                    getattr(ax, 'set_%s' % key)(val)
                except AttributeError:
                    setattr(ax, key, val)
        if 'hours' in self.pargs.get('ylabel', ''):
            ax.get_yaxis().get_major_locator().set_params(
                steps=[1, 2, 4, 8, 12, 24])
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

        # add custom legend for mean
        if rollingmean:
            yoff = 0.01 * float.__div__(*axes[0].get_position().size)
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

        for ax in axes:
            try:
                plot.add_legend(ax=ax, **legendargs)
            except AttributeError:
                pass

        # add extra axes and finalise
        if not plot.colorbars:
            for ax in axes:
                plot.add_colorbar(ax=ax, visible=False)
        if self.state:
            self.add_state_segments(axes[-1])

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
            self.bitmask = map(get_odc_bitmask, self.channels)

    def get_bitmask_channels(self):
        return type(self.channels)(list(map(get_channel, self.bitmask)))

    @property
    def pid(self):
        try:
            return self._pid
        except:
            chans = "".join(map(str, self.channels))
            masks = "".join(map(str, self.get_bitmask_channels()))
            self._pid = hashlib.md5(chans+masks).hexdigest()[:6]
            if self.pargs.get('bits', None):
                self._pid = hashlib.md5(
                    self._pid + str(self.pargs['bits'])).hexdigest()[:6]
            return self.pid

    def draw(self):
        # make font size smaller
        _labelsize = rcParams['ytick.labelsize']
        labelsize = self.pargs.pop('labelsize', 12)
        rcParams['ytick.labelsize'] = labelsize

        # make figure
        (plot, axes) = self.init_plot(plot=SegmentPlot)
        ax = axes[0]
        ax.grid(False, which='both', axis='y')

        # extract plotting arguments
        ax.set_insetlabels(self.pargs.pop('insetlabels', True))
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
                    if not 'int' in str(stateseries.dtype):
                        stateseries = stateseries.astype('uint32')
                    newflags = stateseries.to_dqflags()
                    if flags[type_] is None:
                        flags[type_] = newflags
                    else:
                        for i, flag in newflags.iteritems():
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
        m = ax.build_segment(seg, y=0, facecolor=inmaskcolor, edgecolor='none')
        v = ax.build_segment(seg, y=0, facecolor=maskoncolor,
                             edgecolor=edgecolor)
        x = ax.build_segment(seg, y=0, facecolor=maskoffcolor,
                             edgecolor=edgecolor)
        a = ax.build_segment(seg, y=0, facecolor=activecolor,
                             edgecolor=edgecolor)
        if edgecolor not in [None, 'none']:
            t = ax.build_segment(seg, y=0, facecolor=edgecolor)
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
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if 'ylim' not in self.pargs:
            ax.set_ylim(-nflags+0.5, 0.5)

        # add bit mask axes and finalise
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        if self.state and self.state.name != ALLSTATE:
            self.add_state_segments(ax)
        out = self.finalize()
        rcParams['ytick.labelsize'] = _labelsize
        return out

register_plot(ODCDataPlot)


class SegmentPiePlot(PiePlot, SegmentDataPlot):
    type = 'segment-pie'
    defaults = {
        'legend-loc': 'center left',
        'legend-bbox_to_anchor': (.8, .5),
        'legend-frameon': False,
        'wedge-width': .55,
        'wedge-edgecolor': 'white',
    }

    def init_plot(self, plot=Plot, geometry=(1,1)):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesDataPlot`.
        """
        figsize = self.pargs.pop('figsize', [12, 6])
        self.plot = Plot(figsize=figsize)
        axes = [self.plot.gca()]
        return self.plot, axes

    def parse_wedge_kwargs(self, defaults=dict()):
        wedgeargs = defaults.copy()
        for key in self.pargs.keys():
            if key.startswith('wedge-') or key.startswith('wedge_'):
                wedgeargs[key[6:]] = self.pargs.pop(key)
        return wedgeargs

    def draw(self, outputfile=None):
        (plot, axes) = self.init_plot(plot=Plot)
        ax = axes[0]

        # get labels
        #flags = map(lambda f: str(f).replace('_', r'\_'), self.flags)
        #labels = self.pargs.pop('labels', self.pargs.pop('label', flags))
        #labels = map(lambda s: re_quote.sub('', str(s).strip('\n ')), labels)

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
                                        label_to_latex(str(self.state))))
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
        labels = plotargs.pop('labels')
        patches = ax.pie(data, **plotargs)[0]
        ax.axis('equal')

        # set wedge params
        for wedge in patches:
            for key, val in wedgeargs.iteritems():
                getattr(wedge, 'set_%s' % key)(val)

        # make legend
        legendargs['title'] = self.pargs.pop('title', None)
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
                pclabels.append(label_to_latex(
                    '%s [%1.1f%%]' % (label, pc)).replace(r'\\', '\\'))

        # add time to top
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            extra = Rectangle((0,0), 1, 1, fc='w', fill=False, ec='none',
                              linewidth=0)
        # sort entries
        if legsort:
            patches, pclabels, data = map(list, zip(*sorted(
                 zip(patches, pclabels, data),
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
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)

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
        return self.finalize(outputfile=outputfile, transparent="True",
                             pad_inches=0)

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
    NETWORK_COLOR = {
        'H1': 'red',
        'L1': (0.2, 0.8, 0.2),
        'V1': (0.5, 0., 0.75),
        'G1': 'gray',
        'no': 'black',
        'single': (1.0, 0.7, 0.0),
        'double': (0.0, 0.4, 1.0),
        'triple': 'pink',
        'quadruple': (1.0, 0.4, 0.0),
    }
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
        flags = dict((f[:2],f) for f in self.flags)
        network = ''.join(sorted(set(flags.keys())))
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
                    compound = '!%s' % '!'.join(flags.values())
                else:
                    compound = '&'.join(flags[ifo] for ifo in ifoset)
                segs = get_segments(compound, validity=valid, query=False,
                                    padding=self.padding).coalesce()
                networksegs += segs
            globalv.SEGMENTS[flag] = networksegs - exclude
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
    defaults = {
        'edgecolor': 'white',
        'scale': 'percent',
        'color': GREEN,
        'edgecolor': 'green',
        'alpha': .6,
    }
    SCALE_UNIT = {
        None: 'seconds',
        1: 'seconds',
        'percent': r'\%',
        60: 'minutes',
        3600: 'hours',
    }

    def init_plot(self, plot=Plot, geometry=(1,1)):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesDataPlot`.
        """
        figsize = self.pargs.pop('figsize', [12, 6])
        self.plot = Plot(figsize=figsize)
        axes = [self.plot.gca()]
        return self.plot, axes

    def draw(self, outputfile=None):
        (plot, axes) = self.init_plot(plot=Plot)
        ax = axes[0]

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        label_to_latex(str(self.state))))
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
        labels = plotargs.pop('labels', self.flags)
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
            data, labels = zip(*sorted(
                zip(data, labels), key=lambda x: x[0], reverse=True))

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
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)

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

    def parse_plot_kwargs(self, *args, **kwargs):
        return super(SegmentDataPlot, self).parse_plot_kwargs(*args, **kwargs)

    def draw(self, outputfile=None):
        # make axes
        (plot, axes) = self.init_plot()

        # use state to generate suptitle with GPS span
        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        label_to_latex(str(self.state))))
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
            data.append(map(lambda x: float(abs(x)), segs.active))

        # get range
        if not 'range' in histargs[0]:
            l = axes[0].common_limits(data)
            for d in histargs:
                d['range'] = l

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
            for key, val in self.pargs.iteritems():
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
                plot.add_legend(ax=ax, **legendargs)
        if len(axes) % 2 == 0 and axes[0].get_ylabel():
            label = axes[0].yaxis.label
            ax = axes[int(len(axes) // 2)-1]
            ax.set_ylabel(label.get_text())
            ax.yaxis.label.set_position((0, -.2 / len(axes)))
            if len(axes) != 2:
                label.set_text('')

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
