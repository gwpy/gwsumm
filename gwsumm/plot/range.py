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

"""Definitions for range plots
"""

from __future__ import division

import re
from math import pi
from heapq import nlargest
from itertools import combinations

from six import string_types

import numpy

from matplotlib.ticker import MaxNLocator

from gwpy.segments import (Segment, SegmentList)
from gwpy.timeseries import TimeSeries

from .registry import (get_plot, register_plot)
from .utils import hash
from ..data import (get_range_channel, get_range, get_timeseries)
from ..segments import get_segments
from ..channels import split as split_channels

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


class RangePlotMixin(object):
    data = 'spectrogram'
    _threadsafe = False
    defaults = {
        'snr': 8.0,
        'stride': 60.,
        'fftlength': 8,
        'overlap': 4,
        'fmin': 10,
    }

    def __init__(self, *args, **kwargs):
        super(RangePlotMixin, self).__init__(*args, **kwargs)
        self.rangeparams = {}
        for key in ['mass1', 'mass2', 'snr', 'energy', 'stride', 'fftlength',
                    'overlap', 'fmin', 'fmax']:
            try:
                value = self.pargs.pop(key)
            except KeyError:
                continue
            if not isinstance(value, (tuple, list)):
                value = [value]*len(self.channels)
            self.rangeparams[key] = value

    def draw(self):
        """Read in all necessary data, and generate the figure.
        """
        # generate data
        keys = []
        for i, channel in enumerate(self.channels):
            kwargs = dict((key, self.rangeparams[key][i]) for
                          key in self.rangeparams if
                          self.rangeparams[key][i] is not None)
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            rlist = get_range(channel, valid, query=self.read, **kwargs)
            try:
                keys.append(rlist[0].channel)
            except IndexError:
                keys.append(get_range_channel(channel, **kwargs))

        # reset channel lists and generate time-series plot
        channels = self.channels
        outputfile = self.outputfile
        self.channels = keys
        out = super(RangePlotMixin, self).draw(outputfile=outputfile)
        self.channels = channels
        return out


class RangeDataPlot(RangePlotMixin, get_plot('timeseries')):
    type = 'range'
    defaults = get_plot('timeseries').defaults.copy()
    defaults.update(RangePlotMixin.defaults.copy())
    defaults.update({
        'ylabel': 'Sensitive distance [Mpc]',
    })


register_plot(RangeDataPlot)


class RangeDataHistogramPlot(RangePlotMixin, get_plot('histogram')):
    type = 'range-histogram'
    defaults = get_plot('histogram').defaults.copy()
    defaults.update(RangePlotMixin.defaults.copy())
    defaults.update({
        'xlabel': 'Sensitive distance [Mpc]',
    })


register_plot(RangeDataHistogramPlot)


# -- time-volume --------------------------------------------------------------

class SimpleTimeVolumeDataPlot(get_plot('segments')):
    """Time-series of the time-volume searched by an interferometer
    """
    data = 'timeseries'
    type = 'time-volume'
    DRAW_PARAMS = get_plot('timeseries').DRAW_PARAMS
    defaults = get_plot('timeseries').defaults.copy()

    parse_plot_kwargs = get_plot('timeseries').parse_plot_kwargs

    def __init__(self, sources, *args, **kwargs):
        if isinstance(sources, string_types):
            sources = split_channels(sources)
        channels = sources[::2]
        flags = sources[1::2]
        get_plot('timeseries').__init__(self, channels, *args, **kwargs)
        self._allflags = []
        self.flags = flags

    @classmethod
    def from_ini(cls, *args, **kwargs):
        return get_plot('timeseries').from_ini(cls, *args, **kwargs)

    @property
    def pid(self):
        try:
            return self._pid
        except AttributeError:
            self._pid = hash("".join(map(str, self.channels+self.flags)))
            return self.pid

    @pid.setter
    def pid(self, id_):
        self._pid = str(id_)

    @pid.deleter
    def pid(self):
        del self._pid

    @staticmethod
    def calculate_time_volume(segments, range):
        try:
            ts = TimeSeries(numpy.zeros(range.size), xindex=range.times,
                            unit='s')
        except IndexError:
            ts = TimeSeries(numpy.zeros(range.size), unit='s',
                            x0=range.x0, dx=range.dx)
        dx = range.dx.value

        # override range units
        range.override_unit('Mpc')

        # use float, not LIGOTimeGPS for speed
        segments = type(segments)([type(s)(float(s[0]), float(s[1])) for
                                  s in segments])

        def livetime_(t):
            return float(abs(SegmentList([Segment(t, t+dx)]) & segments))

        livetime = numpy.vectorize(livetime_, otypes=[float])
        ts[:] = livetime(ts.times.value) * ts.unit
        return (4/3. * pi * ts * range ** 3).to('Mpc^3 kyr')

    def combined_time_volume(self, allsegments, allranges):
        try:
            combined_range = TimeSeries(numpy.zeros(allranges[0].size),
                                        xindex=allranges[0].times, unit='Mpc')
        except IndexError:
            combined_range = TimeSeries(
                numpy.zeros(allranges[0].size), unit='Mpc',
                x0=allranges[0].x0, dx=allranges[0].dx)

        # get coincident observing segments
        pairs = list(combinations(allsegments, 2))
        coincident = SegmentList()
        for pair in pairs:
            coincident.extend(pair[0] & pair[1])
        coincident = coincident.coalesce()

        # get effective network range
        values = [r.value for r in allranges]
        values = [min(nlargest(2, x)) for x in zip(*values)]
        size = min([r.size for r in allranges])
        combined_range[:size] = values * combined_range.unit

        # compute time-volume
        return self.calculate_time_volume(coincident, combined_range)

    def draw(self, outputfile=None):
        """Generate the figure for this plot
        """
        plot = self.init_plot()
        ax = plot.axes[0]

        # get plotting arguments
        cumulative = self.pargs.pop('cumulative', False)
        plotargs = self.parse_plot_kwargs()
        legendargs = self.parse_legend_kwargs()

        # set ylabel
        if cumulative:
            self.pargs.setdefault('ylabel',
                                  'Cumulative time-volume [Mpc$^3$ kyr]')
        else:
            self.pargs.setdefault('ylabel', 'Time-volume [Mpc$^3$ kyr]')

        # get data
        allsegs, allranges = ([], [])
        for channel, flag, pargs in zip(self.channels, self.flags, plotargs):
            pad = 0
            if self.state and not self.all_data:
                valid = self.state.active
                pad = numpy.nan
            elif channel.sample_rate:
                valid = SegmentList([self.span.protract(
                    1/channel.sample_rate.value)])
            else:
                valid = SegmentList([self.span])
            data = get_timeseries(
                channel, valid, query=False).join(gap='pad', pad=pad)
            if not data.unit or data.unit.to_string() in ['', 'undef']:
                data.override_unit('Mpc')
            segments = get_segments(flag, valid, query=False)
            timevolume = self.calculate_time_volume(segments.active, data)
            if cumulative:
                ax.plot(timevolume.cumsum(), **pargs)
            else:
                ax.plot(timevolume, **pargs)
            allsegs.append(segments.active)
            allranges.append(data)

        # estimate combined time-volume
        if self.all_data and len(self.channels) > 1:
            pargs = plotargs[-1]
            pargs['color'] = '#000000'
            pargs['label'] = 'Combined'
            pargs['linestyle'] = '--'
            combined_timevolume = self.combined_time_volume(
                allsegs, allranges)
            if cumulative:
                ax.plot(combined_timevolume.cumsum(), **pargs)
            else:
                ax.plot(combined_timevolume, **pargs)

        # add horizontal lines to add
        for yval in self.pargs.get('hline', []):
            try:
                yval = float(yval)
            except ValueError:
                continue
            else:
                ax.plot([self.start, self.end], [yval, yval],
                        linestyle='--', color='red')

        # customise plot
        self.apply_parameters(ax, **self.pargs)
        if (len(self.channels) > 1 or plotargs[0].get('label', None) in
                [re.sub(r'(_|\\_)', r'\_', str(self.channels[0])), None]):
            ax.legend(**legendargs)

        # add extra axes and finalise
        self.add_state_segments(ax)
        self.add_future_shade()
        ax.yaxis.set_major_locator(MaxNLocator(8))
        ticks = ax.get_yticks()
        ax.yaxis.set_ticklabels(ticks)
        return self.finalize(outputfile=outputfile)


register_plot(SimpleTimeVolumeDataPlot)


class GWpyTimeVolumeDataPlot(RangePlotMixin, SimpleTimeVolumeDataPlot):
    """TimeVolumeDataPlot where the range is calculated on-the-fly
    """
    type = 'strain-time-volume'
    _threadsafe = False


register_plot(GWpyTimeVolumeDataPlot)
