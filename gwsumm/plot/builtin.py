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

import hashlib

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from gwpy.plotter import *

from .. import version
from ..utils import (re_cchar, split_channels)
from ..data import (get_channel, get_timeseries, get_spectrogram, get_spectrum)
from ..segments import get_segments
from ..triggers import get_triggers
from ..state import ALLSTATE
from .core import DataSummaryPlot
from .registry import register_plot

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

__all__ = ['TimeSeriesSummaryPlot', 'SegmentSummaryPlot',
           'SpectrumSummaryPlot', 'SpectrogramSummaryPlot',
           'StateVectorSummaryPlot', 'TriggerSummaryPlot']


class TimeSeriesSummaryPlot(DataSummaryPlot):
    """SummaryPlot of some `TimeSeries` data.
    """
    type = 'timeseries'
    defaults = {'logy': False,
                'hline': list()}

    def add_state_segments(self, ax, **kwargs):
        """Add an `Axes` below the given ``ax`` displaying the `SummaryState`
        for this `TimeSeriesSummaryPlot`.

        Parameters
        ----------
        ax : `Axes`
            the set of `Axes` below which to display the state segments.
        **kwargs
            other keyword arguments will be passed to the
            :meth:`~gwpy.plotter.timeseries.TimeSeriesPlot.add_state_segments`
            method.
        """
        kwargs.setdefault('edgecolor', 'black')
        kwargs.setdefault('facecolor', 'green')
        sax = self.plot.add_state_segments(self.state, ax, plotargs=kwargs)
        sax.tick_params(axis='y', which='major', labelsize=12)
        sax.set_epoch(self.start)
        sax.auto_gps_scale(self.end - self.start)
        return sax

    def init_plot(self, plot=TimeSeriesPlot):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesSummaryPlot`.
        """
        self.plot = plot()
        ax = self.plot.gca()
        ax.set_epoch(self.start)
        ax.auto_gps_scale(self.end-self.start)
        ax.set_xlim(self.start, self.end)
        return self.plot, ax

    def process(self):
        """Read in all necessary data, and generate the figure.
        """
        (plot, ax) = self.init_plot()

        # work out labels
        mmmchans = self.get_channel_groups()
        labels = self.pargs.pop('labels', mmmchans.keys())
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        labels = map(lambda s: str(s).strip('\n '), labels)

        # add data
        for label, channels in zip(labels, zip(*mmmchans.items())[1]):
            data = [get_timeseries(c, self.state,
                                   query=False).join(pad=numpy.nan)
                    for c in channels]
            # double-check empty
            if not 'x0' in data[0].metadata:
                data[0].epoch = self.start
            # double-check log scales
            if self.pargs['logy']:
                for ts in data:
                    ts[ts.data == 0] = 1e-100
            # plot groups or single TimeSeries
            if len(channels) > 1:
                ax.plot_timeseries_mmm(*data, label=label.replace('_', r'\_'))
            else:
                ax.plot_timeseries(data[0], label=label.replace('_', r'\_'))

            # allow channel data to set parameters
            if hasattr(data[0].channel, 'amplitude_range'):
                self.pargs.setdefault('ylim',
                                      data[0].channel.amplitude_range)

        # add horizontal lines to add
        for yval in self.pargs['hlines']:
            try:
                yval = float(yval)
            except ValueError:
                continue
            else:
                ax.plot([self.start, self.end], [yval, yval],
                        linestyle='--', color='red')

        # customise plot
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if len(mmmchans) > 1:
            plot.add_legend(ax=ax)

        # add extra axes and finalise
        if not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)
        self.add_state_segments(ax)
        self.finalize()


class SpectrogramSummaryPlot(TimeSeriesSummaryPlot):
    """SummaryPlot a Spectrogram
    """
    type = 'spectrogram'
    defaults = {'ratio': None,
                'format': None,
                'clim': None,
                'logcolor': False,
                'colorlabel': None}

    @property
    def tag(self):
        try:
            return self._tag
        except AttributeError:
            tag = super(SpectrogramSummaryPlot, self).tag
            ratio = self.pargs['ratio']
            if ratio:
                tag += '_%s_RATIO' % re_cchar.sub('_', ratio.upper())
            return tag

    def process(self):
        # initialise
        (plot, ax) = self.init_plot()
        ax.grid(b=True, axis='y', which='major')

        # parse data arguments
        sdform = self.pargs.pop('format')
        clim = self.pargs.pop('clim')
        clog = self.pargs.pop('logcolor')
        clabel = self.pargs.pop('colorlabel')
        ratio = self.pargs.get('ratio')

        # get data
        specgrams = get_spectrogram(self.channels[0], self.state, query=False,
                                    format=sdform)
        # calculate ratio spectrum
        if ratio in ['median', 'mean']:
            allspec = specgrams.join(gap='ignore')
            ratio = getattr(allspec, ratio)(axis=0)
        # plot data
        for specgram in specgrams:
            if ratio is not None:
                specgram = specgram.ratio(ratio)
            ax.plot_spectrogram(specgram)

            # allow channel data to set parameters
            if hasattr(specgram.channel, 'frequency_range'):
                self.pargs.setdefault('ylim', specgram.channel.frequency_range)
            if (ratio is None and sdform in ['amplitude', 'asd'] and
                    hasattr(specgram.channel, 'asd_range') and clim is None):
                clim = specgram.channel.asd_range
            elif (ratio is None and hasattr(specgram.channel, 'psd_range') and
                  clim is None):
                clim = specgram.channel.psd_range

        # add colorbar
        if len(specgrams) == 0:
            ax.scatter([1], [1], c=[1], visible=False)
        plot.add_colorbar(ax=ax, clim=clim, log=clog, label=clabel)

        # customise and finalise
        for key, val in self.pargs.iteritems():
            if key == 'ratio':
                continue
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        self.add_state_segments(ax)
        self.finalize()


class SegmentSummaryPlot(TimeSeriesSummaryPlot):
    """Segment plot of one or more `DataQualityFlags <DataQualityFlag>`.
    """
    type = 'segments'
    defaults = {'mask': None,
                'color': None,
                'on_is_bad': False,
                'add_label': 'inset',
                'edgecolor': 'black'}

    def __init__(self, flags, start, end, state=None, outdir='.', tag=None,
                 **kwargs):
        super(SegmentSummaryPlot, self).__init__([], start, end, state=state,
                                                 outdir=outdir, tag=tag,
                                                 **kwargs)
        self.flags = flags

    @property
    def flags(self):
        return self._flags

    @flags.setter
    def flags(self, flist):
        self._flags = flist

    @property
    def ifos(self):
        """Interferometer set for this `SegmentSummaryPlot`
        """
        allflags = [f for flag in self.flags for f in flag.split(',')]
        return set([f[:2] for f in allflags])

    @property
    def tag(self):
        """File tag for this `DataSummaryPlot`.
        """
        try:
            return self._tag
        except AttributeError:
            state = re_cchar.sub('_',
                                 self.state is None and 'MULTI' or
                                 self.state.name)
            hash_ = hashlib.md5("".join(map(str, self.flags))).hexdigest()[:6]
            return '_'.join(map(str.upper, [state, hash_, self.type]))

    @classmethod
    def from_ini(cls, config, section, start, end, flags=None, state=ALLSTATE,
                 **kwargs):
        new = super(SegmentSummaryPlot, cls).from_ini(config, section, start,
                                                      end, state=state,
                                                      **kwargs)
        # get flags
        if flags is None:
            flags = dict(config.items(section)).pop('flags', '')
        new.flags.append(split_channels(flags))
        return new

    def get_segment_color(self):
        """Parse the configured ``pargs`` and determine the colors for
        active and valid segments.
        """
        color = self.pargs.pop('color')
        onisbad = self.pargs.pop('on_is_bad')
        # allow lazy configuration
        if onisbad is None:
            onisbad = True
        else:
            onisbad = bool(onisbad)

        # choose colors
        good = color or 'green'
        if color is None or color.lower() != 'red':
            bad = 'red'
        else:
            bad = 'blue'
        if onisbad is False:
            return good, bad
        else:
            return bad, good

    def process(self):
        (plot, ax) = self.init_plot(plot=SegmentPlot)

        # get labels
        flags = map(lambda f: str(f).replace('_', r'\_'), self.flags)
        labels = self.pargs.pop('labels', self.pargs.pop('label', flags))
        addlabel = self.pargs.pop('add_label')
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        labels = map(lambda s: str(s).strip('\n '), labels)

        # extract plotting arguments
        mask = self.pargs.pop('mask')
        activecolor, validcolor = self.get_segment_color()
        edgecolor = self.pargs.pop('edgecolor')
        plotargs = {'add_label': addlabel,
                    'facecolor': activecolor,
                    'edgecolor': edgecolor,
                    'valid': {'facecolor': validcolor}}

        # plot segments
        for flag, label in zip(self.flags, labels)[::-1]:
            if self.state:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            segs = get_segments(flag, validity=valid, query=False)
            ax.plot(segs, label=label, **plotargs)

        # customise plot
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if 'ylim' not in self.pargs:
            ax.set_ylim(-0.5, len(self.flags) - 0.5)

        # add bit mask axes and finalise
        if mask is None and not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)
        elif mask is not None:
            plot.add_bitmask(mask, topdown=True)
        self.finalize()


class StateVectorSummaryPlot(TimeSeriesSummaryPlot):
    """SummaryPlot of some `StateVector` data.

    While technically a sub-class of the `TimeSeriesSummaryPlot`, for
    data access and processing reasons, the output shadows that of the
    `SegmentSummaryPlot` more closely.
    """
    type = 'statevector'
    defaults = SegmentSummaryPlot.defaults.copy()

    # copy from SegmentSummaryPlot
    flag = property(fget=SegmentSummaryPlot.flags.__get__,
                    fset=SegmentSummaryPlot.flags.__set__,
                    fdel=SegmentSummaryPlot.flags.__delete__,
                    doc="""List of flags generated for this
                        `StateVectorSummaryPlot`.""")
    get_segment_color = SegmentSummaryPlot.__dict__['get_segment_color']

    def __init__(self, *args, **kwargs):
        super(StateVectorSummaryPlot, self).__init__(*args, **kwargs)
        self.flags = []

    def process(self):
        (plot, ax) = self.init_plot(plot=SegmentPlot)

        # extract plotting arguments
        mask = self.pargs.pop('mask')
        addlabel = self.pargs.pop('add_label')
        activecolor, validcolor = self.get_segment_color()
        edgecolor = self.pargs.pop('edgecolor')
        plotargs = {'add_label': addlabel,
                    'facecolor': activecolor,
                    'edgecolor': edgecolor,
                    'valid': {'facecolor': validcolor}}

        # plot segments
        nflags = 0
        for channel in self.channels[::-1]:
            if self.state:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            data = get_timeseries(str(channel), valid, query=False,
                                  statevector=True)
            flags = None
            for stateseries in data:
                if not stateseries.size:
                    stateseries.epoch = self.start
                    stateseries.dx = 0
                    if channel.sample_rate is not None:
                        stateseries.sample_rate = channel.sample_rate
                stateseries.bitmask = channel.bitmask
                newflags = stateseries.to_dqflags()
                if flags is None:
                    flags = newflags
                else:
                    for i, flag in enumerate(newflags):
                        flags[i] += flag
            nflags += len([m for m in channel.bitmask if m is not None])
            ax.plot(*flags[::-1], **plotargs)

        # customise plot
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if 'ylim' not in self.pargs:
            ax.set_ylim(-0.5, nflags - 0.5)

        # add bit mask axes and finalise
        if mask is None and not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)
        elif mask is not None:
            plot.add_bitmask(mask, topdown=True)
        self.finalize()


class SpectrumSummaryPlot(DataSummaryPlot):
    """Spectrum plot for a `SummaryTab`
    """
    type = 'spectrum'
    defaults = {'logx': True,
                'logy': True,
                'color': [],
                'format': None}

    def process(self):
        """Load all data, and generate this `SpectrumSummaryPlot`
        """
        plot = self.plot = SpectrumPlot()
        ax = plot.gca()

        # get spectrum format: 'amplitude' or 'power'
        sdform = self.pargs.pop('format')

        # get colors
        colors = self.pargs.pop('color', [])
        if isinstance(colors, str):
            colors = colors.split(',')
        while len(colors) < len(self.channels):
            colors.append(None)

        # get labels
        labels = self.pargs.pop('labels', map(str, self.channels))
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        labels = map(lambda s: str(s).strip('\n '), labels)

        # add data
        for label, color, channel in zip(labels, colors, self.channels):
            data = get_spectrum(str(channel), self.state, query=False,
                                format=sdform)

            # anticipate log problems
            if self.pargs['logx']:
                data = [s[1:] for s in data]
            if self.pargs['logy']:
                for sp in data:
                    sp[sp.data == 0] = 1e-100

            if color is not None:
                ax.plot_spectrum_mmm(*data, label=label.replace('_', r'\_'),
                                     color=color)
            else:
                ax.plot_spectrum_mmm(*data, label=label.replace('_', r'\_'))

            # allow channel data to set parameters
            if hasattr(data[0].channel, 'frequency_range'):
                self.pargs.setdefault('xlim', data[0].channel.frequency_range)
            if (sdform in ['amplitude', 'asd'] and
                    hasattr(data[0].channel, 'asd_range')):
                self.pargs.setdefault('ylim', data[0].channel.asd_range)
            elif hasattr(data[0].channel, 'psd_range'):
                self.pargs.setdefault('ylim', data[0].channel.psd_range)

        # customise
        legendargs = {}
        for key, val in self.pargs.iteritems():
            if re.match('legend-', key):
                legendargs[key[7:]] = val
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if len(self.channels) > 1:
            plot.add_legend(ax=ax, **legendargs)
        if not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)

        # quick fix for x-axis labels hitting the axis
        for ax in plot.axes:
            ax.tick_params(axis='x', pad=10)
            ax.xaxis.labelpad = 10

        self.finalize()


class TriggerSummaryPlot(TimeSeriesSummaryPlot):
    type = 'triggers'
    defaults = {'x': 'time',
                'y': 'snr',
                'color': None,
                'edgecolor': None,
                'facecolor': None,
                'clim': None,
                'logcolor': False,
                'colorlabel': None,
                'size_by': None,
                'size_by_log': None,
                'size_range': None}

    def __init__(self, channels, start, end, state=None, outdir='.',
                 tag=None, etg=None, **kwargs):
        super(TriggerSummaryPlot, self).__init__(channels, start, end,
                                                 state=state, outdir=outdir,
                                                 tag=tag, **kwargs)
        self.etg = etg

    @property
    def tag(self):
        """Unique identifier for this `TriggerSummaryPlot`.

        Extends the standard `TimeSeriesSummaryPlot` tag with the ETG
        and each of the column names.
        """
        try:
            return self._tag
        except AttributeError:
            tag = super(TriggerSummaryPlot, self).tag
            tag += '_%s' % re_cchar.sub('_', self.etg).upper()
            for column in ('x', 'y', 'color'):
                if self.pargs[column]:
                    tag += '_%s' % re_cchar.sub('_', self.pargs[column])
            return tag

    def process(self):
        # get columns
        xcolumn = self.pargs.pop('x', 'time')
        ycolumn = self.pargs.pop('y', 'snr')
        ccolumn = self.pargs.get('color', None)

        # initialise figure
        base = xcolumn == 'time' and TimeSeriesPlot or None
        plot = self.plot = EventTablePlot(base=base)
        ax = plot.gca()
        if isinstance(plot, TimeSeriesPlot):
            ax.set_epoch(self.start)
            ax.auto_gps_scale(self.end-self.start)
            ax.set_xlim(self.start, self.end)

        # get colouring params
        clim = self.pargs.pop('clim', None)
        clog = self.pargs.pop('logcolor', False)
        clabel = self.pargs.pop('colorlabel', None)

        # get plot arguments
        plotargs = dict()
        plotargs['edgecolor'] = self.pargs.pop('edgecolor')
        plotargs['facecolor'] = self.pargs.pop('facecolor')
        plotargs['size_by'] = self.pargs.pop('size_by')
        plotargs['size_by_log'] = self.pargs.pop('size_by_log')
        plotargs['size_range'] = self.pargs.pop('size_range')
        if (plotargs['size_range'] is None and
                (plotargs['size_by'] or plotargs['size_by_log']) and
                clim is not None):
            size_range = clim

        # add data
        ntrigs = 0
        for channel in self.channels:
            channel = get_channel(channel)
            table = get_triggers(str(channel), self.etg, self.state,
                                 query=False)
            ntrigs += len(table)
            # access channel parameters for limits
            for c, column in zip(('x', 'y', 'c'), (xcolumn, ycolumn, ccolumn)):
                if not column:
                    continue
                # hack for SnglBurst frequency nonsense
                if column in ['peak_frequency', 'central_freq']:
                    column = 'frequency'
                # set x and y in plotargs
                param = '%s_range' % column
                lim = '%slim' % c
                if hasattr(channel, param) and c in ('x', 'y'):
                    self.pargs.setdefault(lim, getattr(channel, param))
                # set clim separately
                elif hasattr(channel, param):
                    if not clim:
                        clim = getattr(channel, param)
                    if not plotargs['size_range']:
                        plotargs['size_range'] = getattr(channel, param)

            ax.plot_table(table, xcolumn, ycolumn, color=ccolumn, **plotargs)

        # customise plot
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if 'title' not in self.pargs.keys() and len(self.channels) == 1:
            plot.title = '%s (%s)' % (self.channels[0].texname, self.etg)

        # add colorbar
        if ccolumn:
            if not ntrigs:
                ax.scatter([1], [1], c=[1], visible=False, **plotargs)
            plot.add_colorbar(ax=ax, clim=clim, log=clog, label=clabel)
        else:
            plot.add_colorbar(ax=ax, visible=False)

        # add state segments
        if isinstance(plot, TimeSeriesPlot):
            ax = self.add_state_segments(plot, ax)

        # finalise
        self.finalize()


for PlotClass in __all__:
    PlotClass = locals()[PlotClass]
    register_plot(PlotClass.type, PlotClass)
