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

from gwpy.detector import Channel
from gwpy.timeseries import StateVector
from gwpy.plotter import *

from .. import version
from ..utils import re_cchar
from ..data import (get_channel, get_timeseries, get_spectrogram, get_spectrum)
from ..segments import get_segments
from ..triggers import get_triggers
from ..state import ALLSTATE
from .core import TabPlot
from .registry import register_plot

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

__all__ = ['TimeSeriesTabPlot', 'SegmentTabPlot', 'SpectrumTabPlot',
           'SpectrogramPlot', 'StateVectorTabPlot', 'TriggerTabPlot']


class TimeSeriesTabPlot(TabPlot):
    """Request for a time-series plot
    """
    FigureClass = TimeSeriesPlot
    AxesClass = TimeSeriesAxes
    type = 'timeseries'

    @property
    def tag(self):
        try:
            return self._tag
        except AttributeError:
            self._tag = hashlib.md5("".join(map(str,
                                            self.channels))).hexdigest()[:6]
            return self.tag

    def add_data_source(self, source):
        self._channels.append(source)

    def add_state_segments(self, plot, ax, color='green', **kwargs):
        sax = plot.add_state_segments(self.state, ax,
                                     plotargs={'edgecolor': 'black',
                                               'facecolor': 'green'})
        sax.tick_params(axis='y', which='major', labelsize=12)
        sax.set_epoch(self.gpsstart)
        sax.auto_gps_scale(self.gpsend - self.gpsstart)
        return sax

    def process(self):
        plot = self.FigureClass()
        ax = plot._add_new_axes(self.AxesClass.name)
        # add data
        mmmchans = self.find_mmm_channels()
        labels = self.plotargs.pop('labels', mmmchans.keys())
        if isinstance(labels, basestring):
            labels = labels.split(',')
        labels = map(lambda s: str(s).strip('\n '), labels)
        for label, channels in zip(labels, zip(*mmmchans.items())[1]):
            data = [get_timeseries(c, self.state,
                                   query=False).join(pad=numpy.nan) for
                    c in channels]
            if not 'x0' in data[0].metadata:
                data[0].epoch = self.gpsstart
            if 'logy' in self.plotargs and self.plotargs['logy']:
                for ts in data:
                    ts[ts.data == 0] = 1e-100
            if len(channels) > 1:
                ax.plot_timeseries_mmm(*data, label=label.replace('_', r'\_'))
            else:
                ax.plot_timeseries(data[0], label=label.replace('_', r'\_'))

            # allow channel data to set parameters
            if hasattr(data[0].channel, 'amplitude_range'):
                self.plotargs.setdefault('ylim',
                                         data[0].channel.amplitude_range)

        ax.set_epoch(self.gpsstart)
        # find hlines to add
        hlines = self.plotargs.pop('hline', [])
        for yval in hlines:
            try:
                yval = float(yval)
            except ValueError:
                continue
            else:
                ax.plot([self.gpsstart, self.gpsend], [yval, yval],
                        linestyle='--', color='red')
        # customise plot
        for key, val in self.plotargs.iteritems():
            try:
                setattr(plot, key, val)
            except AttributeError:
                getattr(plot, 'get_%s' % key)(val)
        if len(mmmchans) > 1:
            plot.add_legend(ax=ax)
        if not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)
        segax = self.add_state_segments(plot, ax)
        if 'xlim' not in self.plotargs.keys():
            segax.set_xlim(self.gpsstart, self.gpsend)
        plot.save(self.outputfile)
        plot.close()


class SpectrumTabPlot(TabPlot):
    """Spectrum plot for a `SummaryTab`
    """
    FigureClass = SpectrumPlot
    AxesClass = SpectrumAxes
    type = 'spectrum'
    defaults = {'logx': True, 'logy': True}

    @property
    def tag(self):
        try:
            return self._tag
        except AttributeError:
            self._tag = hashlib.md5("".join(map(str,
                                            self.channels))).hexdigest()[:6]
            return self.tag

    def add_data_source(self, source):
        self._channels.append(source)

    def process(self):
        """Load all data, and generate this `SpectrumTabPlot`
        """
        # get spectrum format: 'amplitude' or 'power'
        if 'format' in self.plotargs:
            sdform = self.plotargs['format']
        else:
            sdform = None
        plot = self.FigureClass(figsize=[12, 6])
        ax = plot._add_new_axes(self.AxesClass.name)
        # add data
        labels = self.plotargs.pop('labels', map(str, self.channels))
        if isinstance(labels, basestring):
            labels = labels.split(',')
        labels = map(lambda s: str(s).strip('\n '), labels)
        colors = self.plotargs.pop('color', [])
        if isinstance(colors, str):
            colors = colors.split(',')
        while len(colors) < len(self.channels):
            colors.append(None)
        for label, color, channel in zip(labels, colors, self.channels):
            data = get_spectrum(str(channel), self.state, query=False,
                                format=sdform)

            if 'logx' in self.plotargs and self.plotargs['logx']:
                data = [s[1:] for s in data]
            if 'logy' in self.plotargs and self.plotargs['logy']:
                for sp in data:
                    sp[sp.data == 0] = 1e-100
            if color is not None:
                ax.plot_spectrum_mmm(*data, label=label.replace('_', r'\_'),
                                     color=color)
            else:
                ax.plot_spectrum_mmm(*data, label=label.replace('_', r'\_'))

            # allow channel data to set parameters
            if hasattr(data[0].channel, 'frequency_range'):
                self.plotargs.setdefault('xlim',
                                         data[0].channel.frequency_range)
            if (sdform in ['amplitude', 'asd'] and
                hasattr(data[0].channel, 'asd_range')):
                self.plotargs.setdefault('ylim', data[0].channel.asd_range)
            elif hasattr(data[0].channel, 'psd_range'):
                self.plotargs.setdefault('ylim', data[0].channel.psd_range)

        # set plotting params
        legendargs = {}
        for key, val in self.plotargs.iteritems():
            if re.match('legend-', key):
                legendargs[key[7:]] = val
            try:
                setattr(plot, key, val)
            except AttributeError:
                getattr(plot, 'get_%s' % key)(val)
        if len(self.channels) > 1:
            plot.add_legend(ax=ax, **legendargs)
        if not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)
        for ax in plot.axes:
            ax.tick_params(axis='x', pad=10)
            ax.xaxis.labelpad = 10
        plot.save(self.outputfile)
        plot.close()


class SegmentTabPlot(TabPlot):
    """Segment plot for a `SummaryTab`
    """
    FigureClass = SegmentPlot
    AxesClass = SegmentAxes
    type = 'segments'
    FACECOLOR = {}

    def __init__(self, flags, start, end, state='all', outdir='.', href=None,
                 **kwargs):
        self.flags = flags
        self.gpsstart = start
        self.gpsend = end
        self.state = state
        self.outdir = outdir
        self.plotargs = kwargs
        self.href = href

    @property
    def tag(self):
        try:
            return self._tag
        except AttributeError:
            self._tag = hashlib.md5("".join(map(str,
                                                self.flags))).hexdigest()[:6]
            return self.tag

    def add_data_source(self, source):
        self.flags.append(str(source))

    @property
    def ifos(self):
        """Interferometer set for this `SegmentTabPlot`
        """
        allflags = [f for flag in self.flags for f in flag.split(',')]
        return set([f[:2] for f in allflags])

    @classmethod
    def from_ini(cls, cp, section, state=ALLSTATE, **kwargs):
        return super(SegmentTabPlot, cls).from_ini(cp, section, state=state,
                                                   source='data-quality-flags',
                                                   **kwargs)

    def process(self):
        # separate plot arguments
        mask = self.plotargs.pop('mask', None)
        axargs = self.plotargs.copy()
        flags = map(lambda f: str(f).replace('_', r'\_'), self.flags)
        labels = self.plotargs.pop('labels', self.plotargs.pop('label', flags))
        addlabel = self.plotargs.pop('add_label', None)
        if addlabel is None:
            inset = self.plotargs.pop('inset_labels', True)
            if inset is not False:
                addlabel = 'inset'
            else:
                addlabel = True
        if isinstance(labels, basestring):
            labels = labels.split(',')
        labels = map(lambda s: str(s).strip('\n '), labels)
        plotargs = {}
        color = axargs.pop('color', None)
        if (color is None and len(self.ifos) == 1 and
            list(self.ifos)[0] in self.FACECOLOR):
            color = self.FACECOLOR[list(self.ifos)[0]]
        elif color is None:
            color = 'green'
        plotargs['edgecolor'] = axargs.pop('edgecolor', None)
        plotargs['facecolor'] = color
        plotargs['add_label'] = addlabel

        # make plot
        plot = self.FigureClass()
        ax = plot._add_new_axes(self.AxesClass.name)
        for flag, label in zip(self.flags, labels)[::-1]:
            segs = get_segments(flag, validity=self.state.active, query=False)
            ax.plot(segs, label=label, valid={'facecolor': 'red'}, **plotargs)
        ax.set_epoch(self.gpsstart)
        ax.auto_gps_scale(self.gpsend-self.gpsstart)
        for key, val in axargs.iteritems():
            try:
                setattr(plot, key, val)
            except AttributeError:
                getattr(plot, 'get_%s' % key)(val)
        if 'xlim' not in axargs:
            ax.set_xlim(self.gpsstart, self.gpsend)
        if mask is None and not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)
        if 'ylim' not in axargs:
            ax.set_ylim(-0.5, len(self.flags) - 0.5)
        if mask is not None:
            plot.add_bitmask(mask, topdown=True)
        plot.save(self.outputfile)
        plot.close()


class SpectrogramPlot(TimeSeriesTabPlot):
    """Plot a Spectrogram
    """
    type = 'spectrogram'

    @property
    def tag(self):
        try:
            return self._tag
        except AttributeError:
            self._tag = hashlib.md5("".join(map(str,
                                            self.channels))).hexdigest()[:6]
            if 'ratio' in self.plotargs:
                ratio = self.plotargs['ratio']
                self._tag += '_%s_RATIO' % re_cchar.sub('_', ratio.upper())
            return self.tag

    def process(self):
        # parse arguments
        sdform = self.plotargs.pop('format', None)
        clim = self.plotargs.pop('clim', None)
        clog = self.plotargs.pop('logcolor', False)
        clabel = self.plotargs.pop('colorlabel', None)
        ratio = self.plotargs.pop('ratio', None)

        # make figure
        plot = self.FigureClass()
        ax = plot._add_new_axes(self.AxesClass.name)
        # add data
        specgrams = get_spectrogram(self.channels[0], self.state, query=False,
                                    format=sdform)
        if ratio in ['median', 'mean']:
            allspec = specgrams.join(gap='ignore')
            ratio = getattr(allspec, ratio)(axis=0)
        for specgram in specgrams:
            if ratio is not None:
                specgram = specgram.ratio(ratio)
            ax.plot_spectrogram(specgram)

            # allow channel data to set parameters
            if hasattr(specgram.channel, 'frequency_range'):
                self.plotargs.setdefault('ylim',
                                         specgram.channel.frequency_range)
            if ratio is None and (sdform in ['amplitude', 'asd'] and
                    hasattr(specgram.channel, 'asd_range') and clim is None):
                clim = specgram.channel.asd_range
            elif ratio is None and hasattr(specgram.channel, 'psd_range') and clim is None:
                clim = specgram.channel.psd_range
        if len(specgrams) == 0:
            ax.scatter([1], [1], c=[1], visible=False)
        ax.set_epoch(self.gpsstart)
        ax.auto_gps_scale(self.gpsend-self.gpsstart)
        # add colorbar
        plot.add_colorbar(ax=ax, clim=clim, log=clog, label=clabel)
        for key, val in self.plotargs.iteritems():
            try:
                setattr(plot, key, val)
            except AttributeError:
                getattr(plot, 'get_%s' % key)(val)
        segax = self.add_state_segments(plot, ax)
        if 'xlim' not in self.plotargs.keys():
            segax.set_xlim(self.gpsstart, self.gpsend)
        ax.grid(b=True, axis='y', which='major')
        plot.save(self.outputfile)
        plot.close()


class StateVectorTabPlot(TimeSeriesTabPlot):
    type = 'statevector'
    FigureClass = SegmentPlot
    AxesClass = SegmentAxes
    FACECOLOR = {}

    def process(self):
        # separate plot arguments
        mask = self.plotargs.pop('mask', None)
        axargs = self.plotargs.copy()
        plotargs = {}
        color = axargs.pop('color', None)
        if (color is None and len(self.ifos) == 1 and
            list(self.ifos)[0] in self.FACECOLOR):
            color = self.FACECOLOR[list(self.ifos)[0]]
        elif color is None:
            color = 'green'
        plotargs['edgecolor'] = axargs.pop('edgecolor', None)
        plotargs['facecolor'] = color
        plotargs.setdefault('add_label', 'inset')

        # make figure
        plot = self.FigureClass()
        ax = plot._add_new_axes(self.AxesClass.name)

        # add data
        labels = self.plotargs.pop('labels', self.channels)
        labels = map(lambda s: str(s).strip('\n '), labels)
        nflags = 0
        for label, channel in zip(labels, self.channels[::-1]):
            data = get_timeseries(str(channel), self.state, query=False,
                                  statevector=True)
            flags = None
            for stateseries in data:
                if not stateseries.size:
                    stateseries.epoch = self.gpsstart
                    stateseries.dx = 0
                    if channel.sample_rate is not None:
                        stateseries.sample_rate = channel.sample_rate
                stateseries.bitmask = channel.bitmask
                newflags = stateseries.to_dqflags()
                if flags is None:
                    flags = newflags
                else:
                    for i,flag in enumerate(newflags):
                        flags[i] += flag
            nflags += len([m for m in channel.bitmask if m is not None])
            ax.plot(*flags[::-1], valid={'facecolor': 'red'}, **plotargs)
        ax.set_epoch(self.gpsstart)
        ax.auto_gps_scale(self.gpsend-self.gpsstart)
        for key, val in self.plotargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                getattr(ax, 'set_%s' % key)(val)
        if 'xlim' not in self.plotargs.keys():
            ax.set_xlim(self.gpsstart, self.gpsend)
        if 'ylim' not in self.plotargs.keys():
            ax.set_ylim(-0.5, nflags-0.5)
        if mask is None and not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)
        if mask is not None:
            plot.add_bitmask(mask, topdown=True)
        plot.save(self.outputfile)
        plot.close()


class TriggerTabPlot(TimeSeriesTabPlot):
    type = 'triggers'
    FigureClass = EventTablePlot
    AxesClass = EventTableAxes

    def __init__(self, channels, start, end, state=ALLSTATE, outdir='.',
                 href=None, etg=None, **kwargs):
        super(TriggerTabPlot, self).__init__(channels, start, end, state=state,
                                             outdir=outdir, href=href, **kwargs)
        if etg is None:
            raise ValueError("An 'etg' option must be given in the INI "
                             "section defining this plot")
        self.etg = etg

    @property
    def tag(self):
        """Unique identifier for this `TriggerTabPlot`.

        Extends the standard `TimeSeriesTabPlot` tag with the ETG
        and each of the column names.
        """
        try:
            return self._tag
        except AttributeError:
            self._tag = hashlib.md5("".join(map(str,
                                            self.channels))).hexdigest()[:6]
            self._tag += '_%s' % re_cchar.sub('_', self.etg).upper()
            for column in ('x', 'y', 'color'):
                if column in self.plotargs:
                    self._tag += '_%s' % re_cchar.sub('_',
                                                      self.plotargs[column])
            return self.tag

    def process(self):
        # get plot arguments
        plotargs = dict()
        xcolumn = self.plotargs.pop('x', 'time')
        ycolumn = self.plotargs.pop('y', 'snr')
        ccolumn = self.plotargs.get('color', None)
        plotargs['edgecolor'] = self.plotargs.pop('edgecolor', None)
        plotargs['facecolor'] = self.plotargs.pop('facecolor', None)
        # get colouring params
        clim = self.plotargs.pop('clim', None)
        clog = self.plotargs.pop('logcolor', False)
        clabel = self.plotargs.pop('colorlabel', None)
        size_by = self.plotargs.pop('size_by', None)
        size_by_log = self.plotargs.pop('size_by_log', None)
        size_range = self.plotargs.pop('size_range', None)
        if size_range is None and (size_by or size_by_log) and clim is not None:
            size_range = clim
        # generate figure
        base = xcolumn == 'time' and TimeSeriesPlot or None
        plot = self.FigureClass(base=base)
        ax = plot._add_new_axes(self.AxesClass.name)
        # add data
        ntrigs = 0
        for channel in self.channels:
            channel = get_channel(channel)
            table = get_triggers(str(channel), self.etg, self.state,
                                 query=False)
            ntrigs += len(table)
            # access channel parameters for limits
            for c,column in zip(('x', 'y', 'c'), (xcolumn, ycolumn, ccolumn)):
                if not column:
                    continue
                # hack for SnglBurst frequency nonsense
                if column in ['peak_frequency', 'central_freq']:
                    column = 'frequency'
                # set x and y in plotargs
                param = '%s_range' % column
                lim = '%slim' % c
                if hasattr(channel, param) and c in ('x', 'y'):
                    self.plotargs.setdefault(lim, getattr(channel, param))
                # set clim separately
                elif hasattr(channel, param):
                    if not clim:
                        clim = getattr(channel, param)
                    if not size_range:
                        size_range = getattr(channel, param)

            ax.plot_table(table, xcolumn, ycolumn, color=ccolumn,
                          size_by=size_by, size_by_log=size_by_log,
                          size_range=size_range, **plotargs)

        ax.set_epoch(self.gpsstart)
        ax.auto_gps_scale(self.gpsend-self.gpsstart)
        for key, val in self.plotargs.iteritems():
            try:
                setattr(plot, key, val)
            except AttributeError:
                getattr(plot, 'get_%s' % key)(val)
        if 'title' not in self.plotargs.keys() and len(self.channels) == 1:
            plot.title = '%s (%s)' % (self.channels[0].texname, self.etg)
        if ccolumn:
            if not ntrigs:
                ax.scatter([1], [1], c=[1], visible=False, **plotargs)
            plot.add_colorbar(ax=ax, clim=clim, log=clog, label=clabel)
        else:
            plot.add_colorbar(ax=ax, visible=False)
        if isinstance(plot, TimeSeriesPlot):
            ax = self.add_state_segments(plot, ax)
        if 'xlim' not in self.plotargs.keys():
            ax.set_xlim(self.gpsstart, self.gpsend)
        plot.save(self.outputfile)
        plot.close()


for PlotClass in __all__:
    PlotClass = locals()[PlotClass]
    register_plot(PlotClass.type, PlotClass)
