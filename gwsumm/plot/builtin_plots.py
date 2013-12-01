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

import sys
import hashlib

from gwpy.detector import Channel
from gwpy.plotter import *

from .. import version
from ..utils import re_cchar
from ..data import (get_timeseries, get_spectrogram, get_spectrum)
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
        if isinstance(source, basestring):
            source = Channel(source)
        self.channels.append(source)

    def process(self):
        plot = self.FigureClass()
        ax = plot._add_new_axes(self.AxesClass.name)
        # add data
        mmmchans = self.find_mmm_channels()
        labels = self.plotargs.pop('labels', mmmchans.keys())
        for label, channels in zip(labels, zip(*mmmchans.items())[1]):
            data = [get_timeseries(str(c), self.state, query=False).join() for c in channels]
            if len(channels) > 1:
                ax.plot_timeseries_mmm(*data, label=label.replace('_', r'\_'))
            else:
                ax.plot_timeseries(data[0])
        ax.set_epoch(self.gpsstart)
        for key, val in self.plotargs.iteritems():
            try:
                setattr(plot, key, val)
            except AttributeError:
                getattr(plot, 'get_%s' % key)(val)
        if 'xlim' not in self.plotargs.keys():
            ax.set_xlim(self.gpsstart, self.gpsend)
        if len(mmmchans) > 1:
            plot.add_legend(ax=ax)
        if not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)
        plot.add_state_segments(self.state, plotargs={'color':'green'})
        plot.subplots_adjust(left=0.1, right=0.9)
        plot.save(self.outputfile)


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
        if isinstance(source, basestring):
            source = Channel(source)
        self.channels.append(source)

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
        for label, channel in zip(labels, self.channels):
            data = get_spectrum(str(channel), self.state, query=False,
                                format=sdform)
            if 'logx' in self.plotargs and self.plotargs['logx']:
                data = [s[1:] for s in data]
            if 'logy' in self.plotargs and self.plotargs['logy']:
                for sp in data:
                    sp[sp.data == 0] = 1e-100
            ax.plot_spectrum_mmm(*data, label=label.replace('_', r'\_'))
        for key, val in self.plotargs.iteritems():
            try:
                setattr(plot, key, val)
            except AttributeError:
                getattr(plot, 'get_%s' % key)(val)
        if len(self.channels) > 1:
            plot.add_legend(ax=ax)
        if not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)
        plot.subplots_adjust(left=0.1, right=0.9)
        if sys.platform.startswith('linux'):
            for ax in plot.axes:
                ax.tick_params(axis='x', pad=10)
                ax.xaxis.labelpad = 10
        plot.save(self.outputfile)


class SegmentTabPlot(TabPlot):
    """Segment plot for a `SummaryTab`
    """
    FigureClass = SegmentPlot
    AxesClass = SegmentAxes
    type = 'segments'

    def __init__(self, flags, state='all', outdir='.', href=None, **kwargs):
        self.flags = flags
        self.state = state
        self.outdir = outdir
        self.plotargs = kwargs
        self.href = href

    @property
    def tag(self):
        try:
            return self._tag
        except AttributeError:
            self._tag = hashlib.md5("".join(map(str, self.flags))).hexdigest()[:6]
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
        axargs = self.plotargs.copy()
        plotargs = {}
        plotargs['color'] = axargs.pop('color', None)
        plotargs['edgecolor'] = axargs.pop('edgecolor', None)
        plotargs['facecolor'] = axargs.pop('facecolor', None)
        plotargs.setdefault('add_label', 'inset')

        # make plot
        plot = self.FigureClass()
        ax = plot._add_new_axes(self.AxesClass.name)
        for flag in self.flags:
            segs = get_segments(flag, validity=self.state.active, query=False)
            ax.plot(segs, **plotargs)
        ax.set_epoch(self.gpsstart)
        for key, val in axargs.iteritems():
            try:
                setattr(plot, key, val)
            except AttributeError:
                getattr(plot, 'get_%s' % key)(val)
        if 'xlim' not in axargs:
            ax.set_xlim(self.gpsstart, self.gpsend)
        if not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)
        ax.autoscale(axis='y', tight=True)
        ylim = ax.get_ylim()
        pad = (ylim[1] - ylim[0]) * 0.05
        ax.set_ylim(ylim[0] - pad, ylim[1] + pad)
        plot.subplots_adjust(left=0.1, right=0.9)
        plot.save(self.outputfile)


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
        ax.auto_gps_scale()
        ax.set_epoch(self.gpsstart)
        # add colorbar
        plot.add_colorbar(ax=ax, clim=clim, log=clog, label=clabel)
        for key, val in self.plotargs.iteritems():
            try:
                setattr(plot, key, val)
            except AttributeError:
                getattr(plot, 'get_%s' % key)(val)
        plot.add_state_segments(self.state, plotargs={'edgecolor': 'black',
                                                      'facecolor': 'green'})
        if 'xlim' not in self.plotargs.keys():
            ax.set_xlim(self.gpsstart, self.gpsend)
        plot.subplots_adjust(left=0.1, right=0.9)
        plot.save(self.outputfile)


class StateVectorTabPlot(TimeSeriesTabPlot):
    type = 'statevector'
    FigureClass = SegmentAxes
    AxesClass = SegmentPlot

    def process(self):
        # separate plot arguments
        axargs = self.plotargs.copy()
        plotargs = {}
        plotargs['color'] = axargs.pop('color', None)
        plotargs['edgecolor'] = axargs.pop('edgecolor', None)
        plotargs['facecolor'] = axargs.pop('facecolor', None)
        plotargs.setdefault('add_label', 'inset')

        # make figure
        plot = self.FigureClass()
        ax = plot._add_new_axes(self.AxesClass.name)

        # add data
        labels = self.plotargs.pop('labels', self.channels)
        for label, channels in zip(labels, self.channels):
            data = [get_timeseries(str(c), self.state, query=False).join() for c in channels]
            data.bitmask = c.bitmask
            ax.plot(*data.to_dqflags(), **plotargs)
        ax.set_epoch(self.gpsstart)
        for key, val in self.plotargs.iteritems():
            try:
                setattr(plot, key, val)
            except AttributeError:
                getattr(plot, 'get_%s' % key)(val)
        if 'xlim' not in self.plotargs.keys():
            ax.set_xlim(self.gpsstart, self.gpsend)
        if not plot.coloraxes:
            plot.add_colorbar(ax=ax, visible=False)
        plot.subplots_adjust(left=0.1, right=0.9)
        plot.save(self.outputfile)
        plot.close()


class TriggerTabPlot(TimeSeriesTabPlot):
    type = 'triggers'
    FigureClass = EventTablePlot
    AxesClass = EventTableAxes

    def __init__(self, channels, state=ALLSTATE, outdir='.', href=None,
                 etg=None, **kwargs):
        super(TriggerTabPlot, self).__init__(channels, state=state,
                                             outdir=outdir, href=href, **kwargs)
        if etg is None:
            raise ValueError("An 'etg' option must be given in the INI "
                             "section defining this plot")
        self.etg = etg

    @property
    def tag(self):
        try:
            return self._tag
        except AttributeError:
            self._tag = hashlib.md5("".join(map(str,
                                            self.channels))).hexdigest()[:6]
            self._tag += '_%s' % re_cchar.sub('_', self.etg).upper()
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
        size_range = self.plotargs.pop('size_range', clim)
        # generate figure
        base = xcolumn == 'time' and TimeSeriesPlot or None
        plot = self.FigureClass(base=base)
        ax = plot._add_new_axes(self.AxesClass.name)
        # add data
        ntrigs = 0
        for channel in self.channels:
            table = get_triggers(str(channel), self.etg, self.state,
                                 query=False)
            ntrigs += len(table)
            ax.plot_table(table, xcolumn, ycolumn, color=ccolumn,
                          size_by=size_by, size_by_log=size_by_log,
                          size_range=size_range, **plotargs)
        ax.set_epoch(self.gpsstart)
        ax.auto_gps_scale()
        for key, val in self.plotargs.iteritems():
            try:
                setattr(plot, key, val)
            except AttributeError:
                getattr(plot, 'get_%s' % key)(val)
        if 'title' not in self.plotargs.keys() and len(self.channels) == 1:
            plot.title = '%s (%s)' % (self.channels[0].tex_name, self.etg)
        if 'xlim' not in self.plotargs.keys():
            ax.set_xlim(self.gpsstart, self.gpsend)
        if ccolumn:
            if not ntrigs:
                ax.scatter([1], [1], c=[1], visible=False, **plotargs)
            plot.add_colorbar(ax=ax, clim=clim, log=clog, label=clabel)
        else:
            plot.add_colorbar(ax=ax, visible=False)
        if isinstance(plot, TimeSeriesPlot):
            plot.add_state_segments(self.state, ax=ax,
                                    plotargs={'color':'green'})
        plot.subplots_adjust(left=0.1, right=0.9)
        plot.save(self.outputfile)
        plot.close()


for PlotClass in __all__:
    PlotClass = locals()[PlotClass]
    register_plot(PlotClass.type, PlotClass)
