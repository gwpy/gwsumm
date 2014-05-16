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
from itertools import (izip, cycle)

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from gwpy.plotter import *
from gwpy.plotter.table import get_column_string
from gwpy.plotter.utils import (color_cycle, marker_cycle)
from gwpy.table.rate import (event_rate, binned_event_rates)
from gwpy.table.utils import get_table_column

from .. import (globalv, version)
from ..utils import (re_quote, re_cchar, split_channels)
from ..data import (get_channel, get_timeseries, get_spectrogram, get_spectrum,
                    add_timeseries)
from ..segments import get_segments
from ..triggers import get_triggers
from ..state import ALLSTATE
from .core import DataSummaryPlot
from .registry import register_plot

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


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
        kwargs.setdefault('valid', {'facecolor': 'red'})
        sax = self.plot.add_state_segments(self.state, ax, plotargs=kwargs)
        sax.tick_params(axis='y', which='major', labelsize=12)
        sax.set_epoch(self.start)
        return sax

    def init_plot(self, plot=TimeSeriesPlot):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesSummaryPlot`.
        """
        self.plot = plot()
        ax = self.plot.gca()
        ax.set_epoch(self.start)
        return self.plot, ax

    def parse_kwargs(self):
        plotargs = {}
        for kwarg in ['alpha', 'color', 'drawstyle', 'fillstyle', 'linestyle',
                       'linewidth', 'marker', 'markeredgecolor',
                       'markeredgewidth', 'markerfacecolor',
                       'markerfacecoloralt', 'markersize']:
            try:
                plotargs[kwarg] = self.pargs[kwarg]
            except KeyError:
                pass
        return plotargs

    def finalize(self, outputfile=None):
        plot = self.plot
        ax = plot.axes[0]
        if 'xlim' not in self.pargs:
            ax.set_xlim(self.start, self.end)
        return super(TimeSeriesSummaryPlot, self).finalize(
                   outputfile=outputfile)

    def process(self, outputfile=None):
        """Read in all necessary data, and generate the figure.
        """
        (plot, ax) = self.init_plot()

        plotargs = self.parse_kwargs()

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
                ax.plot_timeseries_mmm(*data, label=label.replace('_', r'\_'),
                                       **plotargs)
            else:
                ax.plot_timeseries(data[0], label=label.replace('_', r'\_'),
                                   **plotargs)

            # allow channel data to set parameters
            if hasattr(data[0].channel, 'amplitude_range'):
                self.pargs.setdefault('ylim',
                                      data[0].channel.amplitude_range)

        # add horizontal lines to add
        for yval in self.pargs['hline']:
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
            plot.add_legend(ax=ax, markerscale=4)

        # add extra axes and finalise
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        self.add_state_segments(ax)
        return self.finalize(outputfile=outputfile)

register_plot(TimeSeriesSummaryPlot)


class SpectrogramSummaryPlot(TimeSeriesSummaryPlot):
    """SummaryPlot a Spectrogram
    """
    type = 'spectrogram'
    defaults = {'ratio': None,
                'format': None,
                'clim': None,
                'cmap': 'jet',
                'logcolor': False,
                'colorlabel': None}

    def __init__(self, *args, **kwargs):
        super(SpectrogramSummaryPlot, self).__init__(*args, **kwargs)
        self.ratio = self.pargs.pop('ratio')

    @property
    def tag(self):
        try:
            return self._tag
        except AttributeError:
            tag = super(SpectrogramSummaryPlot, self).tag
            if self.ratio:
                tag += '_%s_RATIO' % re_cchar.sub('_', self.ratio.upper())
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
        cmap = self.pargs.pop('cmap')
        ratio = self.ratio

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
            ax.plot_spectrogram(specgram, cmap=cmap)

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
            ax.scatter([1], [1], c=[1], visible=False, cmap=cmap)
        plot.add_colorbar(ax=ax, clim=clim, log=clog, label=clabel, cmap=cmap)

        # customise and finalise
        for key, val in self.pargs.iteritems():
            if key == 'ratio':
                continue
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        self.add_state_segments(ax)
        return self.finalize()

register_plot(SpectrogramSummaryPlot)


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
            return '_'.join([state, hash_, self.type]).upper()

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
        labels = map(lambda s: re_quote.sub('', str(s).strip('\n ')), labels)

        # extract plotting arguments
        mask = self.pargs.pop('mask')
        activecolor, validcolor = self.get_segment_color()
        edgecolor = self.pargs.pop('edgecolor')
        plotargs = {'add_label': addlabel,
                    'facecolor': activecolor,
                    'edgecolor': edgecolor,
                    'valid': {'facecolor': validcolor}}
        for key in plotargs:
            if key in self.pargs:
                plotargs[key] = self.pargs.pop(key)

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
        if mask is None and not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        elif mask is not None:
            plot.add_bitmask(mask, topdown=True)
        return self.finalize()

register_plot(SegmentSummaryPlot)


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
            if flags is not None:
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
        if mask is None and not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        elif mask is not None:
            plot.add_bitmask(mask, topdown=True)
        return self.finalize()

register_plot(StateVectorSummaryPlot)


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
        plot = self.plot = SpectrumPlot(figsize=[12, 6])
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
            if re.match('legend[-_]', key):
                legendargs[key[7:]] = val
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if len(self.channels) > 1:
            plot.add_legend(ax=ax, **legendargs)
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)

        return self.finalize()

register_plot(SpectrumSummaryPlot)


class TimeSeriesHistogramPlot(DataSummaryPlot):
    """HistogramPlot from a Series
    """
    type = 'histogram'
    defaults = {'ylabel': 'Rate [Hz]',
                'rwidth': 1}

    def init_plot(self, plot=HistogramPlot):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesSummaryPlot`.
        """
        self.plot = plot(figsize=[12, 6])
        ax = self.plot.gca()
        return self.plot, ax

    def parse_histogram_kwargs(self):
        """Extract the histogram arguments from everything else.
        """
        histargs = {}
        for param in ['bins', 'range', 'normed', 'weights', 'cumulative',
                      'bottom', 'histtype', 'align', 'orientation', 'rwidth',
                      'log', 'color', 'label', 'stacked', 'logbins',
                      'alpha', 'linecolor', 'edgecolor', 'facecolor']:
            try:
                histargs[param] = self.pargs.pop(param)
            except KeyError:
                if param == 'range':
                    try:
                        histargs[param] = self.pargs['xlim']
                    except KeyError:
                        pass
                if param == 'log':
                    histargs['log'] = self.pargs.pop('logy', False)
            else:
                if isinstance(histargs[param], (unicode, str)):
                    try:
                        histargs[param] = eval(histargs[param])
                    except:
                        pass
        # remove logy if already requesting log histogram
        if histargs.get('log', True) and self.pargs.get('logy', True):
            self.pargs.pop('logy', None)
        # set alpha
        if len(self.channels) > 1:
             histargs.setdefault('alpha', 0.7)
        return histargs

    def process(self):
        """Get data and generate the figure.
        """
        # get histogram parameters
        (plot, ax) = self.init_plot()

        # extract histogram arguments
        histargs = self.parse_histogram_kwargs()

        # work out labels
        labels = self.pargs.pop('labels', self.channels)
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        labels = map(lambda s: str(s).strip('\n '), labels)

        # get data
        data = []
        for channel in self.channels:
            data.append(get_timeseries(channel, self.state,
                                       query=False).join(pad=numpy.nan))
            # allow channel data to set parameters
            if hasattr(data[-1].channel, 'amplitude_range'):
                self.pargs.setdefault('xlim', data[-1].channel.amplitude_range)

        # get range
        if not 'range' in histargs:
            histargs['range'] = ax.common_limits(data)

        # plot
        for label, arr in zip(labels, data):
            if arr.size:
                ax.hist(arr, label=label, **histargs)
            else:
                ax.plot([], label=label)

        # customise plot
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if len(self.channels) > 1:
            plot.add_legend(ax=ax)

        # add extra axes and finalise
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        return self.finalize()

register_plot(TimeSeriesHistogramPlot)


class TriggerSummaryPlot(TimeSeriesSummaryPlot):
    type = 'triggers'
    defaults = {'x': 'time',
                'y': 'snr',
                'color': None,
                'edgecolor': 'face',
                'facecolor': None,
                'marker': 'o',
                's': 20,
                'vmin': None,
                'vmax': None,
                'clim': None,
                'logcolor': False,
                'colorlabel': None,
                'cmap': 'jet',
                'size_by': None,
                'size_by_log': None,
                'size_range': None}

    def __init__(self, channels, start, end, state=None, outdir='.',
                 tag=None, etg=None, **kwargs):
        super(TriggerSummaryPlot, self).__init__(channels, start, end,
                                                 state=state, outdir=outdir,
                                                 tag=tag, **kwargs)
        self.etg = etg
        self.columns = [self.pargs.pop(c) for c in ('x', 'y', 'color')]

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
            tag += '_%s' % re_cchar.sub('_', self.etg)
            for column in self.columns:
                if column:
                    tag += '_%s' % re_cchar.sub('_', column)
            return tag.upper()

    def finalize(self, outputfile=None):
        if isinstance(self.plot, TimeSeriesPlot):
            return super(TriggerSummaryPlot, self).finalize(
                       outputfile=outputfile)
        else:
            return super(TimeSeriesSummaryPlot, self).finalize(
                       outputfile=outputfile)

    def process(self):
        # get columns
        xcolumn, ycolumn, ccolumn = self.columns

        # initialise figure
        if 'time' in xcolumn:
            base = TimeSeriesPlot
        elif 'freq' in xcolumn:
            base = SpectrumPlot
        else:
            base = Plot
        plot = self.plot = EventTablePlot(figsize=[12, 6], base=base)
        ax = plot.gca()
        if isinstance(plot, TimeSeriesPlot):
            ax.set_epoch(self.start)
            ax.set_xlim(self.start, self.end)

        # work out labels
        labels = self.pargs.pop('labels', self.channels)
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        labels = map(lambda s: str(s).strip('\n '), labels)

        # get colouring params
        clim = self.pargs.pop('clim', self.pargs.pop('colorlim', None))
        clog = self.pargs.pop('logcolor', False)
        clabel = self.pargs.pop('colorlabel', None)

        # get plot arguments
        plotargs = []
        for i in range(len(self.channels)):
            plotargs.append(dict())
        # get plot arguments
        for key in ['vmin', 'vmax', 'edgecolor', 'facecolor', 'size_by',
                    'size_by_log', 'size_range', 'cmap', 's', 'marker']:
            val = self.pargs.pop(key)
            if key == 'facecolor' and len(self.channels) > 1 and val is None:
                val = color_cycle()
            if key == 'marker' and len(self.channels) > 1 and val is None:
                val = marker_cycle()
            elif (isinstance(val, (list, tuple, cycle)) and
                  key not in ['size_range']):
                val = cycle(val)
            else:
                val = cycle([val] * len(self.channels))
            for i in range(len(self.channels)):
                plotargs[i][key] = val.next()
        # fix size_range
        if (len(self.channels) == 1 and plotargs[0]['size_range'] is None and
                (plotargs[0]['size_by'] or plotargs[0]['size_by_log']) and
                clim is not None):
            plotargs[0]['size_range'] = clim

        # add data
        ntrigs = 0
        for channel, label, pargs in izip(self.channels, labels, plotargs):
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
                    if not pargs['size_range']:
                        pargs['size_range'] = getattr(channel, param)

            ax.plot_table(table, xcolumn, ycolumn, color=ccolumn,
                          label=label, **pargs)

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
                ax.scatter([1], [1], c=[1], visible=False)
            plot.add_colorbar(ax=ax, clim=clim, log=clog, label=clabel)
        else:
            plot.add_colorbar(ax=ax, visible=False)

        if len(self.channels) == 1 and len(table):
            if ccolumn is None:
                ax.add_loudest(table, ycolumn, xcolumn, ycolumn)
            else:
                ax.add_loudest(table, ccolumn, xcolumn, ycolumn)

        if len(self.channels) > 1:
            plot.add_legend(ax=ax, markerscale=4)

        # add state segments
        if isinstance(plot, TimeSeriesPlot):
            self.add_state_segments(ax)

        # finalise
        return self.finalize()

register_plot(TriggerSummaryPlot)


class TriggerTimeSeriesSummaryPlot(TimeSeriesSummaryPlot):
    """Custom time-series plot to handle discontiguous `TimeSeries`.
    """
    type = 'trigger-timeseries'

    def process(self):
        """Read in all necessary data, and generate the figure.
        """
        (plot, ax) = self.init_plot()

        # work out labels
        labels = self.pargs.pop('labels', self.channels)
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        labels = map(lambda s: str(s).strip('\n '), labels)

        # add data
        for label, channel in zip(labels, self.channels):
            label = label.replace('_', r'\_')
            data = get_timeseries(channel, self.state, query=False)
            # handle no timeseries
            if not len(data):
                ax.plot([0], [0], visible=False, label=label)
                continue
            # plot time-series
            color = None
            for ts in data:
                # double-check log scales
                if self.pargs['logy']:
                    ts[ts.data == 0] = 1e-100
                if color is None:
                    line = ax.plot_timeseries(ts, label=label)[0]
                    color = line.get_color()
                else:
                    ax.plot_timeseries(ts, color=color, label=None)

            # allow channel data to set parameters
            if hasattr(data[0].channel, 'amplitude_range'):
                self.pargs.setdefault('ylim',
                                      data[0].channel.amplitude_range)

        # add horizontal lines to add
        for yval in self.pargs['hline']:
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
        if len(self.channels) > 1:
            plot.add_legend(ax=ax)

        # add extra axes and finalise
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        self.add_state_segments(ax)
        return self.finalize()

register_plot(TriggerTimeSeriesSummaryPlot)


class TriggerHistogramPlot(TimeSeriesHistogramPlot):
    """HistogramPlot from a LIGO_LW Table
    """
    type = 'trigger-histogram'

    def init_plot(self, plot=HistogramPlot):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesSummaryPlot`.
        """
        self.plot = plot(figsize=[12, 6])
        ax = self.plot.gca()
        return self.plot, ax

    def process(self):
        """Get data and generate the figure.
        """
        # get histogram parameters
        (plot, ax) = self.init_plot()

        etg = self.pargs.pop('etg')
        column = self.pargs.pop('column')

        # extract histogram arguments
        histargs = self.parse_histogram_kwargs()

        # work out labels
        labels = self.pargs.pop('labels', self.channels)
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        labels = map(lambda s: str(s).strip('\n '), labels)

        # add data
        data = []
        for label, channel in zip(labels, self.channels):
            channel = get_channel(channel)
            table_ = get_triggers(str(channel), etg, self.state,
                                  query=False)
            data.append(get_table_column(table_, column))
            # allow channel data to set parameters
            if hasattr(channel, 'amplitude_range'):
                self.pargs.setdefault('xlim', channel.amplitude_range)

        # get range
        if not 'range' in histargs:
            histargs['range'] = ax.common_limits(data)

        # plot
        for label, arr in zip(labels, data):
            if arr.size:
                ax.hist(arr, label=label, **histargs)
            else:
                ax.plot([], label=label)

        # customise plot
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if len(self.channels) > 1:
            plot.add_legend(ax=ax)

        # add extra axes and finalise
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        return self.finalize()

register_plot(TriggerHistogramPlot)


class TriggerRateSummaryPlot(TimeSeriesSummaryPlot):
    """TimeSeriesSummaryPlot of trigger rate.
    """
    type = 'trigger-rate'
    _threadsafe = False
    defaults = TimeSeriesSummaryPlot.defaults.copy()
    defaults.update({'column': None,
                     'ylabel': 'Rate [Hz]'})

    def __init__(self, *args, **kwargs):
        if not 'stride' in kwargs:
            raise ValueError("'stride' must be configured for all rate plots.")
        if 'column' in kwargs and 'bins' not in kwargs:
            raise ValueError("'bins' must be configured for rate plots if "
                             "'column' is given.")
        super(TriggerRateSummaryPlot, self).__init__(*args, **kwargs)

    def process(self):
        """Read in all necessary data, and generate the figure.
        """
        etg = self.pargs.pop('etg')

        # get rate arguments
        stride = self.pargs.pop('stride')
        column = self.pargs.pop('column', None)
        if column:
            cname = get_column_string(column)
            bins = self.pargs.pop('bins')
            operator = self.pargs.pop('operator', '>=')
        else:
            bins = ['_']

        # work out labels
        labels = self.pargs.pop('labels', None)
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        elif labels is None and column and len(self.channels) > 1:
            labels = []
            for channel, bin in [(c, b) for c in self.channels for b in bins]:
                labels.append(r' '.join([channel, cname, '$%s$' % str(operator),
                                         str(b)]))
        elif labels is None and column:
            labels = [r' '.join([cname,'$%s$' % str(operator), str(b)])
                      for b in bins]
        elif labels is None:
            labels = self.channels
        self.pargs['labels'] = map(lambda s: str(s).strip('\n '), labels)

        # generate data
        keys = []
        for channel in self.channels:
            table_ = get_triggers(str(channel), etg, self.state, query=False)
            if column:
                rates = binned_event_rates(
                    table_, stride, column, bins, operator, self.start,
                    self.end).values()
            else:
                rates = [event_rate(table_, stride, self.start, self.end)]
            for bin, rate in zip(bins, rates):
                keys.append('%s:RATE %s %s' % (channel.ifo, str(channel), bin))
                if keys[-1] not in globalv.DATA:
                    add_timeseries(rate, keys[-1])

        # reset channel lists and generate time-series plot
        channels = self.channels
        outputfile = self.outputfile
        self.channels = keys
        out = super(TriggerRateSummaryPlot, self).process(outputfile=outputfile)
        self.channels = channels
        return out

register_plot(TriggerRateSummaryPlot)
