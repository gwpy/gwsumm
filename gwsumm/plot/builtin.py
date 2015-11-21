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
import warnings
from itertools import cycle

from matplotlib.pyplot import subplots

from astropy.units import Quantity

from gwpy.spectrum import Spectrum
from gwpy.plotter import *
from gwpy.plotter.tex import label_to_latex

from .. import (globalv, mode, version)
from ..utils import (re_quote, re_cchar, split_channels)
from ..data import (get_channel, get_timeseries, get_spectrogram, get_spectrum,
                    add_timeseries)
from ..state import ALLSTATE
from .registry import (get_plot, register_plot)
from .mixins import *

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

DataPlot = get_plot('data')
GREEN = (0.2, 0.8, 0.2)


class TimeSeriesDataPlot(DataLabelSvgMixin, DataPlot):
    """DataPlot of some `TimeSeries` data.
    """
    type = 'timeseries'
    data = 'timeseries'
    defaults = {'logy': False,
                'hline': list()}

    def __init__(self, *args, **kwargs):
        super(TimeSeriesDataPlot, self).__init__(*args, **kwargs)
        for c in self.channels:
            c._timeseries = True

    def add_state_segments(self, ax, **kwargs):
        """Add an `Axes` below the given ``ax`` displaying the `SummaryState`
        for this `TimeSeriesDataPlot`.

        Parameters
        ----------
        ax : `Axes`
            the set of `Axes` below which to display the state segments.
        **kwargs
            other keyword arguments will be passed to the
            :meth:`~gwpy.plotter.timeseries.TimeSeriesPlot.add_state_segments`
            method.
        """
        epoch = ax.get_epoch()
        kwargs.setdefault('edgecolor', 'darkgreen')
        kwargs.setdefault('facecolor', GREEN)
        kwargs.setdefault('known', {'facecolor': 'red', 'edgecolor': 'darkred'})
        sax = self.plot.add_state_segments(self.state, ax, plotargs=kwargs)
        ax.set_epoch(epoch)
        sax.set_epoch(epoch)
        sax.tick_params(axis='y', which='major', labelsize=12)
        sax.yaxis.set_ticks_position('none')
        sax.set_epoch(epoch)
        ax.set_epoch(epoch)
        return sax

    def init_plot(self, plot=TimeSeriesPlot, geometry=(1,1)):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesDataPlot`.
        """
        figsize = self.pargs.pop('figsize', [12, 6])
        self.plot, axes = subplots(
            nrows=geometry[0], ncols=geometry[1], sharex=True,
            subplot_kw={'projection': plot._DefaultAxesClass.name},
            FigureClass=plot, figsize=figsize, squeeze=True)
        if geometry[0] * geometry[1] == 1:
            axes = [axes]
        for ax in axes:
            if mode.MODE_NAME[mode.get_mode()] == 'MONTH':
                ax.set_xscale('days')
                ax.set_xlabel('_auto')
            ax.set_epoch(float(self.start))
            ax.grid(True, which='both')
        return self.plot, axes

    def finalize(self, outputfile=None, close=True, **savekwargs):
        plot = self.plot
        ax = plot.axes[0]
        if 'xlim' not in self.pargs:
            ax.set_xlim(float(self.start), float(self.end))
        return super(TimeSeriesDataPlot, self).finalize(
                   outputfile=outputfile, close=close, **savekwargs)

    def process(self, outputfile=None):
        """Read in all necessary data, and generate the figure.
        """
        (plot, axes) = self.init_plot()
        ax = axes[0]

        plotargs = self.parse_plot_kwargs()
        legendargs = self.parse_legend_kwargs()

        # add data
        channels, groups = zip(*self.get_channel_groups())
        for clist, pargs in zip(groups, plotargs):
            # pad data request to over-fill plots (no gaps at the end)
            if self.state and not self.all_data:
                valid = self.state.active
            elif clist[0].sample_rate:
                valid = SegmentList([self.span.protract(
                    1/clist[0].sample_rate.value)])
            else:
                valid = SegmentList([self.span])
            # get data
            data = [get_timeseries(c, valid, query=False)
                    for c in clist]
            if len(clist) > 1:
                data = [tsl.join(gap='pad', pad=numpy.nan) for tsl in data]
            flatdata = [ts for tsl in data for ts in tsl]
            # validate parameters
            for ts in flatdata:
                # double-check empty
                if (hasattr(ts, 'metadata') and
                        not 'x0' in ts.metadata) or not ts.x0:
                    ts.epoch = self.start
                # double-check log scales
                if self.pargs.get('logy', False):
                    ts.value[ts.value == 0] = 1e-100
            # set label
            try:
                label = pargs.pop('label')
            except KeyError:
                try:
                    label = label_to_latex(flatdata[0].name)
                except IndexError:
                    label = clist[0]
                else:
                    if self.fileformat == 'svg' and not label.startswith(
                            label_to_latex(
                            str(flatdata[0].channel)).split('.')[0]):
                        label += ' [%s]' % (
                            label_to_latex(str(flatdata[0].channel)))
            # plot groups or single TimeSeries
            if len(clist) > 1:
                ax.plot_timeseries_mmm(*data, label=label, **pargs)
            elif len(flatdata) == 0:
                ax.plot_timeseries(
                    data[0].EntryClass([], epoch=self.start, unit='s',
                                       name=label), label=label, **pargs)
            else:
                for ts in data[0]:
                    line = ax.plot_timeseries(ts, label=label, **pargs)[0]
                    label = None
                    pargs['color'] = line.get_color()

            # allow channel data to set parameters
            if len(flatdata):
                chan = get_channel(str(flatdata[0].channel))
            else:
                chan = get_channel(clist[0])
            if getattr(chan, 'amplitude_range', None) is not None:
                self.pargs.setdefault('ylim', chan.amplitude_range)

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
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if (len(channels) > 1 or plotargs[0].get('label', None) in
                [re.sub(r'(_|\\_)', r'\_', channels[0]), None]):
            plot.add_legend(ax=ax, **legendargs)

        # add extra axes and finalise
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        if self.state:
            self.add_state_segments(ax)
        return self.finalize(outputfile=outputfile)

register_plot(TimeSeriesDataPlot)


class SpectrogramDataPlot(TimeSeriesDataPlot):
    """DataPlot a Spectrogram
    """
    type = 'spectrogram'
    data = 'spectrogram'
    defaults = {'ratio': None,
                'format': None,
                'clim': None,
                'logcolor': False,
                'colorlabel': None}

    def __init__(self, *args, **kwargs):
        super(SpectrogramDataPlot, self).__init__(*args, **kwargs)
        self.ratio = self.pargs.pop('ratio')

    @property
    def pid(self):
        try:
            return self._pid
        except AttributeError:
            pid = super(SpectrogramDataPlot, self).pid
            if self.ratio:
                self._pid += '_%s_RATIO' % re_cchar.sub(
                    '_', str(self.ratio).upper())
            return self.pid

    @pid.setter
    def pid(self, id_):
        self._pid = str(id_)

    @pid.deleter
    def pid(self):
        del self._pid

    def process(self):
        # initialise
        (plot, axes) = self.init_plot()
        ax = axes[0]
        ax.grid(b=True, axis='y', which='major')
        channel = self.channels[0]

        # parse data arguments
        sdform = self.pargs.pop('format')
        clim = self.pargs.pop('clim')
        clog = self.pargs.pop('logcolor')
        clabel = self.pargs.pop('colorlabel')
        ratio = self.ratio

        # get cmap
        if ratio in ['median', 'mean']:
            self.pargs.setdefault('cmap', 'Spectral_r')
        cmap = self.pargs.pop('cmap', None)

        # get data
        if self.state and not self.all_data:
            valid = self.state.active
        else:
            valid = SegmentList([self.span])
        specgrams = get_spectrogram(channel, valid, query=False,
                                    format=sdform)
        # calculate ratio spectrum
        if len(specgrams) and (ratio in ['median', 'mean'] or
                               isinstance(ratio, int)):
            try:
                allspec = specgrams.join(gap='ignore')
            except ValueError as e:
                if 'units do not match' in str(e):
                    warnings.warn(str(e))
                    for spec in specgrams[1:]:
                        spec.unit = specgrams[0].unit
                    allspec = specgrams.join(gap='ignore')
                else:
                    raise
            if isinstance(ratio, int):
                ratio = allspec.percentile(ratio)
            else:
                ratio = getattr(allspec, ratio)(axis=0)

        # allow channel data to set parameters
        if getattr(channel, 'frequency_range', None) is not None:
            self.pargs.setdefault('ylim', channel.frequency_range)
            if isinstance(self.pargs['ylim'], Quantity):
                self.pargs['ylim'] = self.pargs['ylim'].value
        if (ratio is None and sdform in ['amplitude', 'asd'] and
                hasattr(channel, 'asd_range') and clim is None):
            clim = channel.asd_range
        elif (ratio is None and hasattr(channel, 'psd_range') and
              clim is None):
            clim = channel.psd_range

        # plot data
        for specgram in specgrams:
            if ratio is not None:
                specgram = specgram.ratio(ratio)
            ax.plot_spectrogram(specgram, cmap=cmap)

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
        if self.state:
            self.add_state_segments(ax)
        return self.finalize()

register_plot(SpectrogramDataPlot)


class SpectrumDataPlot(DataPlot):
    """Spectrum plot for a `SummaryTab`
    """
    type = 'spectrum'
    data = 'spectrum'
    defaults = {'logx': True,
                'logy': True,
                'format': None,
                'alpha': 0.1,
                'zorder': 1,
                'no_percentiles': False,
                'reference_linestyle': '--'}

    def process(self):
        pargs = self.pargs.copy()
        try:
            self._process()
        except OverflowError:
            self.pargs = pargs
            self.pargs['alpha'] = 0.0
            self._process()

    def _process(self):
        """Load all data, and generate this `SpectrumDataPlot`
        """
        plot = self.plot = SpectrumPlot(
            figsize=self.pargs.pop('figsize', [12, 6]))
        ax = plot.gca()
        ax.grid(b=True, axis='both', which='both')

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        label_to_latex(str(self.state))))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.993, va='top')

        # get spectrum format: 'amplitude' or 'power'
        sdform = self.pargs.pop('format')
        use_percentiles = str(
            self.pargs.pop('no_percentiles')).lower() == 'false'

        # parse plotting arguments
        plotargs = self.parse_plot_kwargs()
        legendargs = self.parse_legend_kwargs()

        # get reference arguments
        refs = []
        refkey = 'None'
        for key in sorted(self.pargs.keys()):
            if key == 'reference' or re.match('reference\d+\Z', key):
                refs.append(dict())
                refs[-1]['source'] = self.pargs.pop(key)
                refkey = key
            if re.match('%s[-_]' % refkey, key):
                refs[-1][key[len(refkey)+1:]] = self.pargs.pop(key)

        # add data
        for channel, pargs in zip(self.channels, plotargs):
            if self.state and not self.all_data:
                valid = self.state
            else:
                valid = SegmentList([self.span])
            data = get_spectrum(str(channel), valid, query=False,
                                format=sdform)

            # anticipate log problems
            if self.pargs['logx']:
                data = [s[1:] for s in data]
            if self.pargs['logy']:
                for sp in data:
                    sp.value[sp.value == 0] = 1e-100

            if use_percentiles:
                ax.plot_spectrum_mmm(*data, **pargs)
            else:
                pargs.pop('alpha', None)
                ax.plot_spectrum(data[0], **pargs)

            # allow channel data to set parameters
            if getattr(channel, 'frequency_range', None) is not None:
                self.pargs.setdefault('xlim', channel.frequency_range)
                if isinstance(self.pargs['xlim'], Quantity):
                    self.pargs['xlim'] = self.pargs['xlim'].value
            if (sdform in ['amplitude', 'asd'] and
                    hasattr(channel, 'asd_range')):
                self.pargs.setdefault('ylim', channel.asd_range)
            elif hasattr(channel, 'psd_range'):
                self.pargs.setdefault('ylim', channel.psd_range)

        # display references
        for i, ref in enumerate(refs):
            if 'source' in ref:
                source = ref.pop('source')
                try:
                    refspec = Spectrum.read(source)
                except IOError as e:
                    warnings.warn('IOError: %s' % str(e))
                except Exception as e:
                    # hack for old versions of GWpy
                    # TODO: remove me when GWSumm requires GWpy > 0.1
                    if 'Format could not be identified' in str(e):
                        refspec = Spectrum.read(source, format='dat')
                    else:
                        raise
                else:
                    ref.setdefault('zorder', -len(refs) + 1)
                    if 'filter' in ref:
                        refspec = refspec.filter(*ref.pop('filter'))
                    if 'scale' in ref:
                        refspec *= ref.pop('scale', 1)
                    ax.plot(refspec, **ref)

        # customise
        hlines = list(self.pargs.pop('hline', []))
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)

        # add horizontal lines to add
        if hlines:
            if not isinstance(hlines[-1], float):
                lineparams = hlines.pop(-1)
            else:
                lineparams = {'color':'r', 'linestyle': '--'}
        for yval in hlines:
            try:
                yval = float(yval)
            except ValueError:
                continue
            else:
                ax.plot(ax.get_xlim(), [yval, yval], **lineparams)

        if len(self.channels) > 1 or ax.legend_ is not None:
            plot.add_legend(ax=ax, **legendargs)
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)

        return self.finalize()

register_plot(SpectrumDataPlot)


class TimeSeriesHistogramPlot(DataPlot):
    """HistogramPlot from a Series
    """
    type = 'histogram'
    data = 'timeseries'
    defaults = {'ylabel': 'Rate [Hz]',
                'log': True,
                'histtype': 'stepfilled',
                'rwidth': 1}

    def init_plot(self, plot=HistogramPlot, geometry=None, figsize=[12, 6]):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesDataPlot`.
        """
        if geometry is None and self.pargs.pop('sep', False):
            if self.data == 'segments':
                geometry = (len(self.flags), 1)
            else:
                geometry = (len(self.channels), 1)
        elif geometry is None:
            geometry = (1, 1)
        self.plot, axes = subplots(
            nrows=geometry[0], ncols=geometry[1], sharex=True,
            subplot_kw={'projection': plot._DefaultAxesClass.name},
            FigureClass=plot, figsize=figsize, squeeze=True)
        if isinstance(axes, plot._DefaultAxesClass):
            axes = [axes]
        for ax in axes:
            ax.grid(True, which='both')
        return self.plot, axes

    def parse_plot_kwargs(self, **defaults):
        kwargs = super(TimeSeriesHistogramPlot, self).parse_plot_kwargs(
            **defaults)
        for histargs in kwargs:
            histargs.setdefault('logbins', self.pargs.get('logx', False))
            self.pargs.setdefault('logy', histargs.get('log', False))
            # set range as xlim
            if 'range' not in histargs and 'xlim' in self.pargs:
                histargs['range'] = self.pargs.get('xlim')
            # set alpha
            if len(self.channels) > 1:
                 histargs.setdefault('alpha', 0.7)
        return kwargs

    def process(self, outputfile=None):
        """Get data and generate the figure.
        """
        # get plot and axes
        (plot, axes) = self.init_plot()

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        label_to_latex(str(self.state))))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.993, va='top')

        # extract histogram arguments
        histargs = self.parse_plot_kwargs()

        # get data
        data = []
        for channel in self.channels:
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            data.append(get_timeseries(channel, valid, query=False).join(
                gap='ignore', pad=numpy.nan))
            # allow channel data to set parameters
            if hasattr(data[-1].channel, 'amplitude_range'):
                self.pargs.setdefault('xlim', data[-1].channel.amplitude_range)

        # get range
        if not 'range' in histargs[0]:
            l = axes[0].common_limits(data)
            for d in histargs:
                d['range'] = l

        # plot
        for ax, arr, pargs in zip(cycle(axes), data, histargs):
            if arr.size == 0:
                kwargs = dict(
                    (k, pargs[k]) for k in ['label', 'color'] if pargs.get(k))
                ax.plot([], **kwargs)
            else:
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
            if len(self.channels) > 1:
                    plot.add_legend(ax=ax, **legendargs)
        if len(axes) % 2 == 0 and axes[0].get_ylabel():
            label = axes[0].yaxis.label
            ax = axes[int(len(axes) // 2)-1]
            ax.set_ylabel(label.get_text())
            ax.yaxis.label.set_position((0, -.2 / len(axes)))
            if len(axes) != 2:
                label.set_text('')

        # add extra axes and finalise
        if not plot.colorbars:
            for ax in axes:
                plot.add_colorbar(ax=ax, visible=False)
        return self.finalize(outputfile=outputfile)

register_plot(TimeSeriesHistogramPlot)


class TimeSeriesHistogram2dDataPlot(TimeSeriesHistogramPlot):
    """DataPlot of the 2D histogram of two `TimeSeries`.
    """
    type = 'histogram2d'
    data = 'timeseries'
    defaults = {
        'logy': False,
        'hline': list(),
        'grid': 'both',
        'shading': 'flat',
        'cmap': 'inferno_r',
        'alpha': None,
        'edgecolors': 'None',
        'bins': 100,
        'normed': True
    }

    def __init__(self, *args, **kwargs):
        super(TimeSeriesHistogram2dDataPlot, self).__init__(*args, **kwargs)
        channels = self.channels
        if isinstance(channels, (list, tuple)) and len(channels) > 2:
            raise ValueError("Cannot generate TimeSeriesHistogram2dDataPlot "
                             " plot with more than 2 channels")

    def parse_hist_kwargs(self, **defaults):
        kwargs = {'bins': self.pargs.pop('bins'),
                  'normed': self.pargs.pop('normed')}
        if 'range' in self.pargs:
            ranges = [float(r) for r in self.pargs['range'].split(',')]
            kwargs['range'] = [[ranges[0], ranges[1]],
                               [ranges[2], ranges[3]]]
        elif 'xlim' in self.pargs and 'ylim' in self.pargs:
            xlim = self.pargs['xlim']
            ylim = self.pargs['ylim']
            kwargs['range'] = [[xlim[0], xlim[1]], [ylim[0], ylim[1]]]
        else:
            kwargs['range'] = None
        return kwargs

    def parse_pcmesh_kwargs(self, **defaults):
        kwargs = {
                  'cmap': self.pargs.pop('cmap'),
                  'edgecolors': self.pargs.pop('edgecolors'),
                  'shading': self.pargs.pop('shading'),
                  'alpha': self.pargs.pop('alpha')
                 }
        return kwargs

    def process(self, outputfile=None):
        """Get data and generate the figure.
        """
        # get histogram parameters
        plot, axes = self.init_plot()
        ax = axes[0]

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        label_to_latex(str(self.state))))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.993, va='top')
        # get data
        data = []
        for channel in self.channels:
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            data.append(get_timeseries(channel, valid, query=False).join(
                gap='ignore', pad=numpy.nan))
        if len(data) == 1:
            data.append(data[0])
        # allow channel data to set parameters
        self.pargs.setdefault('xlabel', label_to_latex(data[0].name))
        self.pargs.setdefault('ylabel', label_to_latex(data[1].name))
        if hasattr(data[0].channel, 'amplitude_range'):
            self.pargs.setdefault('xlim', data[0].channel.amplitude_range)
        if hasattr(data[1].channel, 'amplitude_range'):
            self.pargs.setdefault('ylim', data[1].channel.amplitude_range)
        # histogram
        hist_kwargs = self.parse_hist_kwargs()
        h, xedges, yedges = numpy.histogram2d(data[0], data[1],
                                              **hist_kwargs)
        h = numpy.ma.masked_where(h==0, h)
        x, y = numpy.meshgrid(xedges, yedges, copy=False, sparse=True)
        # plot
        pcmesh_kwargs = self.parse_pcmesh_kwargs()
        ax.pcolormesh(x, y, h.T, **pcmesh_kwargs)
        # customise plot
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                if key == 'grid':
                    if val == 'off':
                        ax.grid('off')
                    elif val in ['both', 'major', 'minor']:
                        ax.grid('on', which=val)
                    else:
                        raise ValueError("Invalid ax.grid argument; "
                                         "valid options are: 'off', "
                                         "'both', 'major' or 'minor'")
                else:
                    setattr(ax, key, val)
        # add extra axes and finalise
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        return self.finalize(outputfile=outputfile)

register_plot(TimeSeriesHistogram2dDataPlot)


class TimeSeriesScatterDataPlot(TimeSeriesHistogramPlot):
    """DataPlot of two `TimeSeries` with x & y corresponding to different
     channels and color to time. If requested, joins points with line.
    """
    type = 'scatter'
    data = 'timeseries'
    defaults = {
        'logy': False,
        'hline': list(),
        'grid': 'both',
        'marker': 'o',
        'cmap': 'inferno_r',
        'edgecolors': 'none',
        'alpha': None,
        'linewidths': None,
        'line': False,
        'line_color': 'black',
        'line_alpha': '0.2'
    }

    def __init__(self, *args, **kwargs):
        super(TimeSeriesScatterDataPlot, self).__init__(*args, **kwargs)
        channels = self.channels
        if isinstance(channels, (list, tuple)) and len(channels) > 2:
            raise ValueError("Cannot generate TimeSeriesScatterDataPlot "
                             " plot with more than 2 channels")

    def parse_scatter_kwargs(self, **defaults):
        kwargs = {
                  'marker': self.pargs.pop('marker'),
                  'cmap': self.pargs.pop('cmap'),
                  'edgecolors': self.pargs.pop('edgecolors'),
                  'alpha': self.pargs.pop('alpha'),
                  'linewidths': self.pargs.pop('linewidths')
                 }
        if self.pargs['line'] == 'False':
            self.pargs['line'] = False
        return kwargs

    def parse_line_kwargs(self, **defaults):
        """Pop keyword arguments for `Axes.plot` from the `pargs` for this Plot
        """
        plotargs = defaults.copy()
        plotargs.setdefault('label', self._parse_labels())
        kwargs = ['alpha', 'color', 'drawstyle', 'fillstyle', 'linestyle',
                  'linewidth', 'marker', 'markeredgecolor',
                  'markeredgewidth', 'markerfacecolor',
                  'markerfacecoloralt', 'markersize']
        for kwarg in kwargs:
            try:
                val = self.pargs.pop('line_%s' % kwarg)
            except KeyError:
                try:
                    val = self.pargs.pop('line_%ss' % kwarg)
                except KeyError:
                    val = None
            if val is not None:
                try:
                    val = eval(val)
                except Exception:
                    pass
                plotargs[kwarg] = val
        return plotargs

    def process(self, outputfile=None):
        """Get data and generate the figure.
        """
        plot, axes = self.init_plot()
        ax = axes[0]

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        label_to_latex(str(self.state))))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.993, va='top')
        # get data
        data = []
        for channel in self.channels:
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            data.append(get_timeseries(channel, valid, query=False).join(
                gap='ignore', pad=numpy.nan))
        if len(data) == 1:
            data.append(data[0])
        # allow channel data to set parameters
        self.pargs.setdefault('xlabel', label_to_latex(data[0].name))
        self.pargs.setdefault('ylabel', label_to_latex(data[1].name))
        if hasattr(data[0].channel, 'amplitude_range'):
            self.pargs.setdefault('xlim', data[0].channel.amplitude_range)
        else:
            self.pargs.setdefault('xlim', [min(data[0].value),
                                           max(data[0].value)])
        if hasattr(data[1].channel, 'amplitude_range'):
            self.pargs.setdefault('ylim', data[1].channel.amplitude_range)
        else:
            self.pargs.setdefault('ylim', [min(data[1].value),
                                           max(data[1].value)])
        # plot
        if self.pargs['line']:
            line_kwargs = self.parse_line_kwargs()
            ax.plot(data[0], data[1], zorder=1, **line_kwargs)
        scatter_kwargs = self.parse_scatter_kwargs()
        ax.scatter(data[0], data[1], c=data[0].times, zorder=2,
                   **scatter_kwargs)
        # customise plot
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                if key == 'grid':
                    if val == 'off':
                        ax.grid('off')
                    elif val in ['both', 'major', 'minor']:
                        ax.grid('on', which=val)
                    else:
                        raise ValueError("Invalid ax.grid argument; "
                                         "valid options are: 'off', "
                                         "'both', 'major' or 'minor'")
                else:
                    setattr(ax, key, val)
        # add extra axes and finalise
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        return self.finalize(outputfile=outputfile)

register_plot(TimeSeriesScatterDataPlot)


class SpectralVarianceDataPlot(SpectrumDataPlot):
    """SpectralVariance histogram plot for a `DataTab`
    """
    type = 'variance'
    data = 'spectrogram'
    defaults = {
        'logx': True,
        'logy': True,
        'reference_linestyle': '--',
        'log': True,
        'nbins': 100,
    }

    def __init__(self, channels, *args, **kwargs):
        if isinstance(channels, (list, tuple)) and len(channels) > 1:
            raise ValueError("Cannot generate SpectralVariance plot with "
                             "more than 1 channel")
        super(SpectralVarianceDataPlot, self).__init__(
            channels, *args, **kwargs)

    def parse_variance_kwargs(self):
        varargs = dict()
        for key in ['low', 'high', 'log', 'nbins', 'bins', 'density', 'norm']:
            if key in self.pargs:
                varargs[key] = self.pargs.pop(key)
        return varargs

    def _process(self):
        """Load all data, and generate this `SpectrumDataPlot`
        """
        plot = self.plot = SpectrumPlot(
            figsize=self.pargs.pop('figsize', [12, 6]))
        ax = plot.gca()

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        label_to_latex(str(self.state))))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.993, va='top')

        # parse plotting arguments
        cmap = self.pargs.pop('cmap', None)
        varargs = self.parse_variance_kwargs()
        plotargs = self.parse_plot_kwargs()[0]
        legendargs = self.parse_legend_kwargs()

        # get reference arguments
        refs = []
        refkey = 'None'
        for key in sorted(self.pargs.keys()):
            if key == 'reference' or re.match('reference\d+\Z', key):
                refs.append(dict())
                refs[-1]['source'] = self.pargs.pop(key)
                refkey = key
            if re.match('%s[-_]' % refkey, key):
                refs[-1][key[len(refkey)+1:]] = self.pargs.pop(key)

        # get channel arguments
        if hasattr(self.channels[0], 'asd_range'):
            low, high = self.channels[0].asd_range
            varargs.setdefault('low', low)
            varargs.setdefault('high', high)

        # calculate spectral variance and plot
        # pad data request to over-fill plots (no gaps at the end)
        if self.state and not self.all_data:
            valid = self.state.active
        else:
            valid = SegmentList([self.span])
        livetime = float(abs(valid))

        if livetime:
            plotargs.setdefault('vmin', 1/livetime)
        plotargs.setdefault('vmax', 1.)
        plotargs.pop('label')

        specgram = get_spectrogram(self.channels[0], valid, query=False,
                                    format='asd').join(gap='ignore')

        if specgram.size:
            asd = specgram.median(axis=0)
            asd.name = None
            variance = specgram.variance(**varargs)
            # normalize the variance
            variance /= livetime / specgram.dt.value
            # plot
            ax.plot(asd, color='grey', linewidth=0.3)
            m = ax.plot_variance(variance, cmap=cmap, **plotargs)
        #else:
        #    ax.scatter([1], [1], c=[1], visible=False, vmin=plotargs['vmin'],
        #               vmax=plotargs['vmax'], cmap=plotargs['cmap'])
        #plot.add_colorbar(ax=ax, log=True, label='Fractional time at amplitude')

        # allow channel data to set parameters
        if getattr(self.channels[0], 'frequency_range', None) is not None:
            self.pargs.setdefault('xlim', self.channels[0].frequency_range)
            if isinstance(self.pargs['xlim'], Quantity):
                self.pargs['xlim'] = self.pargs['xlim'].value
        if hasattr(self.channels[0], 'asd_range'):
            self.pargs.setdefault('ylim', self.channels[0].asd_range)

        # display references
        for i, ref in enumerate(refs):
            if 'source' in ref:
                source = ref.pop('source')
                try:
                    refspec = Spectrum.read(source)
                except IOError as e:
                    warnings.warn('IOError: %s' % str(e))
                except Exception as e:
                    # hack for old versions of GWpy
                    # TODO: remove me when GWSumm requires GWpy > 0.1
                    if 'Format could not be identified' in str(e):
                        refspec = Spectrum.read(source, format='dat')
                    else:
                        raise
                else:
                    if 'filter' in ref:
                        refspec = refspec.filter(*ref.pop('filter'))
                    if 'scale' in ref:
                        refspec *= ref.pop('scale', 1)
                    ax.plot(refspec, **ref)

        # customise
        hlines = list(self.pargs.pop('hline', []))
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)

        # add horizontal lines to add
        if hlines:
            if not isinstance(hlines[-1], float):
                lineparams = hlines.pop(-1)
            else:
                lineparams = {'color':'r', 'linestyle': '--'}
        for yval in hlines:
            try:
                yval = float(yval)
            except ValueError:
                continue
            else:
                ax.plot(ax.get_xlim(), [yval, yval], **lineparams)

        # set grid
        ax.grid(b=True, axis='both', which='both')

        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)

        return self.finalize()

register_plot(SpectralVarianceDataPlot)


class RayleighSpectrogramDataPlot(SpectrogramDataPlot):
    """Rayleigh statistic versino of `SpectrogramDataPlot`
    """
    type = 'rayleigh-spectrogram'
    data = 'rayleigh-spectrogram'
    defaults = {'ratio': None,
                'format': 'rayleigh',
                'clim': [0.25, 4],
                'cmap': 'BrBG_r',
                'logcolor': True,
                'colorlabel': 'Rayleigh statistic'}

register_plot(RayleighSpectrogramDataPlot)


class RayleighSpectrumDataPlot(SpectrumDataPlot):
    """Rayleigh statistic versino of `SpectrumDataPlot`
    """
    type = 'rayleigh-spectrum'
    data = 'rayleigh-spectrum'
    defaults = {'format': 'rayleigh',
                'logx': True,
                'logy': True,
                'alpha': 0.1,
                'zorder': 1,
                'no_percentiles': True,
                'reference_linestyle': '--'}

register_plot(RayleighSpectrumDataPlot)
