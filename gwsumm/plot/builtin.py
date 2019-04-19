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

import os.path
import warnings
from itertools import cycle

from six import string_types

import numpy

from matplotlib.colors import LogNorm

from astropy.units import Quantity

from gwpy.plot.colors import tint
from gwpy.plot.gps import GPSTransform
from gwpy.segments import SegmentList

from .. import (globalv, io)
from ..mode import (Mode, get_mode)
from ..utils import re_cchar
from ..data import (get_timeseries, get_spectrogram,
                    get_coherence_spectrogram, get_spectrum,
                    get_coherence_spectrum)
from ..state import ALLSTATE
from .registry import (get_plot, register_plot)
from .mixins import DataLabelSvgMixin
from .utils import usetex_tex

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

DataPlot = get_plot('data')
GREEN = (0.2, 0.8, 0.2)


class TimeSeriesDataPlot(DataLabelSvgMixin, DataPlot):
    """DataPlot of some `TimeSeries` data.
    """
    type = 'timeseries'
    data = 'timeseries'
    defaults = DataPlot.defaults.copy()
    defaults.update({
        'xscale': 'auto-gps',
        'yscale': 'linear',
    })

    def __init__(self, *args, **kwargs):
        super(TimeSeriesDataPlot, self).__init__(*args, **kwargs)
        if self.data == 'timeseries':
            for c in self.channels:
                c._timeseries = True

    def _update_defaults_from_channels(self):
        for chan in self.channels:
            if getattr(chan, 'amplitude_range', None) is not None:
                self.pargs.setdefault('ylim', chan.amplitude_range)
                break

    # -- utilities ------------------------------

    def add_state_segments(self, ax, visible=None, **kwargs):
        """Add an `Axes` below the given ``ax`` displaying the `SummaryState`
        for this `TimeSeriesDataPlot`.

        Parameters
        ----------
        ax : `Axes`
            the set of `Axes` below which to display the state segments.

        visible : `bool`, optional
            whether or not to display the axes, or just make space for them,
            default is `None`, meaning a dynamic choice based on the state

        **kwargs
            other keyword arguments will be passed to the
            :meth:`~gwpy.plot.Plot.add_segments_bar`
            method.
        """
        # allow user to disable the state segments axes
        if self.pargs.pop('no-state-segments', False):
            visible = False
        if visible is None and self.state is not None and (
                self.state.name.lower() != ALLSTATE):
            visible = True
        if visible:
            sax = self.plot.add_segments_bar(self.state, ax, height=.14,
                                             pad=.1,  **kwargs)
            sax.tick_params(axis='y', which='major', labelsize=12)
            sax.yaxis.set_ticks_position('none')
            sax.set_ylim(-.4, .4)
            return sax
        else:
            self.plot.subplots_adjust(bottom=0.17)
            return None

    def add_future_shade(self, gps=None, facecolor='gray', alpha=.1,
                         **kwargs):
        """Shade those parts of the figure that display times in the future
        """
        end = float(self.end)
        # get time 'now'
        if gps is None:
            gps = globalv.NOW
        # allow user to override
        if self.pargs.pop('no-future-shade', False) or end <= gps:
            return
        # shade time axes
        for ax in filter(
                lambda ax: isinstance(ax.xaxis.get_transform(), GPSTransform),
                self.plot.axes):
            ax.axvspan(gps, end, facecolor=facecolor, alpha=alpha, **kwargs)

    # -- init/finalize --------------------------

    def init_plot(self, *args, **kwargs):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesDataPlot`.
        """
        epoch = kwargs.pop('epoch', self.pargs.pop('epoch', None))
        plot = super(TimeSeriesDataPlot, self).init_plot(*args, **kwargs)
        for ax in plot.axes:
            if get_mode() == Mode.month:
                ax.set_xscale('days')
            if isinstance(ax.xaxis.get_transform(), GPSTransform):
                ax.set_epoch(float(epoch if epoch is not None else self.start))
                if ax.get_autoscalex_on():
                    ax.set_xlim(float(self.start), float(self.end))
            ax.grid(True, which='both')
        return plot

    # -- main draw method -----------------------

    def draw(self, outputfile=None):
        """Read in all necessary data, and generate the figure.
        """
        plot = self.init_plot()
        ax = plot.gca()

        plotargs = self.parse_plot_kwargs()
        legendargs = self.parse_legend_kwargs()

        # add data
        channels, groups = list(zip(*self.get_channel_groups()))
        for clist, pargs in list(zip(groups, plotargs)):
            # get data
            valid = self._get_data_segments(clist[0])
            data = [get_timeseries(c, valid, query=False)
                    for c in clist]

            if len(clist) > 1:
                data = [tsl.join(gap='pad', pad=numpy.nan) for tsl in data]
            flatdata = [ts for tsl in data for ts in tsl]
            # validate parameters
            for ts in flatdata:
                # double-check empty
                if ts.x0 is None:
                    ts.epoch = self.start
                # double-check log scales
                if self.logy:
                    ts.value[ts.value == 0] = 1e-100
            # set label
            try:
                label = pargs.pop('label')
            except KeyError:
                try:
                    label = usetex_tex(flatdata[0].name)
                except IndexError:
                    label = clist[0]
                else:
                    if self.fileformat == 'svg' and not label.startswith(
                            usetex_tex(
                            str(flatdata[0].channel)).split('.')[0]):
                        label += ' [%s]' % (
                            usetex_tex(str(flatdata[0].channel)))

            # plot groups or single TimeSeries
            if len(clist) > 1:
                data[1].name = None  # force no labels for shades
                data[2].name = None
                ax.plot_mmm(*data, label=label, **pargs)
            elif len(flatdata) == 0:
                ax.plot(data[0].EntryClass([], epoch=self.start, unit='s',
                                           name=label),
                        label=label, **pargs)
            else:
                for ts in data[0]:
                    line, = ax.plot(ts, label=label, **pargs)
                    label = None
                    pargs['color'] = line.get_color()

        # customise plot
        self.add_hvlines()
        self.apply_parameters(ax, **self.pargs)

        # add legend
        if ax.get_legend_handles_labels()[0]:
            ax.legend(**legendargs)

        self.add_state_segments(ax)
        self.add_future_shade()

        return self.finalize(outputfile=outputfile)

    def _get_data_segments(self, channel):
        """Get data segments for this plot
        """
        if self.state and not self.all_data:
            return self.state.active
        if channel.sample_rate is not None:
            return SegmentList([self.span.protract(
                1/channel.sample_rate.value)])
        return SegmentList([self.span])


register_plot(TimeSeriesDataPlot)


class SpectrogramDataPlot(TimeSeriesDataPlot):
    """DataPlot a Spectrogram
    """
    type = 'spectrogram'
    data = 'spectrogram'
    defaults = TimeSeriesDataPlot.defaults.copy()
    defaults.update({
        'yscale': 'log',
        'ylabel': 'Frequency [Hz]',
        'ratio': None,
        'format': None,
        'rasterized': True,
    })

    def __init__(self, *args, **kwargs):
        super(SpectrogramDataPlot, self).__init__(*args, **kwargs)
        self.ratio = self.pargs.pop('ratio')
        # set default colour-map for median ratio
        if self.ratio in ('median', 'mean'):
            self.pargs.setdefault('cmap', 'Spectral_r')

    def _update_defaults_from_channels(self):
        for channel in self.channels:
            self.pargs.setdefault('ylim', channel.frequency_range)
            if isinstance(self.pargs['ylim'], Quantity):
                self.pargs['ylim'] = self.pargs['ylim'].value

            if self.ratio is None and self.pargs.get('clim') is None:
                if (self.pargs.get('format') in ('amplitude', 'asd') and
                        hasattr(channel, 'asd_range')):
                    self.pargs['clim'] = channel.asd_range
                elif hasattr(channel, 'psd_range'):
                    self.pargs['clim'] = channel.psd_range

    @property
    def pid(self):
        try:
            return self._pid
        except AttributeError:
            super(SpectrogramDataPlot, self).pid
            if (isinstance(self.ratio, string_types) and
                    os.path.isfile(self.ratio)):
                self._pid += '_REFERENCE_RATIO'
            elif self.ratio:
                self._pid += '_%s_RATIO' % re_cchar.sub(
                    '_', str(self.ratio).upper())
            return self.pid

    @pid.setter
    def pid(self, id_):
        self._pid = str(id_)

    @pid.deleter
    def pid(self):
        del self._pid

    def get_ratio(self, specgrams):
        ratio = self.ratio
        # calculate ratio spectrum
        if len(specgrams) and (
                ratio in ['median', 'mean'] or isinstance(ratio, int)):
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
                return allspec.percentile(ratio)
            else:
                return getattr(allspec, ratio)(axis=0)
        elif isinstance(ratio, string_types) and os.path.isfile(ratio):
            try:
                return io.read_frequencyseries(ratio)
            except IOError as e:  # skip if file can't be read
                warnings.warn('IOError: %s' % str(e))
        return ratio

    def draw(self):
        # initialise
        plot = self.init_plot()
        ax = plot.gca()
        ax.grid(b=True, axis='y', which='major')
        channel = self.channels[0]

        # parse data arguments
        sdform = self.pargs.pop('format')

        # parse colorbar arguments
        clabel = self.pargs.pop('colorlabel', '')
        clim = self.pargs.pop('clim', None)
        clog = self.pargs.pop('logcolor', False)

        # parse plotting arguments
        if clim:  # clim -> (vmin, vmax)
            vmin, vmax = clim
            self.pargs.setdefault('vmin', vmin)
            self.pargs.setdefault('vmax', vmax)
        if clog:  # logcolor -> norm
            self.pargs.setdefault('norm', 'log')
        plotargs = self.parse_plot_kwargs()[0]  # only one channel

        # rework norm='log' into a LogNorm object
        # (this is only 'required' in the case of no data, when ax.scatter
        #  gets called manually below)
        if plotargs.get('norm', None) == 'log':
            vmin = self.pargs.get('vmin')
            vmax = self.pargs.get('vmax')
            plotargs['norm'] = LogNorm(vmin=vmin, vmax=vmax)

        # get data
        if self.state and not self.all_data:
            valid = self.state.active
        else:
            valid = SegmentList([self.span])
        if self.type == 'coherence-spectrogram':
            specgrams = get_coherence_spectrogram(self.channels, valid,
                                                  query=False)
        else:
            try:
                specgrams = get_spectrogram(channel, valid, query=False,
                                            format=sdform)
            except ValueError as exc:
                if 'need more than 0 values' not in str(exc):
                    raise
                # attempted to do math but one input has zero size
                specgrams = []

        # get ratio as FrequencySeries
        ratio = self.get_ratio(specgrams)

        # plot data
        for i, specgram in enumerate(specgrams):

            # calculate ratio
            if ratio is not None:
                specgram = specgram.ratio(ratio)

            # undo demodulation and crop frequencies
            ylim = self.pargs.get('ylim', None)
            specgram = undo_demodulation(specgram, channel, ylim)
            if ylim is not None:
                specgram = specgram.crop_frequencies(*ylim)

            # plot
            ax.imshow(specgram, **plotargs)

        # add colorbar
        if len(specgrams) == 0:
            ax.pcolormesh([1, 10], [1, 10], [[1, 10]], visible=False,
                          **plotargs)
        ax.colorbar(label=clabel)

        # customise and finalise
        self.apply_parameters(ax, **self.pargs)
        self.add_state_segments(ax)
        self.add_future_shade()

        return self.finalize()


register_plot(SpectrogramDataPlot)


class CoherenceSpectrogramDataPlot(SpectrogramDataPlot):
    """DataPlot a Spectrogram of the coherence between two channels
    """
    type = 'coherence-spectrogram'
    data = 'coherence-spectrogram'
    defaults = SpectrogramDataPlot.defaults.copy()
    defaults.update({
        'ratio': None,
        'format': None,
        'clim': None,
        'logcolor': False,
        'colorlabel': None,
    })


register_plot(CoherenceSpectrogramDataPlot)


class SpectrumDataPlot(DataPlot):
    """Spectrum plot for a `SummaryTab`
    """
    type = 'spectrum'
    data = 'spectrum'
    defaults = DataPlot.defaults.copy()
    defaults.update({
        'xscale': 'log',
        'yscale': 'log',
        'format': None,
        'zorder': 1,
        'no-percentiles': False,
        'reference-linestyle': '--',
    })

    def _update_defaults_from_channels(self):
        for channel in self.channels:
            if getattr(channel, 'frequency_range', None) is not None:
                self.pargs.setdefault('xlim', channel.frequency_range)
                if isinstance(self.pargs['xlim'], Quantity):
                    self.pargs['xlim'] = self.pargs['xlim'].value
            if (self.pargs.get('format') in ['amplitude', 'asd'] and
                    hasattr(channel, 'asd_range')):
                self.pargs.setdefault('ylim', channel.asd_range)
            elif hasattr(channel, 'psd_range'):
                self.pargs.setdefault('ylim', channel.psd_range)

    def draw(self):
        pargs = self.pargs.copy()
        try:
            self._draw()
        except OverflowError:
            self.pargs = pargs
            self.pargs['alpha'] = 0.0
            self._draw()

    def _draw(self):
        """Load all data, and generate this `SpectrumDataPlot`
        """
        plot = self.init_plot()
        ax = plot.gca()
        ax.grid(b=True, axis='both', which='both')

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        usetex_tex(str(self.state))))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.993, va='top')

        # get spectrum format: 'amplitude' or 'power'
        sdform = self.pargs.pop('format')
        if sdform == 'rayleigh':
            method = 'rayleigh'
        else:
            method = None
        use_percentiles = str(
            self.pargs.pop('no-percentiles')).lower() == 'false'

        # parse plotting arguments
        plotargs = self.parse_plot_kwargs()
        legendargs = self.parse_legend_kwargs()
        use_legend = False

        # get reference arguments
        refs = self.parse_references()

        # add data
        if self.type == 'coherence-spectrum':
            iterator = list(zip(self.channels[0::2], self.channels[1::2],
                                plotargs))
        else:
            iterator = list(zip(self.channels, plotargs))

        for chantuple in iterator:
            channel = chantuple[0]
            channel2 = chantuple[1]
            pargs = chantuple[-1]

            if self.state and not self.all_data:
                valid = self.state
            else:
                valid = SegmentList([self.span])

            if self.type == 'coherence-spectrum':
                data = get_coherence_spectrum([str(channel), str(channel2)],
                                              valid, query=False)
            else:
                try:
                    data = get_spectrum(str(channel), valid, query=False,
                                        format=sdform, method=method)
                except ValueError as exc:
                    # math op failed beacuse one of the datasets is empty
                    if (
                            'could not be broadcast' in str(exc) and
                            '(0,)' in str(exc)
                    ):
                        data = []
                    else:
                        raise

            # undo demodulation
            data = list(data)
            for i, spec in enumerate(data):
                data[i] = undo_demodulation(spec, channel,
                                            self.pargs.get('xlim', None))

            # anticipate log problems
            if self.logx:
                data = [s[1:] for s in data]
            if self.logy:
                for sp in data:
                    sp.value[sp.value == 0] = 1e-100

            if 'label' in pargs:
                use_legend = True

            if data and use_percentiles:
                _, minline, maxline, _ = ax.plot_mmm(*data, **pargs)
                # make min, max lines lighter:
                minline.set_alpha(pargs.get('alpha', .1) * 2)
                maxline.set_alpha(pargs.get('alpha', .1) * 2)
            elif data:
                ax.plot(data[0], **pargs)

        # display references
        for source, refparams in refs.items():
            refspec = io.read_frequencyseries(source)
            refparams.setdefault('zorder', -len(refs) + 1)
            if 'filter' in refparams:
                refspec = refspec.filter(*refparams.pop('filter'))
            if 'scale' in refparams:
                refspec *= refparams.pop('scale', 1)
            if 'label' in refparams:
                use_legend = True
            ax.plot(refspec, **refparams)

        # customise
        self.add_hvlines()
        self.apply_parameters(ax, **self.pargs)
        if use_legend or len(self.channels) > 1 or ax.legend_ is not None:
            ax.legend(**legendargs)

        return self.finalize()

    def parse_references(self, prefix=r'reference(\d+)?\Z'):
        """Parse parameters for displaying one or more reference traces
        """
        return self.parse_list('reference')


register_plot(SpectrumDataPlot)


class CoherenceSpectrumDataPlot(SpectrumDataPlot):
    """Coherence pectrum plot for a `SummaryTab`
    """
    type = 'coherence-spectrum'
    data = 'coherence-spectrogram'
    defaults = SpectrumDataPlot.defaults.copy()
    defaults.update({
        'yscale': 'linear',
        'alpha': 0.1,
    })

    # override this to allow us to set the legend manually
    def _parse_labels(self, defaults=None):
        return self.pargs.pop('labels', defaults)

    def get_channel_groups(self):
        """Hi-jacked method to return pairs of channels

        For the `CoherenceSpectrumDataPlot` this method is only used in
        determining how to separate lists of plotting argument given by
        the user.
        """
        all_ = self.allchannels
        return [(all_[i], all_[i:i+2]) for i in range(0, len(all_), 2)]


register_plot(CoherenceSpectrumDataPlot)


class TimeSeriesHistogramPlot(DataPlot):
    """HistogramPlot from a Series
    """
    type = 'histogram'
    data = 'timeseries'
    defaults = DataPlot.defaults.copy()
    defaults.update({
        'ylabel': 'Rate [Hz]',
        'log': True,
        'histtype': 'stepfilled',
        'rwidth': 1,
        'bottom': 1e-300,
    })

    def _update_defaults_from_channels(self):
        for channel in self.channels:
            if hasattr(channel, 'amplitude_range'):
                self.pargs.setdefault('xlim', channel.amplitude_range)
                break

    def init_plot(self, geometry=None, **kwargs):
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
        plot = super(TimeSeriesHistogramPlot, self).init_plot(
            geometry=geometry, **kwargs)
        for ax in plot.axes:
            ax.grid(True, which='both')
        return plot

    def parse_plot_kwargs(self, **defaults):
        kwargs = super(TimeSeriesHistogramPlot, self).parse_plot_kwargs(
            **defaults)
        for histargs in kwargs:
            histargs.setdefault('logbins', self.logx)
            logy = histargs.get('log', False)
            self.pargs.setdefault('yscale', 'log' if logy else 'linear')
            # set range as xlim
            if 'range' not in histargs and 'xlim' in self.pargs:
                histargs['range'] = self.pargs.get('xlim')
            # set alpha
            if len(self.channels) > 1:
                histargs.setdefault('alpha', 0.7)
        return kwargs

    def draw(self, outputfile=None):
        """Get data and generate the figure.
        """
        # get plot and axes
        plot = self.init_plot()
        axes = plot.axes

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        usetex_tex(str(self.state))))
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

        # plot
        for ax, arr, pargs in zip(cycle(axes), data, histargs):
            # set range if not given
            if pargs.get('range') is None:
                pargs['range'] = self._get_range(
                    data,
                    # use range from first dataset if already calculated
                    range=histargs[0].get('range'),
                    # use xlim if manually set (user or INI)
                    xlim=None if ax.get_autoscalex_on() else ax.get_xlim(),
                )

            # plot histogram
            _, _, patches = ax.hist(arr, **pargs)

            # update edge color of histogram to be tinted version of face
            if pargs.get('histtype', None) == 'stepfilled':
                for p in patches:
                    if not p.get_edgecolor()[3]:
                        p.set_edgecolor(tint(p.get_facecolor(), .7))

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
            if len(self.channels) > 1:
                ax.legend(**legendargs)
        if len(axes) % 2 == 0 and axes[0].get_ylabel():
            label = axes[0].yaxis.label
            ax = axes[int(len(axes) // 2)-1]
            ax.set_ylabel(label.get_text())
            ax.yaxis.label.set_position((0, -.2 / len(axes)))
            if len(axes) != 2:
                label.set_text('')

        # add extra axes and finalise
        return self.finalize(outputfile=outputfile)

    def _get_range(self, data, range=None, xlim=None):
        if range is not None or xlim is not None:
            return range or xlim
        try:
            return numpy.min(data), numpy.max(data)
        except ValueError as exc:
            if not str(exc).startswith('zero-size array'):
                raise
        return None


register_plot(TimeSeriesHistogramPlot)


class TimeSeriesHistogram2dDataPlot(TimeSeriesHistogramPlot):
    """DataPlot of the 2D histogram of two `TimeSeries`.
    """
    type = 'histogram2d'
    data = 'timeseries'
    defaults = TimeSeriesHistogramPlot.defaults.copy()
    defaults.update({
        'yscale': 'linear',
        'grid': 'both',
        'shading': 'flat',
        'cmap': 'inferno_r',
        'alpha': None,
        'edgecolors': 'None',
        'bins': 100,
        'normed': True
    })

    def __init__(self, *args, **kwargs):
        super(TimeSeriesHistogram2dDataPlot, self).__init__(*args, **kwargs)
        channels = self.channels
        if isinstance(channels, (list, tuple)) and len(channels) > 2:
            raise ValueError("Cannot generate TimeSeriesHistogram2dDataPlot "
                             " plot with more than 2 channels")

    def _update_defaults_from_channels(self):
        c1, c2 = self.channels
        self.pargs.setdefault('xlabel', usetex_tex(str(c1)))
        self.pargs.setdefault('ylabel', usetex_tex(str(c2)))
        if hasattr(c1, 'amplitude_range'):
            self.pargs.setdefault('xlim', c1.amplitude_range)
        if hasattr(c2, 'amplitude_range'):
            self.pargs.setdefault('ylim', c2.amplitude_range)

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

    def draw(self, outputfile=None):
        """Get data and generate the figure.
        """
        # get histogram parameters
        plot = self.init_plot()
        ax = plot.gca()

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        usetex_tex(str(self.state))))
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
        # histogram
        hist_kwargs = self.parse_hist_kwargs()
        h, xedges, yedges = numpy.histogram2d(data[0], data[1], **hist_kwargs)
        h = numpy.ma.masked_where(h == 0, h)
        x, y = numpy.meshgrid(xedges, yedges, copy=False, sparse=True)
        # plot
        pcmesh_kwargs = self.parse_pcmesh_kwargs()
        ax.pcolormesh(x, y, h.T, **pcmesh_kwargs)

        # customise plot
        self.apply_parameters(ax, **self.pargs)

        return self.finalize(outputfile=outputfile)


register_plot(TimeSeriesHistogram2dDataPlot)


class SpectralVarianceDataPlot(SpectrumDataPlot):
    """SpectralVariance histogram plot for a `DataTab`
    """
    type = 'variance'
    data = 'spectrogram'
    defaults = SpectrumDataPlot.defaults.copy()
    defaults.update({
        'xscale': 'log',
        'yscale': 'log',
        'reference-linestyle': '--',
        'log': True,
        'nbins': 100,
    })

    def __init__(self, channels, *args, **kwargs):
        if isinstance(channels, (list, tuple)) and len(channels) > 1:
            raise ValueError("Cannot generate SpectralVariance plot with "
                             "more than 1 channel")
        super(SpectralVarianceDataPlot, self).__init__(
            channels, *args, **kwargs)

    def _update_defaults_from_channels(self):
        chan = self.channels[0]

        if getattr(chan, 'frequency_range', None) is not None:
            self.pargs.setdefault('xlim', chan.frequency_range)
            if isinstance(self.pargs['xlim'], Quantity):
                self.pargs['xlim'] = self.pargs['xlim'].value

        if hasattr(chan, 'asd_range'):
            self.pargs.setdefault('ylim', chan.asd_range)

        if hasattr(chan, 'asd_range'):
            low, high = chan.asd_range
            self.pargs.setdefault('low', low)
            self.pargs.setdefault('high', high)

    def parse_variance_kwargs(self):
        varargs = dict()
        for key in ['low', 'high', 'log', 'nbins', 'bins', 'density', 'norm']:
            if key in self.pargs:
                varargs[key] = self.pargs.pop(key)
        return varargs

    def _draw(self):
        """Load all data, and generate this `SpectrumDataPlot`
        """
        plot = self.init_plot()
        ax = plot.gca()

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        usetex_tex(str(self.state))))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.993, va='top')

        # parse plotting arguments
        cmap = self.pargs.pop('cmap', None)
        varargs = self.parse_variance_kwargs()
        plotargs = self.parse_plot_kwargs()[0]

        # get reference arguments
        refs = self.parse_references()

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
            # undo demodulation
            variance = undo_demodulation(variance, self.channels[0],
                                         ax.get_xlim())

            # plot
            ax.plot(asd, color='grey', linewidth=0.3)
            ax.imshow(variance, cmap=cmap, **plotargs)

        # display references
        for source, refparams in refs.items():
            refspec = io.read_frequencyseries(source)
            if 'filter' in refparams:
                refspec = refspec.filter(*refparams.pop('filter'))
            if 'scale' in refparams:
                refspec *= refparams.pop('scale', 1)
            ax.plot(refspec, **refparams)

        # customise
        self.add_hvlines()
        self.apply_parameters(ax, **self.pargs)
        ax.grid(b=True, axis='both', which='both')

        return self.finalize()


register_plot(SpectralVarianceDataPlot)


class RayleighSpectrogramDataPlot(SpectrogramDataPlot):
    """Rayleigh statistic versino of `SpectrogramDataPlot`
    """
    type = 'rayleigh-spectrogram'
    data = 'rayleigh-spectrogram'
    defaults = SpectrogramDataPlot.defaults.copy()
    defaults.update({
        'ratio': None,
        'format': 'rayleigh',
        'clim': [0.25, 4],
        'cmap': 'BrBG_r',
        'colorlabel': 'Rayleigh statistic',
    })


register_plot(RayleighSpectrogramDataPlot)


class RayleighSpectrumDataPlot(SpectrumDataPlot):
    """Rayleigh statistic versino of `SpectrumDataPlot`
    """
    type = 'rayleigh-spectrum'
    data = 'rayleigh-spectrum'
    defaults = {'format': 'rayleigh',
                'xscale': 'log',
                'yscale': 'log',
                'alpha': 0.1,
                'zorder': 1,
                'no-percentiles': True,
                'reference-linestyle': '--'}


register_plot(RayleighSpectrumDataPlot)


def undo_demodulation(spec, channel, limits=None):
    if spec.size == 0:
        return spec
    # undo demodulation
    try:
        demod = channel.demodulation
    except AttributeError:
        return spec
    else:
        spec = spec[:]  # views data with copied metadata
        del spec.frequencies
        spec.f0 = demod
        # if physical frequency-range is below demod, get negative df
        try:
            low, high = channel.frequency_range
        except (AttributeError, TypeError):
            try:
                low, high = limits
            except TypeError:
                return spec
        high = Quantity(high, 'Hz')
        if high < spec.f0:
            if spec.ndim > 1:  # Spectrogram
                spec.value[:] = numpy.fliplr(spec.value)
            else:  # FrequencySeries
                spec.value[:] = spec.value[::-1]
            spec.df *= -1
            spec.frequencies = spec.frequencies[::-1]
        return spec
