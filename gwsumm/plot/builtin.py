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
from itertools import (izip, cycle)

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from gwpy.spectrum import Spectrum
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
from .registry import (get_plot, register_plot)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

DataPlot = get_plot('data')


class TimeSeriesDataPlot(DataPlot):
    """DataPlot of some `TimeSeries` data.
    """
    type = 'timeseries'
    defaults = {'logy': False,
                'hline': list()}

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
        kwargs.setdefault('edgecolor', 'black')
        kwargs.setdefault('facecolor', (0.2, 0.8, 0.2))
        kwargs.setdefault('valid', {'facecolor': 'red'})
        sax = self.plot.add_state_segments(self.state, ax, plotargs=kwargs)
        sax.tick_params(axis='y', which='major', labelsize=12)
        sax.set_epoch(float(self.start))
        return sax

    def init_plot(self, plot=TimeSeriesPlot):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesDataPlot`.
        """
        self.plot = plot()
        ax = self.plot.gca()
        ax.set_epoch(float(self.start))
        return self.plot, ax

    def finalize(self, outputfile=None):
        plot = self.plot
        ax = plot.axes[0]
        if 'xlim' not in self.pargs:
            ax.set_xlim(float(self.start), float(self.end))
        return super(TimeSeriesDataPlot, self).finalize(
                   outputfile=outputfile)

    def process(self, outputfile=None):
        """Read in all necessary data, and generate the figure.
        """
        (plot, ax) = self.init_plot()

        plotargs = self.parse_plot_kwargs()
        legendargs = self.parse_legend_kwargs()

        # add data
        mmmchans = self.get_channel_groups()
        for channels, pargs in zip(mmmchans.values(), plotargs):
            # pad data request to over-fill plots (no gaps at the end)
            if self.state and not self.all_data:
                valid = self.state.active
            elif channels[0].sample_rate:
                valid = SegmentList([self.span.protract(
                    1/channels[0].sample_rate.value)])
            else:
                valid = SegmentList([self.span])
            # get data
            data = [get_timeseries(c, valid, query=False).join(
                        gap='ignore', pad=numpy.nan)
                    for c in channels]
            # double-check empty
            for ts in data:
                if not 'x0' in ts.metadata:
                    ts.epoch = self.start
            # double-check log scales
            if self.pargs['logy']:
                for ts in data:
                    ts[ts.data == 0] = 1e-100
            # plot groups or single TimeSeries
            if len(channels) > 1:
                ax.plot_timeseries_mmm(*data, **pargs)
            else:
                ax.plot_timeseries(data[0], **pargs)

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
        if (len(mmmchans) > 1 or plotargs[0].get('label', None) in
                [re.sub(r'(_|\\_)', r'\_', mmmchans.keys()[0]), None]):
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
    defaults = {'ratio': None,
                'format': None,
                'clim': None,
                'cmap': 'jet',
                'logcolor': False,
                'colorlabel': None}

    def __init__(self, *args, **kwargs):
        super(SpectrogramDataPlot, self).__init__(*args, **kwargs)
        self.ratio = self.pargs.pop('ratio')

    @property
    def tag(self):
        try:
            return self._tag
        except AttributeError:
            tag = super(SpectrogramDataPlot, self).tag
            if self.ratio:
                tag += '_%s_RATIO' % re_cchar.sub('_', str(self.ratio).upper())
            return tag

    @tag.setter
    def tag(self, filetag):
        self._tag = filetag or self.type.upper()

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
        if self.state and not self.all_data:
            valid = self.state.active
        else:
            valid = SegmentList([self.span])
        specgrams = get_spectrogram(self.channels[0], valid, query=False,
                                    format=sdform)
        # calculate ratio spectrum
        if ratio in ['median', 'mean'] or isinstance(ratio, int):
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
        if self.state:
            self.add_state_segments(ax)
        return self.finalize()

register_plot(SpectrogramDataPlot)


class SegmentDataPlot(TimeSeriesDataPlot):
    """Segment plot of one or more `DataQualityFlags <DataQualityFlag>`.
    """
    type = 'segments'
    defaults = {'mask': None,
                'color': None,
                'on_is_bad': False,
                'insetlabels': 'inset',
                'edgecolor': 'black'}

    def __init__(self, flags, start, end, state=None, outdir='.', tag=None,
                 **kwargs):
        super(SegmentDataPlot, self).__init__([], start, end, state=state,
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
        """Interferometer set for this `SegmentDataPlot`
        """
        allflags = [f for flag in self.flags for f in flag.split(',')]
        return set([f[:2] for f in allflags])

    @property
    def tag(self):
        """File tag for this `DataPlot`.
        """
        try:
            return self._tag
        except AttributeError:
            state = re_cchar.sub('_',
                                 self.state is None and 'MULTI' or
                                 self.state.name)
            hash_ = hashlib.md5("".join(map(str, self.flags))).hexdigest()[:6]
            return '_'.join([state, hash_, self.type]).upper()

    @tag.setter
    def tag(self, filetag):
        self._tag = filetag or self.type.upper()

    @classmethod
    def from_ini(cls, config, section, start, end, flags=None, state=ALLSTATE,
                 **kwargs):
        new = super(SegmentDataPlot, cls).from_ini(config, section, start,
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
        ax.set_insetlabels(self.pargs.pop('insetlabels', True))
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        labels = map(lambda s: re_quote.sub('', str(s).strip('\n ')), labels)

        # extract plotting arguments
        mask = self.pargs.pop('mask')
        activecolor, validcolor = self.get_segment_color()
        edgecolor = self.pargs.pop('edgecolor')
        plotargs = {'facecolor': activecolor,
                    'edgecolor': edgecolor,
                    'valid': {'facecolor': validcolor}}
        for key in plotargs:
            if key in self.pargs:
                plotargs[key] = self.pargs.pop(key)

        # plot segments
        for flag, label in zip(self.flags, labels)[::-1]:
            if self.state and not self.all_data:
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

register_plot(SegmentDataPlot)


class StateVectorDataPlot(TimeSeriesDataPlot):
    """DataPlot of some `StateVector` data.

    While technically a sub-class of the `TimeSeriesDataPlot`, for
    data access and processing reasons, the output shadows that of the
    `SegmentDataPlot` more closely.
    """
    type = 'statevector'
    defaults = SegmentDataPlot.defaults.copy()

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

    def process(self):
        (plot, ax) = self.init_plot(plot=SegmentPlot)

        # extract plotting arguments
        mask = self.pargs.pop('mask')
        ax.set_insetlabels(self.pargs.pop('insetlabels', True))
        activecolor, validcolor = self.get_segment_color()
        edgecolor = self.pargs.pop('edgecolor')
        plotargs = {'facecolor': activecolor,
                    'edgecolor': edgecolor,
                    'valid': {'facecolor': validcolor}}

        # plot segments
        nflags = 0
        for channel in self.channels[::-1]:
            if self.state and not self.all_data:
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
                stateseries.bits = channel.bits
                if not 'int' in str(stateseries.dtype):
                    stateseries = stateseries.astype('uint32')
                newflags = stateseries.to_dqflags().values()
                if flags is None:
                    flags = newflags
                else:
                    for i, flag in enumerate(newflags):
                        flags[i] += flag
            nflags += len([m for m in channel.bits if m is not None])
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

register_plot(StateVectorDataPlot)


class SpectrumDataPlot(DataPlot):
    """Spectrum plot for a `SummaryTab`
    """
    type = 'spectrum'
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

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1], self.state))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.99)

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
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            data = get_spectrum(str(channel), valid, query=False,
                                format=sdform)

            # anticipate log problems
            if self.pargs['logx']:
                data = [s[1:] for s in data]
            if self.pargs['logy']:
                for sp in data:
                    sp[sp.data == 0] = 1e-100

            if use_percentiles:
                ax.plot_spectrum_mmm(*data, **pargs)
            else:
                pargs.pop('alpha', None)
                ax.plot_spectrum(data[0], **pargs)

            # allow channel data to set parameters
            if hasattr(data[0].channel, 'frequency_range'):
                self.pargs.setdefault('xlim', data[0].channel.frequency_range)
            if (sdform in ['amplitude', 'asd'] and
                    hasattr(data[0].channel, 'asd_range')):
                self.pargs.setdefault('ylim', data[0].channel.asd_range)
            elif hasattr(data[0].channel, 'psd_range'):
                self.pargs.setdefault('ylim', data[0].channel.psd_range)

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
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if len(self.channels) > 1:
            plot.add_legend(ax=ax, **legendargs)
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)

        return self.finalize()

register_plot(SpectrumDataPlot)


class TimeSeriesHistogramPlot(DataPlot):
    """HistogramPlot from a Series
    """
    type = 'histogram'
    defaults = {'ylabel': 'Rate [Hz]',
                'rwidth': 1}

    def init_plot(self, plot=HistogramPlot):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesDataPlot`.
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

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1], self.state))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle)

        # get spectrum format: 'amplitude' or 'power'
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
        if not 'range' in histargs:
            histargs['range'] = ax.common_limits(data)

        # plot
        for label, arr in zip(labels, data):
            if arr.size:
                ax.hist(arr, label=label, **histargs)
            else:
                ax.plot([], label=label)

        # customise plot
        legendargs = self.parse_legend_kwargs()
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if len(self.channels) > 1:
            plot.add_legend(ax=ax, **legendargs)

        # add extra axes and finalise
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        return self.finalize()

register_plot(TimeSeriesHistogramPlot)


class TriggerDataPlot(TimeSeriesDataPlot):
    _threadsafe = False
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
        super(TriggerDataPlot, self).__init__(channels, start, end,
                                                 state=state, outdir=outdir,
                                                 tag=tag, **kwargs)
        self.etg = etg
        self.columns = [self.pargs.pop(c) for c in ('x', 'y', 'color')]

    @property
    def tag(self):
        """Unique identifier for this `TriggerDataPlot`.

        Extends the standard `TimeSeriesDataPlot` tag with the ETG
        and each of the column names.
        """
        try:
            return self._tag
        except AttributeError:
            tag = super(TriggerDataPlot, self).tag
            tag += '_%s' % re_cchar.sub('_', self.etg)
            for column in self.columns:
                if column:
                    tag += '_%s' % re_cchar.sub('_', column)
            return tag.upper()

    @tag.setter
    def tag(self, filetag):
        self._tag = filetag or self.type.upper()

    def finalize(self, outputfile=None):
        if isinstance(self.plot, TimeSeriesPlot):
            return super(TriggerDataPlot, self).finalize(
                       outputfile=outputfile)
        else:
            return super(TimeSeriesDataPlot, self).finalize(
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
            ax.set_epoch(float(self.start))
            ax.set_xlim(float(self.start), float(self.end))

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
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            table = get_triggers(str(channel), self.etg, valid, query=False)
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
        legendargs = self.parse_legend_kwargs(markerscale=4)
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
            plot.add_legend(ax=ax, **legendargs)

        # add state segments
        if isinstance(plot, TimeSeriesPlot) and self.state:
            self.add_state_segments(ax)

        # finalise
        return self.finalize()

register_plot(TriggerDataPlot)


class TriggerTimeSeriesDataPlot(TimeSeriesDataPlot):
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
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            data = get_timeseries(channel, valid, query=False)
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
        legendargs = self.parse_legend_kwargs()
        for key, val in self.pargs.iteritems():
            try:
                getattr(ax, 'set_%s' % key)(val)
            except AttributeError:
                setattr(ax, key, val)
        if len(self.channels) > 1:
            plot.add_legend(ax=ax, **legendargs)

        # add extra axes and finalise
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        if self.state:
            self.add_state_segments(ax)
        return self.finalize()

register_plot(TriggerTimeSeriesDataPlot)


class TriggerHistogramPlot(TimeSeriesHistogramPlot):
    """HistogramPlot from a LIGO_LW Table
    """
    type = 'trigger-histogram'
    _threadsafe = False

    def __init__(self, *args, **kwargs):
        super(TriggerHistogramPlot, self).__init__(*args, **kwargs)
        self.etg = self.pargs.pop('etg')

    def init_plot(self, plot=HistogramPlot):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesDataPlot`.
        """
        self.plot = plot(figsize=[12, 6])
        ax = self.plot.gca()
        return self.plot, ax

    def process(self):
        """Get data and generate the figure.
        """
        # get histogram parameters
        (plot, ax) = self.init_plot()

        column = self.pargs.pop('column')

        # extract histogram arguments
        histargs = self.parse_histogram_kwargs()
        legendargs = self.parse_legend_kwargs()

        # work out labels
        labels = self.pargs.pop('labels', self.channels)
        if isinstance(labels, (unicode, str)):
            labels = labels.split(',')
        labels = map(lambda s: str(s).strip('\n '), labels)

        # add data
        data = []
        for label, channel in zip(labels, self.channels):
            channel = get_channel(channel)
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            table_ = get_triggers(str(channel), self.etg, valid, query=False)
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
            plot.add_legend(ax=ax, **legendargs)

        # add extra axes and finalise
        if not plot.colorbars:
            plot.add_colorbar(ax=ax, visible=False)
        return self.finalize()

register_plot(TriggerHistogramPlot)


class TriggerRateDataPlot(TimeSeriesDataPlot):
    """TimeSeriesDataPlot of trigger rate.
    """
    type = 'trigger-rate'
    _threadsafe = False
    defaults = TimeSeriesDataPlot.defaults.copy()
    defaults.update({'column': None,
                     'ylabel': 'Rate [Hz]'})

    def __init__(self, *args, **kwargs):
        if not 'stride' in kwargs:
            raise ValueError("'stride' must be configured for all rate plots.")
        if 'column' in kwargs and 'bins' not in kwargs:
            raise ValueError("'bins' must be configured for rate plots if "
                             "'column' is given.")
        super(TriggerRateDataPlot, self).__init__(*args, **kwargs)
        self.etg = self.pargs.pop('etg')

    def process(self):
        """Read in all necessary data, and generate the figure.
        """

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
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            table_ = get_triggers(str(channel), self.etg, valid, query=False)
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
        out = super(TriggerRateDataPlot, self).process(outputfile=outputfile)
        self.channels = channels
        return out

register_plot(TriggerRateDataPlot)
