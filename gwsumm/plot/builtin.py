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
import bisect
from itertools import (izip, cycle)

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from matplotlib.pyplot import subplots

from dateutil.relativedelta import relativedelta

from astropy.units import Quantity

from gwpy.spectrum import Spectrum
from gwpy.plotter import *
from gwpy.plotter.tex import label_to_latex
from gwpy.plotter.utils import (color_cycle, marker_cycle)
from gwpy.time import (from_gps, to_gps)

from .. import (globalv, mode, version)
from ..config import NoOptionError
from ..utils import (re_quote, re_cchar, split_channels, get_odc_bitmask)
from ..data import (get_channel, get_timeseries, get_spectrogram, get_spectrum,
                    add_timeseries)
from ..segments import (get_segments, format_padding)
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
        kwargs.setdefault('edgecolor', 'black')
        kwargs.setdefault('facecolor', GREEN)
        kwargs.setdefault('known', {'facecolor': 'red'})
        sax = self.plot.add_state_segments(self.state, ax, plotargs=kwargs)
        sax.tick_params(axis='y', which='major', labelsize=12)
        sax.set_epoch(float(self.start))
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
                'cmap': 'jet',
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
        cmap = self.pargs.pop('cmap')
        ratio = self.ratio

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


class SegmentDataPlot(SegmentLabelSvgMixin, TimeSeriesDataPlot):
    """Segment plot of one or more `DataQualityFlags <DataQualityFlag>`.
    """
    type = 'segments'
    data = 'segments'
    defaults = {'mask': None,
                'color': None,
                'on_is_bad': False,
                'insetlabels': 'inset',
                'edgecolor': 'black',
                'legend-bbox_to_anchor': (1.01, 1.),
                'legend-loc': 'upper left',
                'legend-borderaxespad': 0,
                'legend-fontsize': 12}

    def __init__(self, flags, start, end, state=None, outdir='.', **kwargs):
        padding = kwargs.pop('padding', None)
        super(SegmentDataPlot, self).__init__([], start, end, state=state,
                                                 outdir=outdir, **kwargs)
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
        if isinstance(f, DataQualityFlag):
            self._flags.append(f)
        else:
            self._flags.append(DataQualityFlag(f))

    @property
    def padding(self):
        return OrderedDict((f.name, f.padding) for f in self._flags)

    @padding.setter
    def padding(self, pad):
        for f, p in format_padding(self._flags, pad).iteritems():
            if isinstance(p, (float, int)):
                f.padding = (p, p)
            else:
                f.padding = p

    @property
    def ifos(self):
        """Interferometer set for this `SegmentDataPlot`
        """
        allflags = [f for flag in self.flags for f in flag.split(',')]
        return set([f.strip('!&-_')[:2] for f in allflags])

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
        active = self.pargs.pop('active', None)
        known = self.pargs.pop('known', None)
        # both defined by user
        if active is not None and known is not None:
            return active, known
        # only active defined by user
        elif isinstance(active, str) and active.lower() != 'red':
            return active, 'red'
        elif active is not None:
            return active, 'blue'
        # only known defined by user
        elif known not in [None, GREEN, 'green', 'g']:
            return GREEN, known
        elif known is not None:
            return 'blue', known
        else:
            onisbad = bool(self.pargs.pop('on_is_bad', True))
            if onisbad:
                return 'red', GREEN
            else:
                return GREEN, 'red'

    def process(self):
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
        activecolor, validcolor = self.get_segment_color()
        if isinstance(activecolor, dict):
            self.pargs.update(activecolor)
        elif isinstance(activecolor, tuple):
            self.pargs['facecolor'] = [activecolor] * len(self.flags)
        else:
            self.pargs['facecolor'] = activecolor
        plotargs = self.parse_plot_kwargs()
        for i, kwdict in enumerate(plotargs):
            if isinstance(validcolor, dict):
                kwdict['known'] = validcolor
            elif (validcolor is None or isinstance(validcolor, str) or
                    isinstance(validcolor[0], (float, int))):
                kwdict['known'] = {'facecolor': validcolor}
            else:
                kwdict['known'] = {'facecolor': validcolor[i]}
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
            ax.plot(segs, y=i, label=label, **pargs)

        # make custom legend
        known = legcolors.pop('known', None)
        if known:
            active = legcolors.pop('facecolor')
            edgecolor = legcolors.pop('edgecolor')
            epoch = ax.get_epoch()
            xlim = ax.get_xlim()
            seg = SegmentList([Segment(self.start - 10, self.start - 9)])
            v = ax.plot(seg, facecolor=known['facecolor'],
                        collection=False)[0][0]
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

    def process(self):
        # make font size smaller
        _labelsize = rcParams['ytick.labelsize']
        labelsize = self.pargs.pop('labelsize', 12)
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
        activecolor, validcolor = self.get_segment_color()
        edgecolor = self.pargs.pop('edgecolor')
        plotargs = {'facecolor': activecolor,
                    'edgecolor': edgecolor}
        if isinstance(validcolor, dict):
            plotargs['known'] = validcolor
        elif (validcolor is None or isinstance(validcolor, str) or
                isinstance(validcolor[0], (float, int))):
            plotargs['known'] = {'facecolor': validcolor}
        else:
            plotargs['known'] = {'facecolor': validcolor[i]}
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
            nflags += len([m for m in bits_ if m is not None])
            labels = pargs.pop('label', [None]*len(flags))
            if isinstance(labels, str):
                labels = [labels]
            while len(labels) < len(flags):
                labels.append(None)
            for flag, label in zip(flags, labels)[::-1]:
                kwargs = pargs.copy()
                kwargs.update(plotargs)
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

        # reset tick size and return
        rcParams['ytick.labelsize'] = _labelsize
        return self.finalize()

register_plot(StateVectorDataPlot)


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

    def init_plot(self, plot=HistogramPlot):
        """Initialise the Figure and Axes objects for this
        `TimeSeriesDataPlot`.
        """
        self.plot = plot(figsize=self.pargs.pop('figsize', [12, 6]))
        ax = self.plot.gca()
        ax.grid(True, which='both')
        return self.plot, ax

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
        # get histogram parameters
        (plot, ax) = self.init_plot()

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        label_to_latex(str(self.state))))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle)

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
            l = ax.common_limits(data)
            for d in histargs:
                d['range'] = l

        # plot
        for arr, pargs in zip(data, histargs):
            if arr.size == 0:
                kwargs = dict(
                    (k, pargs[k]) for k in ['label', 'color'] if pargs.get(k))
                ax.plot([], **kwargs)
            else:
                ax.hist(arr, **pargs)

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
        return self.finalize(outputfile=outputfile)

register_plot(TimeSeriesHistogramPlot)


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

    def get_bins(self):
        """Work out the correct histogram binning for this `DutyDataPlot`
        """
        # if not given anything, work it out from the mode
        if self.bins is None:
            m = mode.MODE_NAME[mode.get_mode()]
            # for day mode, make hourly duty factor
            if m in ['DAY']:
                dt = relativedelta(hours=1)
            # for week and month mode, use daily
            elif m in ['WEEK', 'MONTH']:
                dt = relativedelta(days=1)
            # for year mode, use a month
            elif m in ['YEAR']:
                dt = relativedelta(months=1)
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
            bin = SegmentList([Segment(self.start + sum(bins[:i]),
                                       self.start + sum(bins[:i+1]))])
            d = float(abs(segments & bin))
            if normalized:
                d *= normalized / float(bins[i])
            duty[i] = d
            mean[i] = duty[:i+1].mean()
        if cumulative:
            duty = duty.cumsum()
        return duty, mean

    def process(self, outputfile=None):
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

        # work out times and plot mean for legend
        self.get_bins()
        times = float(self.start) + numpy.concatenate(
                                 ([0], self.bins[:-1].cumsum()))
        if not cumulative:
            axes[0].plot(times[:1], [-1], 'k--', label='Rolling mean')

        try:
            bottom = self.pargs['ylim'][0]
        except KeyError:
            bottom = 0

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
            elif 'label' in pargs and normalized == 'percent':
                if legendargs.get('loc', None) in ['upper left', 2]:
                    pargs['label'] = pargs['label'] + '\n[%.1f\\%%]' % mean[-1]
                else:
                    pargs['label'] = pargs['label'] + r' [%.1f\%%]' % mean[-1]
            color = pargs.pop('color', color)
            now = bisect.bisect_left(times, globalv.NOW)
            if sidebyside:
                t = times + self.bins * (0.1 + .8 * i/len(self.flags))
                b = ax.bar(t[:now], (duty-bottom)[:now], bottom=bottom,
                           width=.8 * self.bins[:now]/len(self.flags),
                           color=color, **pargs)
            else:
                b = ax.bar(times[:now], (duty-bottom)[:now], bottom=bottom,
                           width=self.bins[:now], color=color, **pargs)
            # plot mean
            if not cumulative:
                t = [self.start] + list(times + self.bins/2.) + [self.end]
                mean = [mean[0]] + list(mean) + [mean[-1]]
                ax.plot(t, mean, color=sep and 'k' or color, linestyle='--')

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
        if not cumulative:
            yoff = 0.01 * float.__div__(*axes[0].get_position().size)
            lkwargs = legendargs.copy()
            lkwargs['loc'] = 'lower right'
            lkwargs['bbox_to_anchor'] = (1.0, 1. + yoff)
            lkwargs['fontsize'] = 12
            axes[0].add_artist(axes[0].legend(['Rolling mean'], **lkwargs))
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
        'cmap': 'jet',
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
            plot.suptitle(suptitle, y=0.99)

        # parse plotting arguments
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
        plotargs.setdefault('cmap', 'jet')

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
            m = ax.plot_variance(variance, **plotargs)
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
                'cmap': 'spectral',
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

    def process(self):
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


class SegmentPiePlot(SegmentDataPlot):
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

    def parse_plot_kwargs(self, defaults=dict()):
        """Parse pie() keyword arguments
        """
        plotargs = defaults.copy()
        plotargs.setdefault('labels', self._parse_labels())
        for kwarg in ['explode', 'colors', 'autopct', 'pctdistance', 'shadow',
                      'labeldistance', 'startangle', 'radius', 'counterclock',
                      'wedgeprops', 'textprops']:
            try:
                val = self.pargs.pop(kwarg)
            except KeyError:
                try:
                    val = self.pargs.pop('%ss' % kwarg)
                except KeyError:
                    val = None
            if val is not None:
                try:
                    val = eval(val)
                except Exception:
                    pass
                plotargs[kwarg] = val
        return plotargs

    def parse_wedge_kwargs(self, defaults=dict()):
        wedgeargs = defaults.copy()
        for key in self.pargs.keys():
            if key.startswith('wedge-'):
                wedgeargs[key[6:]] = self.pargs.pop(key)
        return wedgeargs

    def process(self):
        (plot, axes) = self.init_plot(plot=Plot)
        ax = axes[0]

        # get labels
        #flags = map(lambda f: str(f).replace('_', r'\_'), self.flags)
        #labels = self.pargs.pop('labels', self.pargs.pop('label', flags))
        #labels = map(lambda s: re_quote.sub('', str(s).strip('\n ')), labels)

        # extract plotting arguments
        legendargs = self.parse_legend_kwargs()
        wedgeargs = self.parse_wedge_kwargs()
        plotargs = self.parse_plot_kwargs()

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
        tot = float(sum(data))
        pclabels = []
        for d, label in zip(data, labels):
            try:
                pc = d/tot * 100
            except ZeroDivisionError:
                pc = 0.0
            pclabels.append(label_to_latex(
                '%s [%1.1f%%]' % (label, pc)).replace(r'\\', '\\'))
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
        ax.set_position([axpos.x0+offset, .05, axpos.width, .9])

        # add bit mask axes and finalise
        self.pargs['xlim'] = None
        return self.finalize(transparent="True", pad_inches=0)

register_plot(SegmentPiePlot)
