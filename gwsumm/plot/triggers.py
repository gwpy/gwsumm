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

"""Definitions for event trigger plots
"""

from __future__ import division

import re
from collections import OrderedDict
from itertools import cycle

from six import string_types

from numpy import isinf

from astropy.units import Quantity

from gwpy.detector import (Channel, ChannelList)
from gwpy.segments import SegmentList
from gwpy.plot.gps import GPSTransform
from gwpy.plot.utils import (color_cycle, marker_cycle)

from .. import globalv
from ..utils import re_cchar
from ..data import (get_channel, get_timeseries, add_timeseries)
from ..triggers import (get_triggers, get_time_column)
from .registry import (get_plot, register_plot)
from .utils import (get_column_string, hash, usetex_tex)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

TimeSeriesDataPlot = get_plot('timeseries')

LATEX_OPERATOR = {
    '>=': r'\geq',
    '<=': r'\leq',
}


class TriggerPlotMixin(object):
    """Mixin to overwrite `channels` property for trigger plots

    We don't need to get channel data for trigger plots.
    """
    def __init__(self, *args, **kwargs):
        self.filterstr = kwargs.pop('filterstr', None)
        super(TriggerPlotMixin, self).__init__(*args, **kwargs)

    @property
    def allchannels(self):
        """List of all unique channels for this plot
        """
        chans = set([re.split(r'[#@]', str(c), 1)[0] for c in self._channels])
        return ChannelList(map(Channel, chans))

    @property
    def pid(self):
        try:
            return self._pid
        except AttributeError:
            chans = "".join(map(str, self.channels))
            filts = "".join(map(str, [
                getattr(c, 'filter', getattr(c, 'frequency_response', ''))
                for c in self.channels]))
            if self.filterstr:
                filts += self.filterstr
            self._pid = hash(chans + filts)
            return self.pid


class TriggerDataPlot(TriggerPlotMixin, TimeSeriesDataPlot):
    """Standard event trigger plot
    """
    type = 'triggers'
    data = 'triggers'
    defaults = TimeSeriesDataPlot.defaults.copy()
    defaults.update({
        'x': 'time',
        'y': 'snr',
        'color': None,
        'edgecolor': 'face',
        'facecolor': None,
        'marker': 'o',
        's': 20,
        'vmin': None,
        'vmax': None,
        'clim': None,
        'cmap': 'YlGnBu',
        'logcolor': False,
        'colorlabel': None,
    })

    def __init__(self, channels, start, end, state=None, outdir='.',
                 etg=None, **kwargs):
        super(TriggerDataPlot, self).__init__(channels, start, end,
                                              state=state, outdir=outdir,
                                              **kwargs)
        self.etg = etg
        self.columns = [self.pargs.pop(c) for c in ('x', 'y', 'color')]

    @property
    def pid(self):
        """Unique identifier for this `TriggerDataPlot`.

        Extends the standard `TimeSeriesDataPlot` pid with the ETG
        and each of the column names.
        """
        try:
            return self._pid
        except AttributeError:
            super(TriggerDataPlot, self).pid
            self._pid += '_%s' % re_cchar.sub('_', self.etg)
            for column in self.columns:
                if column:
                    self._pid += '_%s' % re_cchar.sub('_', column)
            self._pid = self._pid.upper()
            return self.pid

    @pid.setter
    def pid(self, id_):
        self._pid = str(id_)

    def draw(self):
        # get columns
        xcolumn, ycolumn, ccolumn = self.columns

        # initialise figure
        plot = self.init_plot()
        ax = plot.gca()
        ax.grid(True, which='both')

        # work out labels
        labels = self.pargs.pop('labels', self.channels)
        if isinstance(labels, string_types):
            labels = labels.split(',')
        labels = [str(s).strip('\n ') for s in labels]

        # get colouring params
        cmap = self.pargs.pop('cmap')
        clim = self.pargs.pop('clim', self.pargs.pop('colorlim', None))
        cnorm = 'log' if self.pargs.pop('logcolor', False) else None
        clabel = self.pargs.pop('colorlabel', None)
        no_loudest = self.pargs.pop('no-loudest', False) is not False
        loudest_by = self.pargs.pop('loudest-by', None)

        # get plot arguments
        plotargs = []
        for i in range(len(self.channels)):
            plotargs.append(dict())
        # get plot arguments
        for key in ['vmin', 'vmax', 'edgecolor', 'facecolor', 'cmap', 's',
                    'marker', 'rasterized']:
            try:
                val = self.pargs.pop(key)
            except KeyError:
                continue
            if key == 'facecolor' and len(self.channels) > 1 and val is None:
                val = color_cycle()
            if key == 'marker' and len(self.channels) > 1 and val is None:
                val = marker_cycle()
            elif isinstance(val, (list, tuple, cycle)):
                val = cycle(val)
            else:
                val = cycle([val] * len(self.channels))
            for i in range(len(self.channels)):
                plotargs[i][key] = next(val)

        # add data
        valid = SegmentList([self.span])
        if self.state and not self.all_data:
            valid &= self.state.active
        ntrigs = 0
        for channel, label, pargs in zip(self.channels, labels, plotargs):
            try:
                channel = get_channel(channel)
            except ValueError:
                pass
            if '#' in str(channel) or '@' in str(channel):
                key = '%s,%s' % (str(channel),
                                 self.state and str(self.state) or 'All')
            else:
                key = str(channel)
            table = get_triggers(key, self.etg, valid, query=False)
            if self.filterstr is not None:
                table = table.filter(self.filterstr)
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
                if (getattr(channel, param, None) is not None
                        and c in ('x', 'y')):
                    self.pargs.setdefault(lim, getattr(channel, param))
                    if isinstance(self.pargs[lim], Quantity):
                        self.pargs[lim] = self.pargs[lim].value
                # set clim separately
                elif hasattr(channel, param):
                    if not clim:
                        clim = getattr(channel, param)

            ax.scatter(table[xcolumn], table[ycolumn],
                       c=table[ccolumn] if ccolumn else None,
                       label=label, **pargs)

        # customise plot
        legendargs = self.parse_legend_kwargs(markerscale=3)
        if len(self.channels) == 1:
            self.pargs.setdefault('title', usetex_tex(
                '%s (%s)' % (str(self.channels[0]), self.etg)))
        for axis in ('x', 'y'):  # prevent zeros on log scale
            scale = getattr(ax, 'get_{0}scale'.format(axis))()
            lim = getattr(ax, 'get_{0}lim'.format(axis))()
            if scale == 'log' and lim[0] <= 0 and not ntrigs:
                getattr(ax, 'set_{0}lim'.format(axis))(1, 10)

        self.apply_parameters(ax, **self.pargs)

        # correct log-scale empty axes
        if any(map(isinf, ax.get_ylim())):
            ax.set_ylim(0.1, 10)

        # add colorbar
        if ccolumn:
            if not ntrigs:
                ax.scatter([1], [1], c=[1], visible=False)
            ax.colorbar(cmap=cmap, clim=clim, norm=cnorm, label=clabel)

        if len(self.channels) == 1 and len(table) and not no_loudest:
            columns = [x for x in
                       (loudest_by or ccolumn or ycolumn, xcolumn, ycolumn,
                        ccolumn) if x is not None]
            self.add_loudest_event(ax, table, *columns, fontsize='large')

        if len(self.channels) > 1:
            ax.legend(**legendargs)

        # add state segments
        if isinstance(ax.xaxis.get_transform(), GPSTransform):
            self.add_state_segments(ax)
            self.add_future_shade()

        # finalise
        return self.finalize()

    def add_loudest_event(self, ax, table, rank, *columns, **kwargs):
        # get loudest row
        idx = table[rank].argmax()
        row = table[idx]
        x = float(row[columns[0]])
        y = float(row[columns[1]])

        # clip loudest event to axes limits
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        x1 = max(min(x, xlim[1]), xlim[0])
        y1 = max(min(y, ylim[1]), ylim[0])
        if x1 != x or y1 != y:  # loudest event is out of view
            facecolor = 'pink'
            clipon = False
        else:
            facecolor = 'gold'
            clipon = True

        # mark loudest row with star
        coll = ax.scatter([x1], [y1], marker='*', zorder=1000,
                          facecolor=facecolor, edgecolor='black',
                          s=200, clip_on=clipon)

        # get text
        txt = []
        for col in OrderedDict.fromkeys((rank,) + columns):  # unique ordered
            # format column name
            colstr = get_column_string(col)
            # format row value
            try:
                valstr = '{0:.2f}'.format(row[col]).rstrip('.0')
            except ValueError:  # not float()able
                valstr = str(row[col])
            txt.append('{col} = {val}'.format(col=colstr, val=valstr))

        # get position for new text
        try:
            pos = kwargs.pop('position')
        except KeyError:  # user didn't specify, set default and shunt title
            pos = [0.5, 1.00]
            tpos = ax.title.get_position()
            ax.title.set_position((tpos[0], tpos[1] + 0.05))

        # parse text kwargs
        text_kw = {  # defaults
            'transform': ax.transAxes,
            'verticalalignment': 'bottom',
            'horizontalalignment': 'center',
        }
        text_kw.update(kwargs)
        if 'ha' in text_kw:  # handle short versions or alignment params
            text_kw['horizontalalignment'] = text_kw.pop('ha')
        if 'va' in text_kw:
            text_kw['verticalalignment'] = text_kw.pop('va')

        # add text
        text = ax.text(pos[0], pos[1],
                       'Loudest event: {0}'.format(', '.join(txt)),
                       **text_kw)

        return coll, text


register_plot(TriggerDataPlot)


class TriggerTimeSeriesDataPlot(TimeSeriesDataPlot):
    """Custom time-series plot to handle discontiguous `TimeSeries`.
    """
    type = 'trigger-timeseries'
    data = 'triggers'

    def draw(self):
        """Read in all necessary data, and generate the figure.
        """
        plot = self.init_plot()
        ax = plot.gca()

        # work out labels
        labels = self.pargs.pop('labels', self.channels)
        if isinstance(labels, string_types):
            labels = labels.split(',')
        labels = [str(s).strip('\n ') for s in labels]

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
                if self.logy:
                    ts.value[ts.value == 0] = 1e-100
                if color is None:
                    line = ax.plot(ts, label=label)[0]
                    color = line.get_color()
                else:
                    ax.plot(ts, color=color, label=None)

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
        self.apply_parameters(ax, **self.pargs)
        if len(self.channels) > 1:
            ax.legend(**legendargs)

        # finalise
        self.add_state_segments(ax)
        return self.finalize()


register_plot(TriggerTimeSeriesDataPlot)


class TriggerHistogramPlot(TriggerPlotMixin, get_plot('histogram')):
    """HistogramPlot from a LIGO_LW Table
    """
    type = 'trigger-histogram'
    data = 'triggers'

    def __init__(self, *args, **kwargs):
        super(TriggerHistogramPlot, self).__init__(*args, **kwargs)
        self.etg = self.pargs.pop('etg')
        self.column = self.pargs.pop('column')

    @property
    def pid(self):
        try:
            return self._pid
        except AttributeError:
            etg = re_cchar.sub('_', self.etg).upper()
            self._pid = '%s_%s' % (etg, super(TriggerHistogramPlot, self).pid)
            if self.column:
                self._pid += '_%s' % re_cchar.sub('_', self.column).upper()
            return self.pid

    @pid.setter
    def pid(self, id_):
        self._pid = str(id_)

    def draw(self):
        """Get data and generate the figure.
        """
        # get histogram parameters
        plot = self.init_plot()
        ax = plot.gca()

        # extract histogram arguments
        histargs = self.parse_plot_kwargs()
        legendargs = self.parse_legend_kwargs()

        # add data
        data = []
        livetime = []
        for channel in self.channels:
            try:
                channel = get_channel(channel)
            except ValueError:
                pass
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            if '#' in str(channel) or '@' in str(channel):
                key = '%s,%s' % (str(channel),
                                 self.state and str(self.state) or 'All')
            else:
                key = str(channel)
            table_ = get_triggers(key, self.etg, valid, query=False)
            if self.filterstr is not None:
                table_ = table_.filter(self.filterstr)
            livetime.append(float(abs(table_.meta['segments'])))
            data.append(table_[self.column])
            # allow channel data to set parameters
            if hasattr(channel, 'amplitude_range'):
                self.pargs.setdefault('xlim', channel.amplitude_range)

        # plot
        for arr, d, pargs in zip(data, livetime, histargs):
            # set range if not given
            if pargs.get('range') is None:
                pargs['range'] = self._get_range(
                    d,
                    # use range from first dataset if already calculated
                    range=histargs[0].get('range'),
                    # use xlim if manually set (user or INI)
                    xlim=None if ax.get_autoscalex_on() else ax.get_xlim(),
                )
            pargs.setdefault('label', None)
            if pargs.get('log', True):
                pargs.setdefault('bottom', 1e-200)
            ax.hist(arr, **pargs)

        # tight scale the axes
        try:
            d = pargs.pop('orientation', 'vertical')
        except NameError:
            pass
        else:
            if d == 'vertical':
                ax.autoscale_view(tight=True, scaley=False)
            elif d == 'horizontal':
                ax.autoscale_view(tight=True, scalex=False)

        # customise plot
        self.apply_parameters(ax, **self.pargs)
        if len(self.channels) > 1:
            ax.legend(**legendargs)

        # finalise
        return self.finalize()


register_plot(TriggerHistogramPlot)


class TriggerRateDataPlot(TriggerPlotMixin, TimeSeriesDataPlot):
    """TimeSeriesDataPlot of trigger rate.
    """
    type = 'trigger-rate'
    data = 'triggers'
    defaults = TimeSeriesDataPlot.defaults.copy()
    defaults.update({
        'column': None,
        'legend-bbox_to_anchor': (1., 1.),
        'legend-loc': 'upper left',
        'legend-markerscale': 3,
        'legend-frameon': False,
        'ylabel': 'Rate [Hz]',
    })

    def __init__(self, *args, **kwargs):
        if 'stride' not in kwargs:
            raise ValueError("'stride' must be configured for all rate plots.")
        if 'column' in kwargs and 'bins' not in kwargs:
            raise ValueError("'bins' must be configured for rate plots if "
                             "'column' is given.")
        super(TriggerRateDataPlot, self).__init__(*args, **kwargs)
        self.etg = self.pargs.pop('etg')
        self.column = self.pargs.pop('column')

    @property
    def pid(self):
        try:
            return self._pid
        except AttributeError:
            etg = re_cchar.sub('_', self.etg).upper()
            pid = '%s_%s' % (etg, super(TriggerRateDataPlot, self).pid)
            if self.column:
                self.pid += '_%s' % re_cchar.sub('_', self.column).upper()
            return pid

    @pid.setter
    def pid(self, id_):
        self._pid = str(id_)

    def draw(self):
        """Read in all necessary data, and generate the figure.
        """

        # get rate arguments
        stride = self.pargs.pop('stride')
        if self.column:
            cname = get_column_string(self.column)
            bins = self.pargs.pop('bins')
            operator = self.pargs.pop('operator', '>=')
            try:
                opstr = LATEX_OPERATOR[operator]
            except KeyError:
                opstr = str(operator)
        else:
            bins = ['_']

        # work out labels
        labels = self.pargs.pop('labels', None)
        if isinstance(labels, string_types):
            labels = labels.split(',')
        elif labels is None and self.column and len(self.channels) > 1:
            labels = []
            for channel, bin in [(c, b) for c in self.channels for b in bins]:
                labels.append(r' '.join([channel, '$%s$' % opstr,
                                         str(b)]))
            self.pargs.setdefault('legend-title', cname)
        elif labels is None and self.column:
            labels = [r' '.join(['$%s$' % opstr, str(b)]) for b in bins]
            self.pargs.setdefault('legend-title', cname)
        elif labels is None:
            labels = self.channels
        self.pargs['labels'] = [str(s).strip('\n ') for s in labels]

        # get time column
        tcol = self.pargs.pop('timecolumn', None)

        # generate data
        keys = []
        for channel in self.channels:
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            if '#' in str(channel) or '@' in str(channel):
                key = '%s,%s' % (str(channel),
                                 str(self.state) if self.state else 'All')
            else:
                key = str(channel)
            table_ = get_triggers(key, self.etg, valid, query=False)
            if self.filterstr is not None:
                table_ = table_.filter(self.filterstr)
            tcol_ = tcol or get_time_column(table_, self.etg)
            if self.column:
                rates = list(table_.binned_event_rates(
                    stride, self.column, bins, operator=operator,
                    start=self.start, end=self.end, timecolumn=tcol_).values())
            else:
                rates = [table_.event_rate(stride, start=self.start,
                                           end=self.end, timecolumn=tcol_)]
            for bin, rate in zip(bins, rates):
                rate.channel = channel
                keys.append('%s_%s_EVENT_RATE_%s_%s'
                            % (str(channel), str(self.etg),
                               str(self.column), bin))
                if keys[-1] not in globalv.DATA:
                    add_timeseries(rate, keys[-1])

        # reset channel lists and generate time-series plot
        channels = self.channels
        outputfile = self.outputfile
        self.channels = keys
        out = super(TriggerRateDataPlot, self).draw(outputfile=outputfile)
        self.channels = channels
        return out


register_plot(TriggerRateDataPlot)
