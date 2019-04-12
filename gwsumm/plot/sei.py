# coding=utf-8
# Copyright (C) Duncan Macleod (2013)
#
# This file is part of GWSumm
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
# along with GWSumm.  If not, see <http://www.gnu.org/licenses/>

"""`SummaryTab` for seismic watchdog monitoring
"""

import re
from configparser import NoOptionError

from matplotlib.pyplot import subplots
from matplotlib.ticker import NullLocator

from gwpy.plot import Plot
from gwpy.timeseries import TimeSeriesDict

from ..channels import get_channel
from ..utils import re_quote
from .registry import (get_plot, register_plot)


class SeiWatchDogPlot(get_plot('data')):
    """Plot a specific SEI WatchDog trip
    """
    type = 'watchdog'
    data = 'watchdog'

    def __init__(self, gpstime, chamber, sensor, config, outfile, ifo=None,
                 duration=30, nds=False, datacache=None):
        """Configure a new `SeiWatchDogPlot`.
        """
        super(SeiWatchDogPlot, self).__init__([],
                                              int(gpstime) - duration/2.,
                                              int(gpstime) + duration/2.)

        # get params
        if ifo is None:
            ifo = config.get('DEFAULT', 'IFO')
        self.ifo = ifo
        self.chamber = chamber
        self.sensor = sensor
        self.gpstime = gpstime
        self.duration = duration
        self.outputfile = outfile
        self.use_nds = nds

        system = (sensor.split(' ')[0] == 'HEPI' and
                  'HPI' or sensor.split(' ')[0])

        # get channels
        mapsec = 'sei-wd-map-%s' % sensor
        if not config.has_section(mapsec) and re.match(r'ISI ST\d ', sensor):
            mapsec = ('sei-wd-map-%s'
                      % (' '.join(sensor.split(' ', 2)[::2])))
        stubs = list(zip(*sorted(
            [o for o in config.items(mapsec) if o[0].isdigit()],
            key=lambda x: x[0],
        )))[1]
        if re.search(r'ISI ST\d ', sensor):
            stage = sensor.split(' ')[1]
            channels = [get_channel('%s:%s-%s_%s_%s'
                                    % (ifo, system, chamber, stage, stub))
                        for stub in stubs]
        else:
            channels = [get_channel('%s:%s-%s_%s'
                                    % (ifo, system, chamber, stub))
                        for stub in stubs]

        # set types
        for channel in channels:
            if not hasattr(channel, 'type') or not channel.type:
                channel.ctype = 'adc'

        self.chanlist = channels

        try:
            self.geometry = list(map(int, config.get(mapsec, 'geometry').split(',')))
        except NoOptionError:
            self.geometry = (len(channels), 1)
        if len(self.chanlist) != self.geometry[0] * self.geometry[1]:
            raise ValueError("Geometry does not match number of channels.")

        try:
            self.unit = '[%s]' % re_quote.sub('', config.get(mapsec, 'unit'))
        except NoOptionError:
            self.unit = ''

    @property
    def outputfile(self):
        return self._outputfile

    @outputfile.setter
    def outputfile(self, outfile):
        self._outputfile = outfile

    def draw(self):

        # data span
        start = self.gpstime - self.duration / 2.
        end = self.gpstime + self.duration / 2.

        # get data
        if self.use_nds:
            data = TimeSeriesDict.fetch(self.chanlist, start, end)
        else:
            from glue.datafind import GWDataFindHTTPConnection
            conn = GWDataFindHTTPConnection()
            cache = conn.find_frame_urls(self.ifo[0], '%s_R' % self.ifo,
                                         self.start, self.end, urltype='file')
            if len(cache) == 0:
                data = {}
            else:
                data = TimeSeriesDict.read(cache, self.chanlist, start=start,
                                           end=end)

        # make plot
        plot, axes = subplots(nrows=self.geometry[0], ncols=self.geometry[1],
                              sharex=True, subplot_kw={'xscale': 'auto-gps'},
                              FigureClass=Plot, figsize=[12, 6])
        axes[0, 0].set_xlim(start, end)
        for channel, ax in zip(self.chanlist, axes.flat):
            ax.set_epoch(self.gpstime)
            # plot data
            try:
                ax.plot(data[channel])
            except KeyError:
                ax.text(self.gpstime, 0.5, "No data", va='center', ha='center',
                        transform=ax.transData)
            # plot trip indicator
            ax.axvline(self.gpstime, linewidth=0.5, linestyle='--',
                       color='red')
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_title(channel.texname, fontsize=10)
            ax.xaxis.set_minor_locator(NullLocator())
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(10)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
        plot.text(0.5, 0.02, 'Time [seconds] from trip (%s)' % self.gpstime,
                  ha='center', va='bottom', fontsize=24)
        plot.text(0.01, 0.5, 'Amplitude %s' % self.unit, ha='left',
                  va='center', rotation='vertical', fontsize=24)

        plot.suptitle('%s %s %s watchdog trip: %s'
                      % (self.ifo, self.chamber, self.sensor, self.gpstime),
                      fontsize=24)

        plot.save(self.outputfile)
        plot.close()
        return self.outputfile


register_plot(SeiWatchDogPlot)
