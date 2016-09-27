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

from __future__ import print_function

import os
import re
from multiprocessing import cpu_count

from dateutil import tz

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

import numpy

from matplotlib.pyplot import subplots
from matplotlib.ticker import NullLocator

from gwpy.time import Time
from gwpy.timeseries import TimeSeriesDict
from gwpy.plotter import TimeSeriesPlot

from .registry import (get_tab, register_tab)
from .. import (globalv, html)
from ..config import (NoOptionError, GWSummConfigParser)
from ..data import (get_timeseries_dict, get_channel)
from ..plot.registry import get_plot
from ..utils import (vprint, re_quote)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['SEIWatchDogTab']

base = get_tab('default')
HAMs = ['HAM%d' % i for i in range(1, 7)]
BSCs = ['BS', 'ITMX', 'ITMY', 'ETMX', 'ETMY']

HEPI_GPS_CHANNEL = '%s:HPI-%s_WD_MON_GPS_TIME'
HEPI_LATCH_CHANNEL = '%s:HPI-%s_WD_MON_FIRSTTRIG_LATCH'

HAM_ISI_GPS_CHANNEL = '%s:ISI-%s_WD_MON_GPS_TIME'
HAM_ISI_LATCH_CHANNEL = '%s:ISI-%s_WD_MON_FIRSTTRIG_LATCH'
BSC_ST1_GPS_CHANNEL = '%s:ISI-%s_ST1_WD_MON_GPS_TIME'
BSC_ST1_LATCH_CHANNEL = '%s:ISI-%s_ST1_WD_MON_FIRSTTRIG_LATCH'
BSC_ST2_GPS_CHANNEL = '%s:ISI-%s_ST2_WD_MON_GPS_TIME'
BSC_ST2_LATCH_CHANNEL = '%s:ISI-%s_ST2_WD_MON_FIRSTTRIG_LATCH'

re_no_count = re.compile('(ISI (.*)IOP|(.*) Reset|(.*)from stage \d+)')

class SEIWatchDogTab(base):
    """Summarise the WatchDog trips recorded from the SEI system.
    """
    type = 'seismic-watchdog'
    window = 5

    @classmethod
    def from_ini(cls, config, section, plotdir='plots', **kwargs):
        """Define a new `SEIWatchDogTab`.
        """
        chambers = config.get(section, 'chamber-type')
        new = super(SEIWatchDogTab, cls).from_ini(config, section,
                                                  plotdir=plotdir, **kwargs)
        new.ifo = config.get('DEFAULT', 'ifo')
        new.plotdir = plotdir
        chamber = str(chambers).upper()
        if chamber == 'HAM':
            new.chambers = HAMs
        elif chamber == 'BSC':
            new.chambers = BSCs
        else:
            raise ValueError("Cannot process chamber-type = '%s'" % chamber)
        try:
            new.plot_duration = config.getint(section, 'plot-duration')
        except NoOptionError:
            new.plot_duration = 30
        return new

    def process(self, nds=None, multiprocess=True,
                config=GWSummConfigParser(), datacache=None,
                trigcache=None, datafind_error='raise', **kwargs):
        """Process data for the given state.
        """

        ifo = self.ifo

        for p in self.plots:
            if p.outputfile in globalv.WRITTEN_PLOTS:
                p.new = False

        # --------------------------------------------------------------------
        # work out which channels are needed

        tschannels = [HEPI_GPS_CHANNEL % (ifo, chamber) for
                      chamber in self.chambers]
        svchannels = [HEPI_LATCH_CHANNEL % (ifo, chamber) for
                      chamber in self.chambers]
        if self.chambers == HAMs:
            tschannels.extend([HAM_ISI_GPS_CHANNEL % (ifo, chamber) for
                               chamber in self.chambers[1:]])
            svchannels.extend([HAM_ISI_LATCH_CHANNEL % (ifo, chamber) for
                               chamber in self.chambers[1:]])
        else:
            tschannels.extend([BSC_ST1_GPS_CHANNEL % (ifo, chamber) for
                               chamber in self.chambers])
            svchannels.extend([BSC_ST1_LATCH_CHANNEL % (ifo, chamber) for
                               chamber in self.chambers])
            tschannels.extend([BSC_ST2_GPS_CHANNEL % (ifo, chamber) for
                               chamber in self.chambers])
            svchannels.extend([BSC_ST2_LATCH_CHANNEL % (ifo, chamber) for
                               chamber in self.chambers])

        state = sorted(self.states, key=lambda s: abs(s.active))[0]

        vprint("    %d channels identified for TimeSeries\n" % len(tschannels))
        tripdata = get_timeseries_dict(tschannels, state, config=config,
                                       nds=nds, multiprocess=multiprocess,
                                       datafind_error=datafind_error,
                                       cache=datacache)
        vprint("    All time-series data loaded\n")

        vprint("    %d channels identified as StateVectors\n" % len(svchannels))
        latchdata = get_timeseries_dict(svchannels, state, config=config,
                                        nds=nds, multiprocess=multiprocess,
                                        datafind_error=datafind_error,
                                        statevector=True, cache=datacache)
        vprint("    All state-vector data loaded\n")

        # --------------------------------------------------------------------
        # find the trips

        self.trips = []
        for (gpschannel, latch) in zip(tschannels, svchannels):
            # get channel data
            gpschannel = get_channel(gpschannel)
            latch = get_channel(latch)
            chamber = gpschannel.subsystem
            system = gpschannel.system
            vprint("    Locating WD trips for %s %s %s...\n"
                   % (ifo, chamber, system))

            # find times of a trip
            trips = []
            for i, ts in enumerate(tripdata[gpschannel.name]):
                alltrips = ts.times[((numpy.diff(ts.value) > 0) &
                                     (ts.value[1:] > 1e8)).nonzero()[0] + 1]
                for j, gpstime in enumerate(alltrips.value):
                    trips.append((i, gpstime))

            vprint("        Found %d WD trips.\n" % len(trips))

            # associate cause
            for i, trip in enumerate(trips):
                tsid, t = trip
                gpstime = int(t)
                # extract 1 second of LATCH data
                ldata = latchdata[latch.name][tsid].crop(gpstime,
                                                         gpstime + 1).value
                # find transition value
                try:
                    bits = ldata[sorted(numpy.unique(ldata,
                                                     return_index=True)[1])][1]
                except IndexError:
                    bits = ldata[0]
                except (IndexError, ValueError):
                    bits = None
                # associate cause
                if not bits:
                    if re.match('ST\d', latch.signal):
                        stage = 'ISI %s' % latch.signal.split('_')[0]
                    else:
                        stage = system
                    causes = ['%s Unknown' % stage]
                else:
                    allbits = numpy.nonzero(map(int, bin(int(bits))[2:][::-1]))[0]
                    causes = [latch.bits[b] for b in allbits]
                t2 = Time(t, format='gps', scale='utc')
                vprint("        Trip GPS %s (%s), triggers:\n" % (t, t2.iso))
                for cause in causes:
                    vprint("            %s\n" % cause)
                    # configure plot
                    mapsec = 'sei-wd-map-%s' % cause
                    if (not config.has_section(mapsec) and
                            re.match('ISI ST\d ', cause)):
                        mapsec = ('sei-wd-map-%s'
                                  % (' '.join(cause.split(' ', 2)[::2])))
                    if config.has_section(mapsec):
                        pstart = gpstime - self.plot_duration / 2.
                        if self.chambers == HAMs or system == 'HPI':
                            platform = chamber
                        else:
                            platform = '%s_%s' % (
                                chamber, gpschannel.signal.split('_')[0])
                        p = os.path.join(self.plotdir,
                                         '%s-%s_%s_WATCHDOG_TRIP-%d-%d.png'
                                         % (ifo, system, platform, pstart,
                                            self.plot_duration))
                        self.plots.append(SeiWatchDogPlot(
                                              t, chamber, cause, config,
                                              p, ifo=ifo,
                                              nds=nds is True or False))
                        plot = self.plots[-1]
                    else:
                        plot = None
                    self.trips.append((t, chamber, cause, plot))

        super(SEIWatchDogTab, self).process(
            config=config, nds=nds, multiprocess=multiprocess,
            datacache=datacache, trigcache=trigcache, **kwargs)

    def write_state_html(self, state):
        """Build HTML summary of watchdog trips
        """
        # find one example of each channel, and get the bits
        hepichannel = get_channel(HEPI_LATCH_CHANNEL
                                  % (self.ifo, self.chambers[1]))
        hepimask = hepichannel.bits + ['HEPI Unknown']
        if self.chambers == HAMs:
             isichannels = [get_channel(HAM_ISI_LATCH_CHANNEL
                                        % (self.ifo, self.chambers[0]))]
             isimask = isichannels[0].bits + ['ISI Unknown']
        else:
             isichannels = [get_channel(BSC_ST1_LATCH_CHANNEL
                                        % (self.ifo, self.chambers[0])),
                            get_channel(BSC_ST2_LATCH_CHANNEL
                                        % (self.ifo, self.chambers[0]))]
             isimask = (isichannels[0].bits + ['ISI ST1 Unknown'] +
                        isichannels[1].bits + ['ISI ST2 Unknown'])
        mask = hepimask + isimask

        # count trips
        count = {}
        for _, chamber, trigger, _ in self.trips:
            key = (chamber, trigger)
            try:
                count[key] += 1
            except KeyError:
                count[key] = 1

        page = html.markup.page()

        # build summary table
        page.div(class_='well')
        page.table(class_='table table-condensed table-hover watchdog')
        page.thead()
        # add headers
        page.tr(class_='header')
        for th in ['WD'] + self.chambers + ['Sub-total']:
            page.th(th)
        page.tr.close()
        page.thead.close()

        # add rows
        page.tbody()
        totals = numpy.zeros((len(mask), len(self.chambers) + 1),
                             dtype=int)
        rows = []
        for i, bit in enumerate(mask):
            class_ = []
            if (i == len(hepimask) or
                    i == len(hepimask + isichannels[0].bits)):
                class_.append('tdiv')
            if re_no_count.match(bit):
                class_.append('ignore')
            if class_:
                page.tr(class_=' '.join(class_))
            else:
                page.tr()
            page.th(bit)
            for j, chamber in enumerate(self.chambers):
                try:
                    c = count[(chamber, bit)]
                except KeyError:
                    c = 0
                    pass
                page.td(c or '-')
                # exclude IOP from total
                if not re_no_count.match(bit):
                    totals[i][j] = c
            # add row total
            totals[i][-1] = totals[i].sum()
            page.th(totals[i][-1])
            page.tr.close()
        page.tbody.close()

        # add totals
        page.thead()
        page.tr(class_='header')
        page.th('Totals')
        for i in range(totals.shape[1]):
            t = totals[:,i].sum()
            page.th(t)
        page.tr.close()
        page.thead.close()
        page.table.close()
        page.div.close()

        # build trip groups
        self.trips.sort(key=lambda (x, y, z, p):
                                (x, z in mask and (mask).index(z) or 1000,
                                 p is None))
        groups = OrderedDict()
        j = 0
        for i in xrange(len(self.trips)):
            if i == 0:
                j = i
                groups[j] = []
                continue
            t = self.trips[i][0]
            t2 = self.trips[i-1][0]
            if (t - t2) < self.window:
                groups[j].append(i)
            else:
                j = i
                groups[j] = []

        # build trip table
        page.h1('Trip list')
        page.div(class_='well')

        page.p('In the following table, individual watchdog trips are '
               'considered \'associated\' if they fall within %s seconds '
               'of each other.' % self.window)
        utc = tz.gettz('UTC')
        if self.ifo in ['H1', 'C1', 'P1']:
            localzone = tz.gettz('America/Los_Angeles')
        elif self.ifo in ['L1']:
            localzone = tz.gettz('America/Chicago')
        else:
            localzone = tz.gettz('Europe/Berlin')
        headers = ['GPS time', 'UTC time', 'Local time', 'Chamber', 'Trigger',
                   'Plot', 'Associated']
        rows = []
        for id in groups:
            t, chamber, trigger, plot = self.trips[id]
            t2 = Time(t, format='gps', scale='utc')
            tlocal = Time(t2.datetime.replace(tzinfo=utc).astimezone(localzone),
                          format='datetime', scale='utc')
            rows.append([t, t2.iso, tlocal.iso, chamber, trigger])
            if plot:
                rows[-1].append(html.markup.oneliner.a(
                    '[Click here]', href=plot.href, class_='fancybox plot',
                    **{'data-fancybox-group': '1'}))
            else:
                rows[-1].append('-')
            assoc = []
            for id2 in groups[id]:
                t2, chamber2, trigger2, plot2 = self.trips[id2]
                dt = t2 - t
                tag = '%s %s (+%.2fs)' % (chamber2, trigger2, dt)
                if plot2:
                    assoc.append(html.markup.oneliner.a(
                                     tag, href=plot2.href,
                                     class_='fancybox plot',
                                     **{'data-fancybox-group': '1'}))
                else:
                    assoc.append(tag)
            if assoc:
                rows[-1].append('<br>'.join(assoc))
            else:
                rows[-1].append('-')
        page.add(str(html.data_table(headers, rows)))

        wdp = []
        for i, p in enumerate(self.plots):
            if 'WATCHDOG_TRIP' in p.href:
                wdp.append(i)
        for idx in wdp[::-1]:
            self.plots.pop(idx)

        # write trips to data file
        tripfile = os.path.join(self.path, 'trips.dat')
        with open(tripfile, 'w') as f:
            for id in groups:
                t, chamber, cause, _ = self.trips[id]
                if cause in hepimask:
                    stage = 'HEPI'
                elif self.chambers == HAMs:
                    stage = 'ISI'
                elif cause in isichannels[0].bits:
                    stage = 'ISI1'
                else:
                    stage = 'ISI2'
                cause = cause.replace(' ', '_')
                print(t, chamber, stage, cause, file=f)
        page.p()
        page.add('The list of trips can be downloaded ')
        page.a('here.', href=tripfile, alt=os.path.basename(tripfile),
               title='Trip data')
        page.p.close()

        page.div.close()

        # write to file
        idx = self.states.index(state)
        with open(self.frames[idx], 'w') as fobj:
            fobj.write(str(page))
        return self.frames[idx]

register_tab(SEIWatchDogTab)


class SeiWatchDogPlot(get_plot('data')):
    """Plot a specific SEI WatchDog trip
    """
    name = 'watchdog'
    data = 'watchdog'

    def __init__(self, gpstime, chamber, sensor, config, outfile, ifo=None,
                 duration=30, nds=False, multiprocess=False, datacache=None):
        """Configure a new `SeiWatchDogPlot`.
        """
        super(SeiWatchDogPlot, self).__init__([],
                                              int(gpstime) - duration/2.,
                                              int(gpstime) + duration/2.)

        # get params
        if ifo is None:
            ifo = config.get('DEFAULT', 'ifo')
        self.ifo = ifo
        self.chamber = chamber
        self.sensor = sensor
        self.gpstime = gpstime
        self.duration = duration
        self.outputfile = outfile
        self.use_nds = nds
        if isinstance(multiprocess, int):
            self.nproc = multiprocess
        elif multiprocess:
            self.nproc = cpu_count() - 1
        else:
            self.nproc = False

        system = (sensor.split(' ')[0] == 'HEPI' and
                  'HPI' or sensor.split(' ')[0])

        # get channels
        mapsec = 'sei-wd-map-%s' % sensor
        if (not config.has_section(mapsec) and
                re.match('ISI ST\d ', sensor)):
            mapsec = ('sei-wd-map-%s'
                      % (' '.join(sensor.split(' ', 2)[::2])))
        stubs = zip(*sorted([o for o in config.items(mapsec) if o[0].isdigit()],
                             key=lambda (x, y): x))[1]
        if re.search('ISI ST\d ', sensor):
            stage = sensor.split(' ')[1]
            channels = [get_channel('%s:%s-%s_%s_%s'
                                    % (ifo, system, chamber, stage, stub))
                        for stub in stubs]
        else:
            channels = [get_channel('%s:%s-%s_%s' % (ifo, system, chamber, stub))
                        for stub in stubs]

        # set types
        for channel in channels:
            if not hasattr(channel, 'type') or not channel.type:
                channel.ctype = 'adc'

        self.chanlist = channels

        try:
            self.geometry = map(int, config.get(mapsec, 'geometry').split(','))
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

    def process(self):

        # data span
        start = self.gpstime - self.duration / 2.
        end = self.gpstime + self.duration / 2.

        # get data
        if self.use_nds:
            data = TimeSeriesDict.fetch(self.chanlist, start, end)
        else:
            from glue.datafind import GWDataFindHTTPConnection
            conn = GWDataFindHTTPConnection()
            cache = conn.find_frame_urls(self.ifo[0], '%s_C' % self.ifo,
                                         self.start, self.end, urltype='file')
            if len(cache) == 0:
                data = {}
            else:
                data = TimeSeriesDict.read(cache, self.chanlist, start=start,
                                           end=end, nproc=self.nproc)

        # make plot
        plot, axes = subplots(nrows=self.geometry[0], ncols=self.geometry[1],
                              sharex=True,
                              subplot_kw={'projection': 'timeseries'},
                              FigureClass=TimeSeriesPlot, figsize=[12, 6])
        axes[0,0].set_xlim(start, end)
        for channel, ax in zip(self.chanlist, axes.flat):
            ax.set_epoch(self.gpstime)
            # plot data
            try:
                ax.plot(data[channel])
            except KeyError:
                ax.text(self.gpstime, 0.5, "No data", va='center', ha='center',
                        transform=ax.transData)
            # plot trip indicator
            ylim = ax.get_ylim()
            ax.plot([self.gpstime, self.gpstime], ylim, linewidth=0.5,
                    linestyle='--', color='red')
            ax.set_ylim(*ylim)
            ax.set_xlabel('')
            ax.set_title(channel.texname, fontsize=10)
            ax.xaxis.set_minor_locator(NullLocator())
            for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(10)
            for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(16)
        plot.text(0.5, 0.04, 'Time [seconds] from trip (%s)' % self.gpstime,
                  ha='center', va='bottom', fontsize=24)
        plot.text(0.01, 0.5, 'Amplitude %s' % self.unit, ha='left', va='center',
                  rotation='vertical', fontsize=24)

        plot.suptitle('%s %s %s watchdog trip: %s'
                      % (self.ifo, self.chamber, self.sensor, self.gpstime),
                      fontsize=24)

        plot.save(self.outputfile)
        plot.close()
        return self.outputfile
