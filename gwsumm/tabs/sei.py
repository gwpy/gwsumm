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

import os
import re

import numpy

from gwpy.time import Time

from .registry import (get_tab, register_tab)
from .. import (globalv, version, html)
from ..config import GWSummConfigParser
from ..data import (get_timeseries_dict, get_channel)
from ..utils import vprint

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

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


class SEIWatchDogTab(base):
    """Summarise the WatchDog trips recorded from the SEI system.
    """
    type = 'seismic-watchdog'
    window = 5

    @classmethod
    def from_ini(cls, config, section, plotdir='plots', base=''):
        """Define a new `SEIWatchDogTab`.
        """
        chambers = config.get(section, 'chamber-type')
        new = super(SEIWatchDogTab, cls).from_ini(
                  config, section, plotdir=plotdir, base=base)
        new.ifo = config.get('DEFAULT', 'ifo')
        chamber = str(chambers).upper()
        if chamber == 'HAM':
            new.chambers = HAMs
        elif chamber == 'BSC':
            new.chambers = BSCs
        else:
            raise ValueError("Cannot process chamber-type = '%s'" % chamber)
        return new

    def process(self, nds='guess', multiprocess=True,
                config=GWSummConfigParser(), datacache=None,
                trigcache=None):
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
                                       cache=datacache)
        vprint("    All time-series data loaded\n")

        vprint("    %d channels identified as StateVectors\n" % len(svchannels))
        latchdata = get_timeseries_dict(svchannels, state, config=config,
                                        nds=nds, multiprocess=multiprocess,
                                        statevector=True, cache=datacache)
        vprint("    All state-vector data loaded\n")

        # --------------------------------------------------------------------
        # find the trips

        self.trips = []
        for (gpschannel, latch) in zip(tschannels, svchannels):
            # get channel data
            gpschannel = get_channel(gpschannel)
            gpsts = tripdata[gpschannel.name].join(gap='pad', pad=0.0)
            latch = get_channel(latch)
            latchts = latchdata[latch.name].join(gap='pad', pad=0.0)
            chamber = gpschannel.subsystem
            system = gpschannel.system
            vprint("    Locating WD trips for %s %s %s...\n"
                   % (ifo, chamber, system))

            # find times of a trip
            trips = []
            for i, ts in enumerate(tripdata[gpschannel.name]):
                alltrips = ts.times[numpy.nonzero(numpy.diff(ts))[0] + 1]
                for j, gpstime in enumerate(alltrips.data):
                    if j and (gpstime - alltrips.data[j-1]) < self.window:
                        continue
                    trips.append((i, gpstime))

            vprint("        Found %d WD trips.\n" % len(trips))

            # associate cause
            for i, trip in enumerate(trips):
                tsid, t = trip
                gpstime = int(t)
                # extract 1 second of LATCH data
                ldata = latchdata[latch.name][tsid].crop(gpstime,
                                                         gpstime + 1).data
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
                    cause = 'No cause found'
                else:
                    firstbit = map(int, bin(int(bits))[2:][::-1]).index(1)
                    t2 = Time(t, format='gps', scale='utc')
                    vprint("        Trip GPS %d (%s): cause bit %d.\n"
                           % (t, t2.iso, firstbit))
                    cause = latch.bitmask[firstbit]
                self.trips.append((t, chamber, cause, None))

        super(SEIWatchDogTab, self).process(
            config=config, nds=nds, multiprocess=multiprocess,
            datacache=datacache, trigcache=trigcache)

    def build_inner_html(self, state):
        """Build HTML summary of watchdog trips
        """
        # find one example of each channel, and get the bitmask
        hepichannel = get_channel(HEPI_LATCH_CHANNEL
                                  % (self.ifo, self.chambers[1]))
        hepimask = hepichannel.bitmask
        if self.chambers == HAMs:
             isichannels = [get_channel(HAM_ISI_LATCH_CHANNEL
                                        % (self.ifo, self.chambers[0]))]
        else:
             isichannels = [get_channel(BSC_ST1_LATCH_CHANNEL
                                        % (self.ifo, self.chambers[0])),
                            get_channel(BSC_ST2_LATCH_CHANNEL
                                        % (self.ifo, self.chambers[0]))]
        isimask = [bit for c in isichannels for bit in c.bitmask]

        # count trips
        count = {}
        for _, chamber, trigger, _ in self.trips:
            key = (chamber, trigger)
            try:
                count[key] += 1
            except KeyError:
                count[key] = 1

        page = html.markup.page()
        page.div(class_='row')
        page.div(class_='col-md-12')

        # build summary table
        page.h1('Summary of watch-dog trips')
        page.table(class_='watchdog data')
        # add headers
        page.tr(class_='header')
        for th in ['WD'] + self.chambers + ['Sub-total']:
            page.th(th)
        page.tr.close()

        # add rows
        totals = numpy.zeros((len(hepimask+isimask), len(self.chambers) + 1),
                             dtype=int)
        rows = []
        for i, bit in enumerate(hepimask + isimask):
            class_ = []
            if i == len(hepimask):
                class_.append('tdiv')
            if re.match('(ISI (.*)IOP|(.*) Reset|ISI Stage 1)', bit):
                class_.append('IOP')
            if class_:
                page.tr(class_=' '.join(class_))
            else:
                page.tr()
            page.th(bit)
            row = []
            for j, chamber in enumerate(self.chambers):
                try:
                    c = count[(chamber, bit)]
                except KeyError:
                    c = 0
                    pass
                page.td(c or '-')
                # exclude IOP from total
                if not re.match('(ISI (.*)IOP|(.*) Reset|ISI Stage 1)', bit):
                    totals[i][j] = c
            # add row total
            totals[i][-1] = totals[i].sum()
            page.th(totals[i][-1])
            page.tr.close()
        # add totals
        page.tr(class_='header')
        page.th('Totals')
        for i in range(totals.shape[1]):
            t = totals[:,i].sum()
            page.th(t)
        page.tr.close()
        page.table.close()

        # build trip table
        page.h1('Trip list')
        self.trips.sort(key=lambda (x, y, z, p): x)
        headers = ['GPS time', 'UTC time', 'Chamber', 'Trigger', 'Plot']
        rows = []
        for (t, chamber, trigger, plot) in self.trips:
            rows.append([t, Time(t, format='gps', scale='utc').iso,
                         chamber, trigger])
            if plot:
                row.append(html.markup.oneliner.a('[Click here]', href=plot.href,
                                             class_='fancybox plot',
                                             **{'data-fancybox-group': '1'}))
            else:
                rows[-1].append('-')
        page.add(str(html.data_table(headers, rows, table='data')))

        # close tabs
        page.div.close()
        page.div.close()

        # call super just in case
        page.add(str(super(SEIWatchDogTab, self).build_inner_html(state)))

        return page

register_tab(SEIWatchDogTab)
