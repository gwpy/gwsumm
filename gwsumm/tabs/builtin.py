# -*- coding: utf-8 -*-
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

"""Definition of a summary SummaryTab.
"""

import re
import warnings
from StringIO import StringIO
from multiprocessing import Process

from numpy import isclose

from glue.segmentsUtils import tosegwizard

from .registry import (get_tab, register_tab)
from ..plot import (PlotList, registry as plotregistry)
from .. import globalv
from ..data import (get_channel, get_timeseries_dict, get_spectrogram,
                    get_spectrum)
from ..segments import get_segments
from ..triggers import get_triggers
from ..utils import *
from ..config import *
from .. import html

from gwsumm import version
__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

Tab = get_tab('basic')
StateTab = get_tab('state')
re_channel = re.compile('[A-Z]\d:[A-Z]+-[A-Z0-9_]+\Z')


class SimpleStateTab(StateTab):
    """A simple `StateTab` with plots and data source summaries.

    This is the 'default' `Tab` that all configuration tab- sections
    will be formatted for unless otherwise specified.
    """
    type = 'default'

    def __init__(self, name, longname=None, parent=None, children=list(),
                 base='', span=None, states=None, layout=None):
        """Initialise a new :class:`SimpleStateTab`
        """
        super(SimpleStateTab, self).__init__(name, longname=longname,
                                             parent=parent, children=children,
                                             base=base, span=span,
                                             states=states)
        self.plots = PlotList()
        self.layout = layout

    @property
    def layout(self):
        """List of how many plots to display on each row in the output.

        By default this is ``1`` if the tab contains only 1 or 3 plots,
        or ``2`` if otherwise.
        The final number given in the list will be repeated as necessary.

        :type: `list` of `ints <int>`
        """
        return self._layout

    @layout.setter
    def layout(self, l):
        if isinstance(l, (str, unicode)):
            l = eval(l)
        if l is None:
            self._layout = None
        else:
            self._layout = map(int, l)

    @property
    def channels(self):
        """Set of all data channels used by this tab
        """
        TimeSeriesPlot = plotregistry.get_plot('timeseries')
        SpectrogramPlot = plotregistry.get_plot('spectrogram')
        SpectrumPlot = plotregistry.get_plot('spectrum')
        StateVectorPlot = plotregistry.get_plot('statevector')
        out = set()
        for plot in self.plots:
            if isinstance(plot, (TimeSeriesPlot, SpectrogramPlot,
                                 SpectrumPlot, StateVectorPlot)):
                out.update(plot.channels)
        return sorted(out, key=lambda ch: ch.name)

    @property
    def timeseries(self):
        TimeSeriesPlot = plotregistry.get_plot('timeseries')
        SpectrogramPlot = plotregistry.get_plot('spectrogram')
        SpectrumPlot = plotregistry.get_plot('spectrum')
        out = set()
        for plot in self.plots:
            if isinstance(plot, (TimeSeriesPlot, SpectrogramPlot,
                                 SpectrumPlot)):
                out.update(plot.channels)
        return sorted(out, key=lambda ch: ch.name)

    @property
    def statevectors(self):
        StateVectorPlot = plotregistry.get_plot('statevector')
        out = set()
        for plot in self.plots:
            if isinstance(plot, StateVectorPlot):
                out.update(plot.channels)
        return sorted(out, key=lambda ch: ch.name)

    @property
    def spectra(self):
        SpectrumPlot = plotregistry.get_plot('spectrum')
        out = set()
        for plot in self.plots:
            if isinstance(plot, SpectrumPlot):
                out.update(plot.channels)
        return sorted(out, key=lambda ch: ch.name)

    @property
    def spectrograms(self):
        SpectrogramPlot = plotregistry.get_plot('spectrogram')
        SpectrumPlot = plotregistry.get_plot('spectrum')
        out = set()
        for plot in self.plots:
            if isinstance(plot, (SpectrogramPlot, SpectrumPlot)):
                out.update(plot.channels)
        return sorted(out, key=lambda ch: ch.name)

    @property
    def triggers(self):
        """Set of all trigger channels used by this tab
        """
        TriggerPlot = plotregistry.get_plot('triggers')
        out = set()
        for plot in self.plots:
            if isinstance(plot, TriggerPlot):
                out.update(plot.channels)
        return out

    @property
    def dataqualityflags(self):
        """Set of all data-quality flags used by this tab.

        This does not include those used in state information only.
        """
        dqflags = set()
        re_flag = re.compile('[&!-,]')
        for plot in self.plots:
            if hasattr(plot, 'flags'):
                pflags = [f for pflag in plot.flags for
                          f in re_flag.split(pflag)]
                dqflags.update(pflags)
        return dqflags

    # -------------------------------------------
    # SummaryTab configuration parser

    @classmethod
    def from_ini(cls, cp, section, plotdir='plots', **kwargs):
        """Define a new `SummaryTab` from the given section of the
        `ConfigParser`.

        Parameters
        ----------
        cp : :class:`~gwsumm.config.GWConfigParser`
            customised configuration parser containing given section
        section : `str`
            name of section to parse
        plotdir : `str`, optional, default: ``'plots'``
            output path for plots, relative to current directory

        Returns
        -------
        tab : :class:`SummaryTab`
            a new tab defined from the configuration
        """
        # get [start, stop) job interval
        start = cp.getint('general', 'gps-start-time')
        end = cp.getint('general', 'gps-end-time')
        # get tab name
        if cp.has_option(section, 'name'):
            # name given explicitly
            name = cp.get(section, 'name')
        else:
            # otherwise strip 'tab-' from section name
            name = section[4:]
        if cp.has_option(section, 'longname'):
            longname = cp.get(section, 'longname')
        else:
            longname = name
        # get parent:
        #     if parent is not given, this assumes a top-level tab
        if cp.has_option(section, 'parent'):
            parent = cp.get(section, 'parent')
            if parent == 'None':
                parent = None
        else:
            parent = None
        # parse states and retrieve their definitions
        if cp.has_option(section, 'states'):
            # states listed individually
            statenames = [re_quote.sub('', s) for s in
                          cp.get(section, 'states').split(',')]
        else:
            # otherwise use 'all' state - full span with no gaps
            statenames = ['All']
        states = [globalv.STATES[s] for s in statenames]

        # get layout
        if cp.has_option(section, 'layout'):
            try:
                layout = eval(cp.get(section, 'layout'))
            except NameError:
                raise ValueError("Cannot parse 'layout' for '%s' tab. Layout "
                                 "should be given as a comma-separated list "
                                 "of integers")
            if isinstance(layout, int):
                layout = [layout]
            for l in layout:
                if not l in [1, 2, 3, 4, 6, 12]:
                    raise ValueError("Cannot print more than %d plots in a "
                                     "single row. The chosen layout value for "
                                     "each row must be a divisor of 12 to fit "
                                     "the Bootstrap scaffolding. For details "
                                     "see http://getbootstrap.com/2.3.2/"
                                     "scaffolding.html")
        else:
            layout = None

        # define new job
        job = cls(name, parent=parent, states=states, span=[start, end],
                  longname=longname, layout=layout, **kwargs)
        job._config = cp._sections[section]

        # -------------------
        # parse plot requests
        #    All config entries whose key is a single integer is
        #    interpreted as a requested plot.

        # find and order the plots
        requests = sorted([(int(opt), val) for (opt, val) in
                           cp.nditems(section) if opt.isdigit()],
                          key=lambda a: a[0])

        # parse plot definition
        for index, definition in requests:
            # find plot customisations within this section
            mods = {}
            for key, val in cp.nditems(section):
                if key.startswith('%d-' % index):
                    opt = re_cchar.sub('_', key.split('-', 1)[1].lower())
                    try:
                        mods[opt] = eval(val)
                    except (NameError, SyntaxError):
                        mods[opt] = val

            # parse definition for section references
            try:
                pdef, sources = [s[::-1] for s in
                                 re.split('[\s,]', definition[::-1], 1)]
            except ValueError:
                pdef = definition
                sources = []
            else:
                if not re_channel.match(sources) and cp.has_section(sources):
                    sources = cp.get(sources, 'channels')
                sources = split_channels(sources)

            # if pdef refers to another config section, it must have a type
            if cp.has_section(pdef):
                type_ = cp.get(pdef, 'type')
                PlotClass = plotregistry.get_plot(type_)
            else:
                type_ = None
                PlotClass = plotregistry.get_plot(pdef)
            # if the plot definition declares multiple states
            if 'all_states' in mods:
                if type_:
                    plot = PlotClass.from_ini(cp, pdef, start, end, sources,
                                              state=None, outdir=plotdir,
                                              **mods)
                else:
                    plot = PlotClass(sources, start, end, state=None,
                                     outdir=plotdir, **mods)
                job.plots.append(plot)
            # otherwise define individually for multiple states
            else:
                for state in job.states:
                    if type_:
                        plot = PlotClass.from_ini(cp, pdef, start, end,
                                                  sources, state=state,
                                                  outdir=plotdir, **mods)
                    else:
                        plot = PlotClass(sources, start, end, state=state,
                                         outdir=plotdir, **mods)
                    job.plots.append(plot)

        return job

    # -------------------------------------------
    # SummaryTab processing

    def process_state(self, state, nds='guess', multiprocess=True,
                      config=GWSummConfigParser()):
        """Process data for this tab in a given state
        """
        vprint("Processing '%s' state\n" % state.name)
        TriggerPlot = plotregistry.get_plot('triggers')

        # --------------------------------------------------------------------
        # process time-series

        # find channels that need a TimeSeries
        if len(self.timeseries):
            vprint("    %d channels identified for TimeSeries\n"
                   % len(self.timeseries))
        get_timeseries_dict(self.timeseries, state, config=config, nds=nds,
                            multiprocess=multiprocess, return_=False)
        if len(self.timeseries):
            vprint("    All time-series data loaded\n")

        # find channels that need a StateVector
        if len(self.statevectors):
            vprint("    %d channels identified as StateVectors\n"
                   % len(self.statevectors))
        get_timeseries_dict(self.statevectors, state, config=config, nds=nds,
                            multiprocess=multiprocess, statevector=True,
                            return_=False)
        if len(self.statevectors):
            vprint("    All state-vector data loaded\n")

        # --------------------------------------------------------------------
        # process spectrograms

        # find FFT parameters
        try:
            fftparams = dict(config.nditems('fft'))
        except NoSectionError:
            fftparams = {}

        for channel in sorted(set(self.spectra + self.spectrograms),
                              key=lambda s: str(s)):
            get_spectrogram(channel, state, config=config, return_=False,
                            multiprocess=multiprocess, **fftparams)

        # --------------------------------------------------------------------
        # process spectra

        spectrumchannels = set()
        for plot in self.plots.spectra:
            spectrumchannels.update(plot.channels)
        for channel in sorted(spectrumchannels, key=lambda s: str(s)):
            get_spectrum(channel, state, config=config, return_=False,
                         **fftparams)

        # --------------------------------------------------------------------
        # process segments

        # find flags that need a DataQualityFlag
        if len(self.dataqualityflags):
            vprint("    %d data-quality flags identified for SegDB query\n"
                   % len(self.dataqualityflags))
            get_segments(self.dataqualityflags, state, config=config)

        # --------------------------------------------------------------------
        # process triggers

        tchannels = set()
        for plot in self.plots:
            if isinstance(plot, TriggerPlot):
                for channel in plot.channels:
                    tchannels.add((plot.etg, channel))
        for etg, channel in tchannels:
            get_triggers(channel, etg, state.active, config=config)

        # make plots
        vprint("    Plotting... \n")
        new_plots = [p for p in self.plots if
                     p.state is None or p.state.name == state.name and
                     not p.outputfile in globalv.WRITTEN_PLOTS]
        new_mp_plots = [p for p in new_plots if not isinstance(p, TriggerPlot)]
        processes = []
        for plot in sorted(new_plots,
                           key=lambda p: isinstance(p,
                                                    TriggerPlot) and 2 or 1):
            globalv.WRITTEN_PLOTS.append(plot.outputfile)
            if (multiprocess and len(new_mp_plots) > 1 and not
                    isinstance(plot, TriggerPlot)):
                p = Process(target=plot.process)
                processes.append(p)
                p.start()
            else:
                plot.process()
                vprint("        %s written\n" % plot.outputfile)
        if len(processes):
            vprint("        %d plot processes spawned, waiting"
                   % len(processes))
        for process in processes:
            process.join()
            vprint(".")
        if len(processes):
            vprint("\n")
        vprint("    Done.\n")

    # -------------------------------------------------------------------------
    # HTML operations

    def build_inner_html(self, state):
        """Write the '#main' HTML content for this tab.

        For now, this function just links all the plots in a 2-column
        format.
        """
        plots = [p for p in self.plots if
                 p.state is None or p.state.name == state.name]
        page = html.markup.page()

        # get layout
        if self.layout:
            layout = list(self.layout)
        else:
            layout = len(plots) == 1 and [1] or [2]
        while sum(layout) < len(plots):
            layout.append(layout[-1])
        l = i = 0
        for j, plot in enumerate(plots):
            # start new row
            if i == 0:
                page.div(class_='row')
            # make plot in its own column
            try:
                page.div(class_='col-md-%d' % (12 // layout[l]))
            except IndexError:
                warnings.warn("Something went wrong with the layout. "
                              "Tried to access element %d of ths following "
                              "layout (%d plots): %s" % (l, len(plots), layout))
                page.div(class_='col-md-%d' % 12 // layout[-1])
            page.a(href=plot.outputfile, class_='fancybox plot',
                   **{'data-fancybox-group': '1'})
            page.img(src=plot.outputfile)
            page.a.close()
            page.div.close()
            # detect end of row
            if (i + 1) == layout[l]:
                i = 0
                l += 1
                page.div.close()
            # detect last plot
            elif j == (len(plots) - 1):
                page.div.close()
                break
            # or move to next column
            else:
                i += 1
        # link data
        page.hr(class_='row-divider')
        page.div(class_='row')
        page.div(class_='col-md-12')
        if len(self.channels):
            page.h1('Channel information')
            page.add("The following channels were used to generate the above "
                     "data")
            headers = ['Channel', 'Type', 'Sample rate', 'Units']
            data = []
            for channel in self.channels:
                channel = get_channel(channel)
                # format CIS url and type
                if re.search('.[a-z]+\Z', channel.name):
                    name, ctype = channel.name.rsplit('.', 1)
                    c2 = get_channel(name)
                    cype = ctype in ['rms'] and ctype.upper() or ctype.title()
                else:
                    c2 = channel
                    ctype = 'Raw'
                if c2.url:
                    link = html.markup.oneliner.a(str(channel),
                                                  href=c2.url,
                                                  target='_blank')
                else:
                    link = str(channel)

                # format sameple rate
                if (channel.sample_rate is not None and
                        isclose(channel.sample_rate.value, 1/60.)):
                    rate = '1/60 %s' % channel.sample_rate.unit
                else:
                    rate = str(channel.sample_rate)
                # format unit
                if channel.unit:
                    unit = str(channel.unit)
                else:
                    unit = 'Unknown'
                data.append([link, ctype, rate, unit])
            page.add(str(html.data_table(headers, data, table='data')))
        if len(self.dataqualityflags):
            page.h1('Data-quality flag information')
            page.add("The following data-quality flags were used to generate "
                     "the above data. This list does not include state "
                     "information")
            # make summary table
            headers = ['IFO', 'Name', 'Version', 'Defined duration',
                       'Active duration']
            data = []
            pc = abs(state.active) / 100.
            for flag in sorted(self.dataqualityflags):
                flag = globalv.SEGMENTS[flag].copy()
                flag.valid &= state.active
                v = flag.version and str(flag.version) or ''
                try:
                    valid = '%.2f (%.2f%%)' % (abs(flag.valid),
                                               abs(flag.valid) / pc)
                except ZeroDivisionError:
                    valid = '0.00 (0.00%)'
                    active = '0.00 (0.00%)'
                else:
                    active = '%.2f (%.2f%%)' % (abs(flag.active),
                                                abs(flag.active) / pc)
                data.append([flag.ifo, flag.name, v, valid, active])
            page.add(str(html.data_table(headers, data, table='data')))
            # print segment lists
            page.div(class_='panel-group', id="accordion")
            for i,flag in enumerate(self.dataqualityflags):
                flag = globalv.SEGMENTS[flag].copy()
                flag.valid &= state.active
                n = ':'.join([flag.ifo, flag.name])
                if flag.version:
                    n += ':%d' % flag.version
                page.div(class_='panel panel-default')
                page.a(href='#flag%d' % i, **{'data-toggle': 'collapse',
                                              'data-parent': '#accordion'})
                page.div(class_='panel-heading')
                page.h4(n, class_='panel-title')
                page.div.close()
                page.a.close()
                page.div(id_='flag%d' % i, class_='panel-collapse collapse')
                page.div(class_='panel-body')
                # write segment summary
                page.p('This flag was defined and had a known state during '
                       'the following segments:')
                segwizard = StringIO()
                tosegwizard(segwizard, flag.valid)
                page.pre(segwizard.getvalue())
                segwizard.close()
                # write segment table
                page.p('This flag was active during the following segments:')
                segwizard = StringIO()
                flag.write(segwizard, format='segwizard')
                page.pre(segwizard.getvalue())
                segwizard.close()
                page.div.close()
                page.div.close()
                page.div.close()
            page.div.close()
        page.div.close()
        page.div.close()
        return page

register_tab(SimpleStateTab)


class AboutTab(Tab):
    type = 'about'

    @staticmethod
    def build_inner_html(config=None):
        return html.about_this_page(config=config)

register_tab(AboutTab)
