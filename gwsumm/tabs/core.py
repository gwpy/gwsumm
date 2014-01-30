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

import os
import re
import warnings
from StringIO import StringIO
from multiprocessing import Process

from numpy import (cumsum, isclose)

from lal import gpstime

from glue.segmentsUtils import tosegwizard

from gwpy.segments import Segment
from gwpy.io import nds as ndsio

from ..plot import (PlotList, registry as plotregistry)
from .. import globalv
from ..data import (get_channel, get_timeseries_dict, get_spectrogram,
                    get_spectrum)
from ..segments import get_segments
from ..triggers import get_triggers
from ..utils import (re_cchar, vprint)
from ..config import *
from .. import html

from gwsumm import version
__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

TriggerPlot = plotregistry.get_plot('triggers')
re_channel = re.compile('[A-Z]\d:[A-Z]+-[A-Z0-9_]+\Z')


class SummaryTab(object):
    """A single, self-enclosed summary webpage.

    This page will be linked to through the HTML navigation bar, hence
    the name 'tab', and can be grouped with other tabs through
    parent/child connections.
    """
    type = 'default'

    def __init__(self, name, parent=None, children=list(), states=None,
                 base='', layout=None, span=None, longname=None):
        """Initialise a new :class:`SummaryTab`
        """
        self.name = name
        self.longname = longname or self.name
        self.parent = parent
        self.children = children
        self.base = base
        self.plots = PlotList()
        self.states = states
        self.layout = layout
        if span is None:
            self.span = None
        else:
            self.span = Segment(*span)

    # -------------------------------------------
    # SummaryTab properties

    @property
    def name(self):
        """Short name for this tab, to be displayed in the HTML header
        and the navigation bar.

        :type: `str`
        """
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    @property
    def parent(self):
        """Short name of the parent page for this tab.

        A given tab can either be a parent for a set of child tabs, or can
        have a parent, it cannot be both. In this system, the `parent`
        attribute defines the heading under which this tab will be linked
        in the HTML navigation bar.

        :type: `str`
        """
        return self._parent

    @parent.setter
    def parent(self, p):
        self._parent = p

    @property
    def children(self):
        """List of child tabs for this tab.

        If this tab is given children, it cannot also have a parent, as it
        will define its own dropdown menu in the HTML navigation bar, linking
        to itself and its children.

        :type: `list` of `tabs <SummaryTab>`
        """
        return self._children

    @children.setter
    def children(self, clist):
        if self.parent and clist:
            raise ValueError("A SummaryTab cannot have both a parent, and a"
                             "set of children.")
        self._children = list(clist)

    @property
    def states(self):
        """The set of :class:`states <gwsumm.state.SummaryState>` over
        whos times this tab's data will be processed.

        The set of states will be linked in the given order with a switch
        on the far-right of the HTML navigation bar.
        """
        return self._states

    @states.setter
    def states(self, slist):
        if slist is None:
            self._states = None
        else:
            self._states = list(slist)

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
    def href(self):
        """HTML href attribute for this tab, relative to some base.

        By default all tabs are written into their own directory, so that the
        HTML href is just a directory reference inside the overall directory
        hierarchy.

        :type: `str`
        """
        try:
            return self._href
        except AttributeError:
            self._href = self.path + os.sep
            return self.href

    @href.setter
    def href(self, url):
        self._href = os.path.normpath(url) + os.sep

    @property
    def index(self):
        """URL to the base index for this tab.

        By default the index is the ``index.html`` file inside the href
        directory.

        :type: `str`
        """
        return os.path.normpath(os.path.join(self.href, 'index.html'))

    @property
    def path(self):
        """Path of this tab's directory relative to the --output-dir
        """
        if self.name.lower() == 'summary':
            p = ''
        else:
            p = re_cchar.sub('_', self.name).lower()
        tab_ = self
        while tab_.parent:
            p = os.path.join(re_cchar.sub('_', tab_.parent.name).lower(), p)
            tab_ = tab_.parent
        if self.base:
            return os.path.normpath(os.path.join(self.base, p))
        else:
            return os.path.normpath(p)

    @property
    def channels(self):
        """Set of all data channels used by this tab
        """
        out = set()
        for plot in self.plots:
            if (not isinstance(plot, TriggerPlot) and
                    hasattr(plot, 'channels')):
                out.update(plot.channels)
        return sorted(out, key=lambda ch: ch.name)

    @property
    def trigchannels(self):
        """Set of all trigger channels used by this tab
        """
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

    @property
    def hrefs(self):
        # write page for each state
        statelinks = []
        for i, state in enumerate(self.states):
            statelinks.append(os.path.join(
                self.href, '%s.html' % re_cchar.sub('_', state.name.lower())))
        return statelinks

    # -------------------------------------------
    # SummaryTab methods

    def add_child(self, tab):
        """Add a child to this `SummaryTab`

        Parameters
        ----------
        tab : `SummaryTab`
            child tab to record
        """
        self.children.append(tab)

    def get_child(self, name):
        """Find a child tab of this `SummaryTab` by name

        Parameters
        ----------
        name : `str`
            string identifier of child tab to use in search

        Returns
        -------
        child : `SummaryTab`
            the child tab found by name

        Raises
        ------
        RuntimeError
            if no child tab can be found matching the given ``name``
        """
        names = [c.name for c in self.children]
        try:
            idx = names.index(name)
        except ValueError:
            raise RuntimeError("This tab has no child named '%s'." % name)
        else:
            return self.children[idx]

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
            statenames = cp.get(section, 'states').split(',')
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

            # define one copy of this plot for each state
            for state in job.states:
                # if pdef refers to another config section, look there
                # for a detailed definition
                if cp.has_section(pdef):
                    type_ = cp.get(pdef, 'type')
                    PlotClass = plotregistry.get_plot(type_)
                    plot = PlotClass.from_ini(cp, pdef, state=state,
                                              outdir=plotdir, **mods)
                    for source in sources:
                        plot.add_data_source(source)
                # otherwise pdef should be a registered plot type
                else:
                    PlotClass = plotregistry.get_plot(pdef)
                    plot = PlotClass(sources, start, end, state=state,
                                     outdir=plotdir, **mods)
                job.plots.append(plot)

        return job

    # -------------------------------------------
    # SummaryTab processing

    def finalize_states(self, config=GWSummConfigParser()):
        """Fetch the segments for each state for this `SummaryTab`
        """
        for state in self.states:
            state.fetch(config=config)
        self.states.sort(key=lambda s: abs(s.active), reverse=True)

    def process(self, nds='guess', multiprocess=True,
                config=GWSummConfigParser()):
        """Process data for this tab
        """
        vprint("\n-------------------------------------------------\n")
        if self.parent:
            vprint("Processing %s/%s\n" % (self.parent.name, self.name))
        else:
            vprint("Processing %s\n" % self.name)
        self.finalize_states(config=config)
        vprint("States finalised\n")
        for state in sorted(self.states, key=lambda s: abs(s.active),
                            reverse=True):
            self.process_state(state, nds=nds, multiprocess=multiprocess,
                               config=config)

    def process_state(self, state, nds='guess', multiprocess=True,
                      config=GWSummConfigParser()):
        """Process data for this tab in a given state
        """
        vprint("Processing '%s' state\n" % state.name)

        # --------------------------------------------------------------------
        # process time-series

        # find channels that need a TimeSeries
        if len(self.channels):
            vprint("    %d channels identified for TimeSeries\n"
                   % len(self.channels))
        get_timeseries_dict(self.channels, state, config=config, nds=nds,
                            multiprocess=multiprocess, return_=False)
        if len(self.channels):
            vprint("    All time-series data loaded\n")

        # --------------------------------------------------------------------
        # process spectrograms

        # find FFT parameters
        try:
            fftparams = dict(config.nditems('fft'))
        except NoSectionError:
            fftparams = {}

        spectrochannels = set()
        for plot in self.plots.spectrograms + self.plots.spectra:
            spectrochannels.update(plot.channels)
        for channel in sorted(spectrochannels, key=lambda s: str(s)):
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
                    tchannels.update(set([(plot.etg, channel)]))
        for etg, channel in tchannels:
            get_triggers(channel, etg, state.active, config=config)

        # make plots
        vprint("    Plotting... \n")
        new_plots = [p for p in self.plots if p.state.name == state.name and
                     not p.outputfile in globalv.WRITTEN_PLOTS and
                     not isinstance(p, TriggerPlot)]
        new_mp_plots = [p for p in new_plots if not isinstance(p, TriggerPlot)]
        processes = []
        for plot in sorted(new_plots, key=lambda p:
                                        isinstance(p, TriggerPlot) and 2 or 1):
            globalv.WRITTEN_PLOTS.append(plot.outputfile)
            if (multiprocess and len(new_mp_plots) and not
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
            p.join()
            vprint(".")
        if len(processes):
            vprint("\n")
        vprint("    Done.\n")

    # -------------------------------------------------------------------------
    # HTML operations

    def write_calendar_link(self):
        if globalv.MODE < 4:
            date = gpstime.gps_to_utc(self.states[0].extent[0])
            cal = str(html.calendar(date, mode=globalv.MODE % 3))
        else:
            start, end = self.states[0].extent
            cal = html.markup.oneliner.p('%d-%d' % (start, end),
                                         class_='navbar-text')
        return cal

    def build_navigation_links(self, tabs=[]):
        navlinks = []
        for tab in tabs:
            if len(tab.children):
                navlinks.append([tab.name, []])
                links = [(child.name, child.href) for child in tab.children]
                if self in tab.children:
                    active = tab.children.index(self)
                else:
                    active = None
                if tab.children[0].name == 'Summary':
                    links.insert(1, None)
                    if active > 0:
                        active += 1
                navlinks[-1][1].extend(links)
                navlinks[-1].append(active)
            else:
                navlinks.append((tab.name, tab.href))
        return navlinks

    def build_html(self, title=None, subtitle=None, tabs=list(), css=list(),
                   js=list(), about=None, writedata=True, writehtml=True):
        """Construct the HTML page for this tab.
        """
        # find relative base path
        n = len(self.index.split(os.path.sep)) - 1
        base = os.sep.join([os.pardir] * n)
        # work title as Parent Name/Tab Name
        if not title and not subtitle:
            if self.parent:
                title = self.longname
                tab_ = self
                while tab_.parent:
                    title = '%s/%s' % (tab_.parent.longname, title)
                    tab_ = tab_.parent
                title, subtitle = title.rsplit('/', 1)
            elif self.name == 'Summary':
                title = 'IFO summary'
                subtitle = '%d-%d' % (self.span[0], self.span[1])
            else:
                title = self.longname
                subtitle = None
        # get calendar and write navigation bar
        brand = self.write_calendar_link()
        navlinks = self.build_navigation_links(tabs)
        # initialise page
        page = html.markup.page()
        page.init(doctype=html.DOCTYPE, css=css, script=js, title=title,
                  base=base, metainfo=html.META)
        page.div(id_='wrap')
        page.add(str(html.banner(title, subtitle=subtitle)))
        if len(self.states) > 1:
            statebtn = html.state_switcher(zip(self.states, self.hrefs), 0)
        else:
            statebtn = False
        page.add(str(html.navbar(navlinks, brand=brand, states=statebtn)))
            # write state pages and include first one
        # write page for each state
        if writedata:
            for i, (state, shtml) in enumerate(zip(self.states, self.hrefs)):
                with open(shtml, 'w') as fobj:
                    fobj.write(str(self.write_html(state)))
        page.add(str(html.wrap_content('')))
        if len(self.states) > 1:
            page.add(str(html.load_state(self.hrefs[0])))
        else:
            page.add(str(html.load(self.hrefs[0])))
        # add footer
        page.div.close()
        page.add(str(html.footer(about=about)))
        # write
        if writehtml:
            with open(self.index, 'w') as fobj:
                fobj.write(str(page))
        return

    def write_html(self, state):
        """Write the '#main' HTML content for this tab.

        For now, this function just links all the plots in a 2-column
        format.
        """
        plots = [p for p in self.plots if p.state.name == state.name]
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
                flag = globalv.SEGMENTS[flag]
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
                flag = globalv.SEGMENTS[flag]
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


class AboutTab(SummaryTab):

    def build_html(self, tabs=list(), css=list(), js=list(), config=None):
        """Construct the HTML page for this `AboutTab`.
        """
        # find relative base path
        n = len(self.index.split(os.path.sep)) - 1
        base = os.sep.join([os.pardir] * n)
        # get calendar and write navigation bar
        brand = self.write_calendar_link()
        navlinks = self.build_navigation_links(tabs)
        # write page for each state
        # initialise page
        page = html.markup.page()
        page.init(doctype=html.DOCTYPE, css=css, script=js, title=self.name,
                  base=base)
        page.div(id_='wrap')
        page.add(str(html.banner(self.name)))
        page.add(str(html.navbar(navlinks, brand=brand)))
        # write state pages and include first one
        page.add(str(html.wrap_content(self.write_html(config=config))))
        # add footer
        page.div.close()
        page.add(str(html.footer(about=os.path.join(base, self.index))))
        page.body.close()
        page.html.close()
        # write
        with open(self.index, 'w') as fobj:
            fobj.write(str(page))
        return

    def write_html(self, config=None):
        return html.about_this_page(config=config)


def split_channels(channelstring):
    """Split a comma-separated list of channels that may, or may not
    contain NDS2 channel types as well
    """
    out = []
    while True:
        channelstring = channelstring.strip(' \n')
        if ',' not in channelstring:
            break
        for nds2type in ndsio.NDS2_CHANNEL_TYPE.keys() + ['']:
            if nds2type and ',%s' % nds2type in channelstring:
                try:
                    channel, ctype, channelstring = channelstring.split(',', 2)
                except ValueError:
                    channel, ctype = channelstring.split(',')
                    channelstring = ''
                out.append('%s,%s' % (channel, ctype))
                break
            elif nds2type == '' and ',' in channelstring:
                channel, channelstring = channelstring.split(',', 1)
                out.append(channel)
                break
    if channelstring:
        out.append(channelstring)
    return out
