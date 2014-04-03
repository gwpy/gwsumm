# coding=utf-8
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

"""DOCSTRING
"""

import os
import multiprocessing

from lal import gpstime

from gwpy.segments import Segment

from .. import (globalv, html, version)
from ..utils import (re_cchar, vprint, count_free_cores)
from ..config import *
from .registry import register_tab

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


class Tab(object):
    """A configurable output HTML page.

    This class defines a stub from which actual output format tabs
    should be sub-classed.

    All sub-classes must define the following methods

    .. autosummary::
       :nosignatures:

       Tab.from_ini
       Tab.process
       Tab.build_inner_html
    """
    type = 'basic'

    def __init__(self, name, longname=None, parent=None, children=list(),
                 base='', span=None):
        """Initialise a new `Tab`.
        """
        self.name = name
        self.longname = longname or self.name
        self.parent = parent
        self.children = children
        self.base = base
        if span is None:
            self.span = None
        else:
            self.span = Segment(*span)

    # -------------------------------------------
    # Tab properties

    @property
    def name(self):
        """Short name for this `Tab`.

        This will be displayed in the navigation bar.

        :type: `str`
        """
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    @property
    def parent(self):
        """Short name of the parent page for this `Tab`.

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
        """List of child tabs for this `Tab`.

        If this tab is given children, it cannot also have a parent, as it
        will define its own dropdown menu in the HTML navigation bar, linking
        to itself and its children.

        :type: `list` of `tabs <Tab>`
        """
        return self._children

    @children.setter
    def children(self, clist):
        if self.parent and clist:
            raise ValueError("A Tab cannot have both a parent, and a"
                             "set of children.")
        self._children = list(clist)

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

    # -------------------------------------------
    # Tab instance methods

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
    def from_ini(cls, *args, **kwargs):
        """Define a new `Tab` from the given section of a `ConfigParser`.

        This method must be defined by all sub-classes.
        """
        raise NotImplementedError("This method should be defined by all "
                                  "sub-classes.")

    # -------------------------------------------
    # SummaryTab processing

    def process(self, *args, **kwargs):
        """Process data for this `Tab`.

        This method must be defined by all sub-classes.
        """
        raise NotImplementedError("This method should be defined by all "
                                  "sub-classes.")

    # -------------------------------------------------------------------------
    # HTML operations

    def initialise_html_page(self, title=None, subtitle=None, css=list(),
                             js=list()):
        """Initialise a new HTML `~markup.page` with a headline.

        Parameters
        ----------
        title : `str`, optional, default: {parent.longname}
            level 1 heading for this `Tab`.
        subtitle : `str`, optional, default: {self.longname}
            level 2 heading for this `Tab`.
        css : `list`, optional
            list of resolvable URLs for CSS files
        js : `list`, optional
            list of resolvable URLs for javascript files

        Returns
        -------
        page : :class:`markup.page`
            initialised HTML markup page
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
        # initialise page
        page = html.markup.page()
        page.init(doctype=html.DOCTYPE, css=css, script=js, title=title,
                  base=base, metainfo=html.META)
        page.div(id_='wrap')
        page.add(str(html.banner(title, subtitle=subtitle)))
        return page

    def build_html_navbar(self, ifo=None, ifomap=None, tabs=list()):
        """Build the navigation bar for this `Tab`.

        Parameters
        ----------
        ifo : `str`, optional
            prefix for this IFO.
        ifomap : `dict`, optional
            `dict` of (ifo, {base url}) pairs to map to summary pages for
            other IFOs.

            .. note::
               The ``ifo`` and ``ifomap`` parameters should be given
               together.

        tabs : `list`, optional
            list of parent tabs (each with a list of children) to include
            in the navigation bar.

        Returns
        -------
        page : `markup.page`
            a markup page containing the navigation bar.
        """
        # build interferometer cross-links
        if ifo is not None and ifomap is not None:
            ifolinks = str(html.base_map_dropdown(ifo, **ifomap))
        else:
            ifolinks = ''
        # build calendar
        if globalv.MODE < 4:
            date = gpstime.gps_to_utc(self.span[0])
            cal = str(html.calendar(date, mode=globalv.MODE % 3))
        else:
            start, end = self.span
            cal = html.markup.oneliner.p('%d-%d' % (start, end),
                                         class_='navbar-text')
        # combine as 'brand'
        brand = ifolinks + cal
        # build navbar links
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
        # combine and return
        return html.navbar(navlinks, brand=brand)

    def write_html(self, outfile=None, title=None, subtitle=None, tabs=list(),
                   ifo=None, ifomap=None, css=list(), js=list(), about=None,
                   writedata=True, writehtml=True, **inargs):
        """Write the HTML page for this `Tab`.

        Parameters
        ----------
        outfile : `str`, optional
            path of file to write HTML into
        title : `str`, optional, default: {parent.longname}
            level 1 heading for this `Tab`.
        subtitle : `str`, optional, default: {self.longname}
            level 2 heading for this `Tab`.
        tabs: `list`, optional
            list of top-level tabs (with children) to populate navbar
        ifo : `str`, optional
            prefix for this IFO.
        ifomap : `dict`, optional
            `dict` of (ifo, {base url}) pairs to map to summary pages for
            other IFOs.

            .. note::
               The ``ifo`` and ``ifomap`` parameters should be given
               together.

        css : `list`, optional
            list of resolvable URLs for CSS files
        js : `list`, optional
            list of resolvable URLs for javascript files
        about : `str`, optional
            href for the 'About' page
        writedata : `bool`, optional, default: `True`
            if `True`, actually write the state frames
        writehtml : `bool`, optional, default: `True`
            if `True`, actually write the outer-frame HTML
        **inargs
            other keyword arguments to pass to the
            :meth:`~Tab.build_inner_html` method
        """
        # initialise HTML page
        page = self.initialise_html_page(title=title, subtitle=subtitle,
                                         css=css, js=js)
        # add navigation
        page.add(str(self.build_html_navbar(ifo=ifo, ifomap=ifomap,
                                            tabs=tabs)))
        # write state pages and include first one
        page.add(str(html.wrap_content(self.build_inner_html(**inargs))))

        # add footer
        page.div.close()
        page.add(str(html.footer(about=about)))
        page.body.close()
        page.html.close()
        # write
        if outfile is None:
            outfile = self.index
        with open(outfile, 'w') as fobj:
            fobj.write(str(page))
        return

    def build_inner_html(self, *args, **kwargs):
        """Build the '#main' HTML content for this `Tab`.

        This method should be defined by all sub-classes.
        """
        raise NotImplementedError("This method should be defined by all "
                                  "sub-classes.")

register_tab(Tab)


class StateTab(Tab):
    """`Tab` with distinct `States <~gwsumm.states.core.SummaryState>`.

    This is another stub to make generate other `Tab` classes with
    multiple states easier.

    This class defines the :meth:`~StateTab.process` method, leaving
    all sub-classes to define the following methods

    .. autosummary::
       :nosignatures:

       Tab.from_ini
       Tab.process_state
       Tab.build_inner_html

    The latter two should operate on a single `SummaryState.
    """
    type = 'state'

    def __init__(self, name, longname=None, parent=None, children=list(),
                 base='', span=None, states=None):
        """Initialise a new `Tab`.
        """
        super(StateTab, self).__init__(name, longname=longname, parent=parent,
                                       children=children, base=base, span=span)
        self.states = states

    # -------------------------------------------
    # StateTab properties

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
    def hrefs(self):
        # write page for each state
        statelinks = []
        for i, state in enumerate(self.states):
            statelinks.append(os.path.join(
                self.href, '%s.html' % re_cchar.sub('_', state.name.lower())))
        return statelinks

    # -------------------------------------------
    # StateTab methods

    def finalize_states(self, config=ConfigParser()):
        """Fetch the segments for each state for this `SummaryTab`
        """
        for state in self.states:
            state.fetch(config=config)

    def process(self, config=ConfigParser(), multiprocess=True, **stateargs):
        """Process data for this `StateTab`.

        Parameters
        ----------
        config : `ConfigParser.ConfigParser`, optional
            job configuration to pass to :math:`~StateTab.finalize_states`
        **stateargs
            all other keyword arguments are passed directly onto the
            :meth:`~StateTab.process_state` method.
        """
        vprint("\n-------------------------------------------------\n")
        if self.parent:
            vprint("Processing %s/%s\n" % (self.parent.name, self.name))
        else:
            vprint("Processing %s\n" % self.name)
        # load state segments
        self.finalize_states(config=config)
        vprint("States finalised\n")
        # setup plotting queue
        if isinstance(multiprocess, int):
            queue = multiprocessing.JoinableQueue(
                count_free_cores(multiprocess))
        elif multiprocess:
            queue = multiprocessing.JoinableQueue(count_free_cores())
        else:
            queue = None
        # process each state
        for state in sorted(self.states, key=lambda s: abs(s.active),
                            reverse=True):
            self.process_state(state, config=config, multiprocess=multiprocess,
                               plotqueue=queue, **stateargs)

        # consolidate child processes
        if queue is not None:
            vprint("Waiting for plotting processes to complete.. ")
            queue.close()
            queue.join()
            vprint('done.\n')
        vprint("%s/%s complete!\n" % (self.parent.name, self.name))

    def process_state(self, state, **kwargs):
        """Process data for this tab in a given state.

        This method should be provided by all sub-classes, in order to
        process all relevant data and figures required for this Tab's
        HTML output.

        Parameters
        ----------
        state : `~gwsumm.state.SummaryState`
            `SummaryState` over which to generate inner HTML
        **kwargs
            any and all other keyword arguments
        """
        raise NotImplementedError("This method should be provided by all "
                                  "sub-classes.")

    def build_html_navbar(self, ifo=None, ifomap=None, tabs=list()):
        """Build the navigation bar for this `Tab`.

        Parameters
        ----------
        ifo : `str`, optional
            prefix for this IFO.
        ifomap : `dict`, optional
            `dict` of (ifo, {base url}) pairs to map to summary pages for
            other IFOs.

            .. note::
               The ``ifo`` and ``ifomap`` parameters should be given
               together.

        tabs : `list`, optional
            list of parent tabs (each with a list of children) to include
            in the navigation bar.

        Returns
        -------
        page : `markup.page`
            a markup page containing the navigation bar.
        """
        # build interferometer cross-links
        if ifo is not None and ifomap is not None:
            ifolinks = str(html.base_map_dropdown(ifo, **ifomap))
        else:
            ifolinks = ''
        # build calendar
        if globalv.MODE < 4:
            date = gpstime.gps_to_utc(self.span[0])
            cal = str(html.calendar(date, mode=globalv.MODE % 3))
        else:
            start, end = self.span
            cal = html.markup.oneliner.p('%d-%d' % (start, end),
                                         class_='navbar-text')
        # combine as 'brand'
        brand = ifolinks + cal
        # build navbar links
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
        # build state switch
        if len(self.states) > 1:
            statebtn = html.state_switcher(zip(self.states, self.hrefs), 0)
        else:
            statebtn = False
        # combine and return
        return html.navbar(navlinks, brand=brand, states=statebtn)

    def write_html(self, outfile=None, title=None, subtitle=None, tabs=list(),
                   ifo=None, ifomap=None, css=list(), js=list(), about=None,
                   writedata=True, writehtml=True, **inargs):
        """Construct the HTML page for this `StateTab`.

        Parameters
        ----------
        outfile : `str`, optional
            path of file to write HTML into
        title : `str`, optional, default: {parent.longname}
            level 1 heading for this `Tab`.
        subtitle : `str`, optional, default: {self.longname}
            level 2 heading for this `Tab`.
        tabs: `list`, optional
            list of top-level tabs (with children) to populate navbar
        ifo : `str`, optional
            prefix for this IFO.
        ifomap : `dict`, optional
            `dict` of (ifo, {base url}) pairs to map to summary pages for
            other IFOs.

            .. note::
               The ``ifo`` and ``ifomap`` parameters should be given
               together.

        css : `list`, optional
            list of resolvable URLs for CSS files
        js : `list`, optional
            list of resolvable URLs for javascript files
        about : `str`, optional
            href for the 'About' page
        writedata : `bool`, optional, default: `True`
            if `True`, actually write the state frames
        writehtml : `bool`, optional, default: `True`
            if `True`, actually write the outer-frame HTML
        """
        # initialise HTML page
        page = self.initialise_html_page(title=title, subtitle=subtitle,
                                         css=css, js=js)
        # add navigation
        page.add(str(self.build_html_navbar(ifo=ifo, ifomap=ifomap,
                                            tabs=tabs)))
        # write page for each state
        if writedata:
            for i, (state, shtml) in enumerate(zip(self.states, self.hrefs)):
                inner = self.build_inner_html(state, **inargs)
                with open(shtml, 'w') as fobj:
                    fobj.write(str(inner))
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
            if outfile is None:
                outfile = self.index
            with open(outfile, 'w') as fobj:
                fobj.write(str(page))
        return

    def build_inner_html(self, state, **kwargs):
        """Build the '#main' HTML content for this `Tab`.

        This method should be defined by all sub-classes to execute on a
        single state, with the following argument format

        Parameters
        ----------
        state : `~gwsumm.state.SummaryState`
            `SummaryState` over which to generate inner HTML
        **kwargs
            any and all other keyword arguments
        """
        raise NotImplementedError("This method should be defined by all "
                                  "sub-classes.")

register_tab(StateTab)
