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

"""This module defines a number of `Tab` subclasses.
"""

from gwpy.segments import (Segment, DataQualityFlag)
from gwpy.time import (from_gps, to_gps)

from .core import (Tab, SummaryArchiveMixin)
from .registry import (get_tab, register_tab)
from ..plot import SummaryPlot
from ..utils import *
from ..config import *
from .. import html

from gwsumm import version
__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


class ExternalTab(Tab):
    """A simple tab to link HTML from an external source

    Parameters
    ----------
    name : `str`
        name of this tab (required)
    url : `str`
        URL of the external content to be linked into this tab.
    index : `str`
        HTML file in which to write. By default each tab is written to
        an index.html file in its own directory. Use :attr:`~Tab.index`
        to find out the default index, if not given.
    shortname : `str`
        shorter name for this tab to use in the navigation bar. By
        default the regular name is used
    parent : :class:`~gwsumm.tabs.Tab`
        parent of this tab. This is used to position this tab in the
        navigation bar.
    children : `list`
        list of child :class:`Tabs <~gwsumm.tabs.Tab>` of this one. This
        is used to position this tab in the navigation bar.
    group : `str`
        name of containing group for this tab in the navigation bar
        dropdown menu. This is only relevant if this tab has a parent.
    base : `str`
        path for HTML base attribute for this tab.
    """
    type = 'external'

    def __init__(self, name, url, **kwargs):
        """Initialise a new `ExternalTab`.
        """
        super(ExternalTab, self).__init__(name, **kwargs)
        self.url = url

    @property
    def url(self):
        return self._url

    @url.setter
    def url(self, link):
        self._url = link

    def build_html_content(self, content, class_='container', id_='main'):
        page = html.markup.page()
        page.div(content, class_=class_, id_=id_)
        page.add(str(html.load(self.url, id_=id_)))
        return page

    def write_html(self, **kwargs):
        """Write the HTML page for this tab.

        See Also
        --------
        gwsumm.tabs.Tab.write_html : for details of all valid keyword
        arguments
        """
        link = html.markup.given_oneliner.a('click here to view the original',
                                            class_='reference',
                                            href=self.url.split()[0])
        kwargs.setdefault('footer', 'This page contains data from an external '
                                    'source, %s.' % link)
        return super(ExternalTab, self).write_html('', **kwargs)

register_tab(ExternalTab)


class ArchivedExternalTab(SummaryArchiveMixin, ExternalTab):
    """An archivable externally-linked tab.
    """
    type = 'archived-external'

    def __init__(self, name, url, start, end, mode=None, **kwargs):
        super(ArchivedExternalTab, self).__init__(name, url, **kwargs)
        self.span = (start, end)
        self.mode = mode

register_tab(ArchivedExternalTab)


class PlotTab(Tab):
    """A simple tab to layout some figures in the #main div.

    Parameters
    ----------
    name : `str`
        name of this tab (required)
    plots : `list`
        list of plots to display on this tab. More plots can be added
        at any time via :meth:`PlotTab.add_plot`
    layout : `int`, `list`
        the number of plots to display in each row, or a list of numbers
        to define each row individually. If the number of plots defined
        by the layout is less than the total number of plots, the layout
        for the final row will be repeated as necessary.

        For example ``layout=[1, 2, 3]`` will display a single plot on
        the top row, two plots on the second, and 3 plots on each row
        thereafter.
    index : `str`
        HTML file in which to write. By default each tab is written to
        an index.html file in its own directory. Use :attr:`~Tab.index`
        to find out the default index, if not given.
    shortname : `str`
        shorter name for this tab to use in the navigation bar. By
        default the regular name is used
    parent : :class:`~gwsumm.tabs.Tab`
        parent of this tab. This is used to position this tab in the
        navigation bar.
    children : `list`
        list of child :class:`Tabs <~gwsumm.tabs.Tab>` of this one. This
        is used to position this tab in the navigation bar.
    group : `str`
        name of containing group for this tab in the navigation bar
        dropdown menu. This is only relevant if this tab has a parent.
    base : `str`
        path for HTML base attribute for this tab.
    """
    type = 'plots'

    def __init__(self, name, plots=[], layout=None, **kwargs):
        """Initialise a new :class:`PlotTab`.
        """
        super(PlotTab, self).__init__(name, **kwargs)
        self.plots = []
        for p in plots:
            self.add_plot(p)
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
        self._layout = l

    @classmethod
    def from_ini(cls, cp, section, base=''):
        """Define a new tab from a :class:`~gwsumm.config.GWConfigParser`

        Parameters
        ----------
        cp : :class:`~gwsumm.config.GWConfigParser`
            customised configuration parser containing given section
        section : `str`
            name of section to parse

        Returns
        -------
        tab : `PlotTab`
            a new tab defined from the configuration
        """
        # get tab name
        try:
            # name given explicitly
            name = re_quote.sub('', cp.get(section, 'name'))
        except NoOptionError:
            # otherwise strip 'tab-' from section name
            name = section[4:]
        try:
            shortname = re_quote.sub('', cp.get(section, 'shortname'))
        except NoOptionError:
            shortname = name
        # get parent:
        #     if parent is not given, this assumes a top-level tab
        try:
            parent = re_quote.sub('', cp.get(section, 'parent'))
        except NoOptionError:
            parent = None
        else:
            if parent == 'None':
                parent = None
        # get group
        try:
            group = cp.get(section, 'group')
        except NoOptionError:
            group = None
        # get HTML file
        try:
            index = cp.get(section, 'index')
        except NoOptionError:
            index = None
        # get GPS start time
        start = cp.getfloat(section, 'gps-start-time')
        end = cp.getfloat(section, 'gps-end-time')
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
                if isinstance(l, (tuple, list)):
                    l = l[0]
                if not l in [1, 2, 3, 4, 6, 12]:
                    raise ValueError("Cannot print more than %d plots in a "
                                     "single row. The chosen layout value for "
                                     "each row must be a divisor of 12 to fit "
                                     "the Bootstrap scaffolding. For details "
                                     "see http://getbootstrap.com/2.3.2/"
                                     "scaffolding.html")
        else:
            layout = None

        # build and return tab
        new = cls(name, start, end, index=index, shortname=shortname,
                  parent=parent, group=group, base=base, layout=layout)
        return new

    def add_plot(self, plot):
        """Add a plot to this tab.

        Parameters
        ----------
        plot : `str`, :class:`~gwsumm.plot.SummaryPlot`
            either the URL of a plot to embed, or a formatted `SummaryPlot`
            object.
        """
        if isinstance(plot, str):
            plot = SummaryPlot(href=plot)
            plot.new = False
        if not isinstance(plot, SummaryPlot):
            raise TypeError("Cannot append plot of type %r" % type(plot))
        self.plots.append(plot)

    def scaffold_plots(self, state=None):
        """Build a grid of plots using bootstrap's scaffolding.

        Returns
        -------
        page : :class:`~gwsumm.html.markup.page`
            formatted markup with grid of plots
        """
        page = html.markup.page()

        if state:
            plots = [p for p in self.plots if p.state in [state, None]]
        else:
            plots = self.plots

        # get layout
        if self.layout:
            layout = list(self.layout)
        else:
            layout = len(plots) == 1 and [1] or [2]
        for i, l in enumerate(layout):
            if isinstance(l, (list, tuple)):
                layout[i] = l
            elif isinstance(l, int):
                layout[i] = (l, None)
            else:
                raise ValueError("Cannot parse layout element '%s'." % l)
        while sum(zip(*layout)[0]) < len(plots):
            layout.append(layout[-1])
        l = i = 0
        for j, plot in enumerate(plots):
            # start new row
            if i == 0:
                page.div(class_='row')
            # determine relative size
            if layout[l][1]:
                colwidth = 12 // int(layout[l][1])
                remainder = 12 - colwidth * layout[l][0]
                if remainder % 2:
                    raise ValueError("Cannot center column of width %d in a "
                                     "12-column format" % colwidth)
                else:
                    offset = remainder / 2
                page.div(class_='col-md-%d col-md-offset-%d'
                                % (colwidth, offset))
            else:
                colwidth = 12 // int(layout[l][0])
                page.div(class_='col-md-%d' % colwidth)
            page.a(href=plot.href, class_='fancybox plot',
                   **{'data-fancybox-group': '1'})
            page.img(src=plot.href)
            page.a.close()
            page.div.close()
            # detect end of row
            if (i + 1) == layout[l][0]:
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
        return page

    def build_html_content(self, sandwich, class_='container', id_='main'):
        page = html.markup.page()
        page.div(class_=class_, id_=id_)
        if isinstance(sandwich, (list, tuple)):
            pre, post = sandwich
        else:
            pre = sandwich
            post = None
        if pre:
            page.add(str(pre))
        page.add(str(self.scaffold_plots()))
        if post:
            page.add(str(post))
        page.div.close()
        return page

    def write_html(self, pre=None, post=None, **kwargs):
        """Write the HTML page for this tab.

        Parameters
        ----------
        pre : `str`, :class:`~gwsumm.html.markup.page`
            content to place above the plot grid
        post : `str`, :class:`~gwsumm.html.markup.page`
            content to place below the plot grid
        **kwargs
            other keyword arguments to be passed to :meth:`~Tab.write_html`

        See Also
        --------
        gwsumm.tabs.Tab.write_html : for details of all valid unnamed
        keyword arguments
        """
        if post:
            return super(PlotTab, self).write_html((pre, post), **kwargs)
        else:
            return super(PlotTab, self).write_html(pre, **kwargs)

register_tab(PlotTab)


class ArchivedPlotTab(SummaryArchiveMixin, PlotTab):

    """An archivable externally-linked tab.
    """
    type = 'archived-plots'

    def __init__(self, name, start, end, mode=None, plots=[], **kwargs):
        super(ArchivedPlotTab, self).__init__(name, plots=plots, **kwargs)
        self.span = (start, end)
        self.mode = mode

register_tab(ArchivedPlotTab)


class StateTab(PlotTab):
    """Tab with multiple content pages defined via 'states'

    Each state is printed to its own HTML file which is loaded via
    javascript upon request into the #main div of the index for the tab.

    Parameters
    ----------
    name : `str`
        name of this tab (required)
    states : `list`
        a list of states for this tab. Each state can take any form,
        but must be castable to a `str` in order to be printed.
    plots : `list`
        list of plots to display on this tab. More plots can be added
        at any time via :meth:`PlotTab.add_plot`
    layout : `int`, `list`
        the number of plots to display in each row, or a list of numbers
        to define each row individually. If the number of plots defined
        by the layout is less than the total number of plots, the layout
        for the final row will be repeated as necessary.

        For example ``layout=[1, 2, 3]`` will display a single plot on
        the top row, two plots on the second, and 3 plots on each row
        thereafter.
    index : `str`
        HTML file in which to write. By default each tab is written to
        an index.html file in its own directory. Use :attr:`~Tab.index`
        to find out the default index, if not given.
    shortname : `str`
        shorter name for this tab to use in the navigation bar. By
        default the regular name is used
    parent : :class:`~gwsumm.tabs.Tab`
        parent of this tab. This is used to position this tab in the
        navigation bar.
    children : `list`
        list of child :class:`Tabs <~gwsumm.tabs.Tab>` of this one. This
        is used to position this tab in the navigation bar.
    group : `str`
        name of containing group for this tab in the navigation bar
        dropdown menu. This is only relevant if this tab has a parent.
    base : `str`
        path for HTML base attribute for this tab.
    """
    type = 'state'

    def __init__(self, name, states=[], **kwargs):
        """Initialise a new `Tab`.
        """
        super(StateTab, self).__init__(name, **kwargs)
        # process states
        self._states = []
        if not isinstance(states, (tuple, list)):
            states = [states]
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
    def states(self, statelist):
        self._states = []
        for state in statelist:
            self.add_state(state)

    def add_state(self, state):
        """Add a `SummaryState` to this `StateTab`.

        Parameters
        ----------
        state : `str`, :class:`~gwsumm.state.SummaryState`
            either the name of a state, or a `SummaryState`.
        register : `bool`, default: `False`
            automatically register all new states
        """
        self._states.append(state)

    @property
    def frames(self):
        # write page for each state
        statelinks = []
        for i, state in enumerate(self.states):
            statelinks.append(os.path.join(
                self.path, '%s.html' % re_cchar.sub('_', str(state).lower())))
        return statelinks

    # -------------------------------------------
    # StateTab methods

    @classmethod
    def from_ini(cls, cp, section, base=''):
        # parse core Tab information
        new = super(StateTab, cls).from_ini(cp, section, base=base)
        # parse states and retrieve their definitions
        if cp.has_option(section, 'states'):
            # states listed individually
            statenames = [re_quote.sub('', s).strip() for s in
                          cp.get(section, 'states').split(',')]
        else:
            # otherwise use 'all' state - full span with no gaps
            statenames = ['All']
        new.states = statenames
        return new

    # ------------------------------------------------------------------------
    # HTML methods

    def build_html_navbar(self, brand=None, ifo=None, ifomap=dict(),
                          tabs=list()):
        """Build the navigation bar for this `Tab`.

        The navigation bar will consist of a switch for this page linked
        to other interferometer servers, followed by the navbar brand,
        then the full dropdown-based navigation menus configured for the
        given ``tabs`` and their descendents.

        Parameters
        ----------
        brand : `str`, :class:`~gwsumm.html.markup.page`
            content for navbar-brand
        ifo : `str`, optional
            prefix for this IFO.
        ifomap : `dict`, optional
            `dict` of (ifo, {base url}) pairs to map to summary pages for
            other IFOs.

        tabs : `list`, optional
            list of parent tabs (each with a list of children) to include
            in the navigation bar.

        Returns
        -------
        page : `~gwsumm.html.markup.page`
            a markup page containing the navigation bar.
        """
        # build interferometer cross-links
        if ifo is not None:
            brand_ = html.base_map_dropdown(ifo, id_='ifos', **ifomap)
        else:
            brand_ = html.markup.page()
        # build HTML brand
        if isinstance(brand, html.markup.page):
            brand_.add(str(brand))
        elif brand:
            brand_.div(str(brand), class_='navbar-brand')
        # build state switch
        if len(self.states) > 1:
            statebtn = html.state_switcher(zip(self.states, self.frames), 0)
        else:
            statebtn = False

        # combine and return
        return html.navbar(self._build_nav_links(tabs), brand=brand_,
                           states=statebtn,
                           dropdown_class=['hidden-xs visible-lg', 'hidden-lg'])

    @staticmethod
    def build_html_content(frame, class_='container', id_='main'):
        """Build the #main div for this tab.

        In this construction, the <div id="id_"> is empty, with a
        javascript hook to load the given frame into the div when ready.
        """
        page = html.markup.page()
        page.div('', class_=class_, id_=id_)
        page.add(str(html.load(frame, id_=id_)))
        return page

    def write_state_html(self, state, pre=None, post=None, plots=True):
        """Write the frame HTML for the specific state of this tab

        Parameters
        ----------
        state : `~gwsumm.state.SummaryState`
            `SummaryState` over which to generate inner HTML
        content : `str`, :class:`~gwsumm.html.markup.page`
            simple string content, or a structured `page` of markup to
            embed as the content of the #main div.
        class_ : `str`
            name of HTML class attribute to assign to enclosing div
        id_ : `str`
            name of HTML id attribute to assign to enclosing div
        """
        # build page
        page = html.markup.page()
        if pre:
            page.add(str(pre))
        if plots:
            page.add(str(self.scaffold_plots(state=state)))
        if post:
            page.add(str(post))
        # write to file
        idx = self.states.index(state)
        with open(self.frames[idx], 'w') as fobj:
            fobj.write(str(page))
        return self.frames[idx]

    def write_html(self, title=None, subtitle=None, tabs=list(), ifo=None,
                   ifomap=dict(), brand=None, css=html.CSS, js=html.JS,
                   about=None, footer=None, **inargs):
        """Write the HTML page for this state Tab.

        Parameters
        ----------
        maincontent : `str`, :class:`~gwsumm.html.markup.page`
            simple string content, or a structured `page` of markup to
            embed as the content of the #main div.
        title : `str`, optional, default: {parent.name}
            level 1 heading for this `Tab`.
        subtitle : `str`, optional, default: {self.name}
            level 2 heading for this `Tab`.
        tabs: `list`, optional
            list of top-level tabs (with children) to populate navbar
        ifo : `str`, optional
            prefix for this IFO.
        ifomap : `dict`, optional
            `dict` of (ifo, {base url}) pairs to map to summary pages for
            other IFOs.
        brand : `str`, :class:`~gwsumm.html.markup.page`, optional
            non-menu content for navigation bar, defaults to calendar
        css : `list`, optional
            list of resolvable URLs for CSS files. See `gwsumm.html.CSS` for
            the default list.
        js : `list`, optional
            list of resolvable URLs for javascript files. See
            `gwumm.html.JS` for the default list.
        about : `str`, optional
            href for the 'About' page
        footer : `str`, `~gwsumm.html.markup.page`
            user-defined content for the footer (placed below everything else)
        **inargs
            other keyword arguments to pass to the
            :meth:`~Tab.build_inner_html` method
        """
        return super(StateTab, self).write_html(
            self.frames[0], title=title, subtitle=subtitle, tabs=tabs, ifo=ifo,
            ifomap=ifomap, brand=brand, css=css, js=js, about=about,
            footer=footer, **inargs)

register_tab(StateTab)


class ArchivedStateTab(SummaryArchiveMixin, StateTab):
    """An archivable externally-linked tab.
    """
    type = 'archived-state'

    def __init__(self, name, start, end, mode=None, states=[], **kwargs):
        super(ArchivedStateTab, self).__init__(name, states=states, **kwargs)
        self.span = (start, end)
        self.mode = mode

register_tab(ArchivedStateTab)


class AboutTab(SummaryArchiveMixin, Tab):
    type = 'about'

    def __init__(self, start, end, name='About', mode=None, **kwargs):
        super(AboutTab, self).__init__(name, **kwargs)
        self.span = (start, end)
        self.mode = mode

    def write_html(self, config=[], **kwargs):
        return super(AboutTab, self).write_html(
            html.about_this_page(config=config), **kwargs)

register_tab(AboutTab)
