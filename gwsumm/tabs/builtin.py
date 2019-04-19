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

The `builtin` classes provide interfaces for simple operations including

- `ExternalTab`: embedding an existing webpage
- `PlotTab`: scaffolding a collection of existing images
- `StateTab`: scaffolding plots split into a number of states, with
  a button to switch between states in the HTML

"""

import os.path
import warnings
from configparser import NoOptionError

from six import string_types

from MarkupPy import markup

from .registry import (get_tab, register_tab)
from ..plot import get_plot
from ..utils import (re_quote, re_cchar)
from ..config import GWSummConfigParser
from ..state import (ALLSTATE, SummaryState, get_state)
from .. import html

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['ExternalTab', 'PlotTab', 'StateTab']

Tab = get_tab('basic')
SummaryPlot = get_plot(None)
DataPlot = get_plot('data')


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
    path : `str`
        base output directory for this tab (should be the same directory
        for all tabs in this run).

    Configuration
    -------------

    """
    type = 'external'

    def __init__(self, name, url, error=True, success=None, **kwargs):
        """Initialise a new `ExternalTab`.
        """
        super(ExternalTab, self).__init__(name, **kwargs)
        self.url = url
        self.error = error
        self.success = success

    @property
    def url(self):
        return self._url

    @url.setter
    def url(self, link):
        self._url = link

    @classmethod
    def from_ini(cls, cp, section, *args, **kwargs):
        """Configure a new `ExternalTab` from a `ConfigParser` section

        Parameters
        ----------
        cp : :class:`~gwsumm.config.ConfigParser`
            configuration to parse.
        section : `str`
            name of section to read

        See Also
        --------
        Tab.from_ini :
            for documentation of the standard configuration
            options

        Notes
        -----
        On top of the standard configuration options, the `ExternalTab` can
        be configured with the ``url`` option, specifying the URL of the
        external content to be included:

        .. code-block:: ini

           [tab-external]
           name = External data
           type = external
           url = https://www.example.org/index.html


        """
        url = cp.get(section, 'url')
        if cp.has_option(section, 'error'):
            kwargs.setdefault(
                'error', re_quote.sub('', cp.get(section, 'error')))
        if cp.has_option(section, 'success'):
            kwargs.setdefault(
                'success', re_quote.sub('', cp.get(section, 'success')))
        return super(ExternalTab, cls).from_ini(cp, section, url,
                                                *args, **kwargs)

    def html_content(self, content):
        wrappedcontent = html.load(self.url, id_='content', error=self.error,
                                   success=self.success)
        return super(ExternalTab, self).html_content(wrappedcontent)

    def write_html(self, **kwargs):
        """Write the HTML page for this tab.

        See Also
        --------
        gwsumm.tabs.Tab.write_html : for details of all valid keyword
        arguments
        """
        if not kwargs.pop('writehtml', True):
            return
        link = markup.given_oneliner.a('click here to view the original',
                                       class_='reference',
                                       href=self.url.split()[0])
        kwargs.setdefault('footer', 'This page contains data from an external '
                                    'source, %s.' % link)
        return super(ExternalTab, self).write_html('', **kwargs)


register_tab(ExternalTab)


class PlotTab(Tab):
    """A simple tab to layout some figures in the #main div.

    Parameters
    ----------
    name : `str`
        name of this tab (required)
    plots : `list`, optional
        list of plots to display on this tab. More plots can be added
        at any time via :meth:`PlotTab.add_plot`
    layout : `int`, `list`, optional
        the number of plots to display in each row, or a list of numbers
        to define each row individually. If the number of plots defined
        by the layout is less than the total number of plots, the layout
        for the final row will be repeated as necessary.

        For example ``layout=[1, 2, 3]`` will display a single plot on
        the top row, two plots on the second, and 3 plots on each row
        thereafter.
    foreword : `~MarkupPy.markup.page`, `str`, optional
        content to include in the #main HTML before the plots
    afterword : `~MarkupPy.markup.page`, `str`, optional
        content to include in the #main HTML after the plots
    index : `str`, optional
        HTML file in which to write. By default each tab is written to
        an index.html file in its own directory. Use :attr:`~Tab.index`
        to find out the default index, if not given.
    shortname : `str`, optional
        shorter name for this tab to use in the navigation bar. By
        default the regular name is used
    parent : :class:`~gwsumm.tabs.Tab`, optional
        parent of this tab. This is used to position this tab in the
        navigation bar.
    children : `list`, optional
        list of child :class:`Tabs <~gwsumm.tabs.Tab>` of this one. This
        is used to position this tab in the navigation bar.
    group : `str`, optional
        name of containing group for this tab in the navigation bar
        dropdown menu. This is only relevant if this tab has a parent.
    path : `str`, optional,
        base output directory for this tab (should be the same directory
        for all tabs in this run).
    """
    type = 'plots'

    def __init__(self, name, plots=list(), layout=None, foreword=None,
                 afterword=None, **kwargs):
        """Initialise a new :class:`PlotTab`.
        """
        super(PlotTab, self).__init__(name, **kwargs)
        self.plots = []
        for p in plots:
            self.add_plot(p)
        self.set_layout(layout)
        self.foreword = foreword
        self.afterword = afterword

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
    def layout(self, layout_):
        warnings.warn("Use of `Tab.layout = ...` has been deprecated, "
                      "please switch to using `Tab.set_layout(...)`",
                      DeprecationWarning)
        self.set_layout(layout_)

    def set_layout(self, layout):
        """Set the plot scaffolding layout for this tab

        Parameters
        ----------
        l : `int`, `list` of `int`
            the desired scaffold layout, one of

            - an `int`, indicating the number of plots on every row, or
            - a `list`, indicating the number of plots on each row, with
              the final `int` repeated for any remaining rows; each entry
              should be an `int` or a pair of `int` indicating the
              number of plots on this row AND the desired the scaling of
              this row, see the examples...

        Examples
        --------
        To layout 2 plots on each row

        >>> tab.set_layout(2)

        or

        >>> tab.set_layout([2])

        To layout 2 plots on the first row, and 3 on all other rows

        >>> tab.set_layout((2, 3))

        To layout 2 plots on the first row, and 1 on the second row BUT
        have it the same size as plots on a 2-plot row

        >>> tab.set_layout((2, (1, 2))
        """
        # shortcut None
        if layout is None:
            self._layout = None
            return
        # parse a single int
        if isinstance(layout, int) or (
                isinstance(layout, string_types) and
                layout.isdigit()
        ):
            self._layout = [int(layout)]
            return
        # otherwise parse as a list of ints or pairs of ints
        if isinstance(layout, string_types):
            layout = layout.split(',')
        self._layout = []
        for item in layout:
            if isinstance(item, int):
                self._layout.append(item)
                continue
            if isinstance(item, string_types):
                item = item.strip('([').rstrip(')]').split(',')
            if not isinstance(item, (list, tuple)) or not len(item) == 2:
                raise ValueError("Cannot parse layout element %r (%s)"
                                 % (item, type(item)))
            self._layout.append(tuple(map(int, item)))

    @property
    def foreword(self):
        """HTML content to be included before the plots
        """
        return self._pre

    @foreword.setter
    def foreword(self, content):
        if isinstance(content, markup.page) or content is None:
            self._pre = content
        else:
            self._pre = markup.page()
            if not str(content).startswith('<'):
                self._pre.p(str(content))
            else:
                self._pre.add(str(content))

    @property
    def afterword(self):
        """HTML content to be included after the plots
        """
        return self._post

    @afterword.setter
    def afterword(self, content):
        if isinstance(content, markup.page) or content is None:
            self._post = content
        else:
            self._post = markup.page()
            if not str(content).startswith('<'):
                self._post.p(str(content))
            else:
                self._post.add(str(content))

    @classmethod
    def from_ini(cls, cp, section, *args, **kwargs):
        """Define a new tab from a :class:`~gwsumm.config.GWConfigParser`

        Parameters
        ----------
        cp : :class:`~ConfigParser.GWConfigParser`
            customised configuration parser containing given section
        section : `str`
            name of section to parse

        Returns
        -------
        tab : `PlotTab`
            a new tab defined from the configuration
        """
        cp = GWSummConfigParser.from_configparser(cp)

        kwargs.setdefault('path', '')

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
            for item in layout:
                if isinstance(item, (tuple, list)):
                    item = item[0]
                if item > 12:
                    raise ValueError("Cannot print more than 12 plots in a "
                                     "single row. The chosen layout value for "
                                     "each row must be a divisor of 12 to fit "
                                     "the Bootstrap scaffolding. For details "
                                     "see http://getbootstrap.com/2.3.2/"
                                     "scaffolding.html")
        else:
            layout = None
        kwargs.setdefault('layout', layout)

        # get plots
        if 'plots' not in kwargs:
            kwargs['plots'] = []
            for idx, url in sorted([(int(opt), url) for (opt, url) in
                                    cp.nditems(section) if opt.isdigit()],
                                   key=lambda x: x[0]):
                plot = SummaryPlot(href=url)
                plot.new = False  # this plot does not need to be generated
                # get caption
                try:
                    plot.caption = re_quote.sub(
                        '', cp.get(section, '%d-caption' % idx))
                except NoOptionError:
                    pass
                kwargs['plots'].append(plot)

        # get content
        try:
            kwargs.setdefault('foreword', cp.get(section, 'foreword'))
        except NoOptionError:
            pass
        try:
            kwargs.setdefault('afterword', cp.get(section, 'afterword'))
        except NoOptionError:
            pass

        # build and return tab
        return super(PlotTab, cls).from_ini(cp, section, *args, **kwargs)

    def add_plot(self, plot):
        """Add a plot to this tab.

        Parameters
        ----------
        plot : `str`, :class:`~gwsumm.plot.SummaryPlot`
            either the URL of a plot to embed, or a formatted `SummaryPlot`
            object.
        """
        if isinstance(plot, string_types):
            plot = SummaryPlot(href=plot)
            plot.new = False
        if not isinstance(plot, SummaryPlot):
            raise TypeError("Cannot append plot of type %r" % type(plot))
        self.plots.append(plot)

    def scaffold_plots(self, plots=None, state=None, layout=None,
                       aclass='fancybox plot', **fancyboxargs):
        """Build a grid of plots using bootstrap's scaffolding.

        Returns
        -------
        page : :class:`~MarkupPy.markup.page`
            formatted markup with grid of plots
        """
        page = markup.page()
        page.div(class_='scaffold well')

        if plots is None:
            if state:
                plots = [p for p in self.plots if not
                         isinstance(p, DataPlot) or p.state in [state, None]]
            else:
                plots = self.plots

        # get layout
        if layout is None:
            if self.layout:
                layout = list(self.layout)
            else:
                layout = len(plots) == 1 and [1] or [2]
        for i, item in enumerate(layout):
            if isinstance(item, (list, tuple)):
                layout[i] = item
            elif isinstance(item, int):
                layout[i] = (item, None)
            else:
                raise ValueError("Cannot parse layout element '%s'." % item)
        while sum(list(zip(*layout))[0]) < len(plots):
            layout.append(layout[-1])
        l = i = 0
        fancyboxargs.setdefault('data-fancybox-group', 1)

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
            if plot.src.endswith('svg'):
                fbkw = fancyboxargs.copy()
                fbkw['data-fancybox-type'] = 'iframe'
                page.a(href='%s?iframe' % plot.href.replace('.svg', '.html'),
                       class_=aclass, **fbkw)
            else:
                fbkw = fancyboxargs.copy()
                fbkw['title'] = plot.caption
                page.a(href=plot.href, class_=aclass, **fbkw)
            if plot.src.endswith('.pdf'):
                page.img(class_='img-responsive',
                         src=plot.src.replace('.pdf', '.png'))
            else:
                page.img(class_='img-responsive', src=plot.src)
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

        page.div.close()
        return page

    def html_content(self, content):
        page = markup.page()
        if self.foreword:
            page.add(str(self.foreword))
        if content:
            page.add(str(content))
        page.add(str(self.scaffold_plots()))
        if self.afterword:
            page.add(str(self.afterword))
        return Tab.html_content(str(page))

    def write_html(self, foreword=None, afterword=None, **kwargs):
        """Write the HTML page for this tab.

        Parameters
        ----------
        foreword : `str`, :class:`~MarkupPy.markup.page`, optional
            content to place above the plot grid, defaults to
            :attr:`PlotTab.foreword`
        afterword : `str`, :class:`~MarkupPy.markup.page`, optional
            content to place below the plot grid, defaults to
            :attr:`PlotTab.afterword`
        **kwargs
            other keyword arguments to be passed through
            :meth:`~Tab.write_html`

        See Also
        --------
        gwsumm.tabs.Tab.write_html : for details of all valid unnamed
                                     keyword arguments
        """
        if not kwargs.pop('writehtml', True):
            return
        if foreword is not None:
            self.foreword = foreword
        if afterword is not None:
            self.afterword = self.afterword
        return super(PlotTab, self).write_html(None, **kwargs)


register_tab(PlotTab)


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
    path : `str`
        base output directory for this tab (should be the same directory
        for all tabs in this run).
    """
    type = 'state'

    def __init__(self, name, states=list(), **kwargs):
        """Initialise a new `Tab`.
        """
        if kwargs.get('mode', None) is None:
            raise ValueError("%s() needs keyword argument 'mode'"
                             % type(self).__name__)
        super(StateTab, self).__init__(name, **kwargs)
        # process states
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
    def states(self, statelist, default=None):
        self._states = []
        for state in statelist:
            # allow default indication by trailing asterisk
            if state == default:
                default_ = True
            elif (default is None and isinstance(state, string_types) and
                    state.endswith('*')):
                state = state[:-1]
                default_ = True
            else:
                default_ = False
            self.add_state(state, default=default_)

    def add_state(self, state, default=False):
        """Add a `SummaryState` to this tab

        Parameters
        ----------
        state : `str`, :class:`~gwsumm.state.SummaryState`
            either the name of a state, or a `SummaryState`.
        register : `bool`, default: `False`
            automatically register all new states
        """
        if not isinstance(state, SummaryState):
            state = get_state(state)
        self._states.append(state)
        if default:
            self.defaultstate = state

    @property
    def defaultstate(self):
        try:
            return self._defaultstate
        except AttributeError:
            self._defaultstate = self.states[0]
            return self._defaultstate

    @defaultstate.setter
    def defaultstate(self, state):
        self._defaultstate = state

    @property
    def frames(self):
        # write page for each state
        statelinks = []
        outdir = os.path.split(self.index)[0]
        for i, state in enumerate(self.states):
            statelinks.append(os.path.join(
                outdir, '%s.html' % re_cchar.sub('_', str(state).lower())))
        return statelinks

    # -------------------------------------------
    # StateTab methods

    @classmethod
    def from_ini(cls, cp, section, *args, **kwargs):
        # parse states and retrieve their definitions
        if cp.has_option(section, 'states'):
            # states listed individually
            kwargs.setdefault(
                'states', [re_quote.sub('', s).strip() for s in
                           cp.get(section, 'states').split(',')])
        else:
            # otherwise use 'all' state - full span with no gaps
            kwargs.setdefault('states', ['All'])
        # parse core Tab information
        return super(StateTab, cls).from_ini(cp, section, *args, **kwargs)

    # ------------------------------------------------------------------------
    # HTML methods

    def html_navbar(self, brand='', **kwargs):
        """Build the navigation bar for this `Tab`.

        The navigation bar will consist of a switch for this page linked
        to other interferometer servers, followed by the navbar brand,
        then the full dropdown-based navigation menus configured for the
        given ``tabs`` and their descendents.

        Parameters
        ----------
        brand : `str`, :class:`~MarkupPy.markup.page`
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
        page : `~MarkupPy.markup.page`
            a markup page containing the navigation bar.
        """
        brand_ = markup.page()
        # build state switch
        if len(self.states) > 1 or str(self.states[0]) != ALLSTATE:
            default = self.states.index(self.defaultstate)
            brand_.add(str(html.state_switcher(
                list(zip(self.states, self.frames)), default)))
        # build HTML brand
        if brand:
            brand_.add(str(brand))
        # combine and return
        return super(StateTab, self).html_navbar(brand=brand_, **kwargs)

    @staticmethod
    def html_content(frame):
        r"""Build the #main div for this tab.

        In this construction, the <div id="id\_"> is empty, with a
        javascript hook to load the given frame into the div when ready.
        """
        wrappedcontent = html.load(frame, id_='content')
        return Tab.html_content(str(wrappedcontent))

    def write_state_html(self, state, pre=None, post=None, plots=True):
        """Write the frame HTML for the specific state of this tab

        Parameters
        ----------
        state : `~gwsumm.state.SummaryState`
            `SummaryState` over which to generate inner HTML
        """
        # build page
        page = markup.page()
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
        maincontent : `str`, :class:`~MarkupPy.markup.page`
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
        brand : `str`, :class:`~MarkupPy.markup.page`, optional
            non-menu content for navigation bar, defaults to calendar
        css : `list`, optional
            list of resolvable URLs for CSS files. See `gwsumm.html.CSS` for
            the default list.
        js : `list`, optional
            list of resolvable URLs for javascript files. See
            `gwumm.html.JS` for the default list.
        about : `str`, optional
            href for the 'About' page
        footer : `str`, `~MarkupPy.markup.page`
            user-defined content for the footer (placed below everything else)
        **inargs
            other keyword arguments to pass to the
            :meth:`~Tab.build_inner_html` method
        """
        default = self.states.index(self.defaultstate)
        return super(PlotTab, self).write_html(
            self.frames[default], title=title, subtitle=subtitle,
            tabs=tabs, ifo=ifo, ifomap=ifomap, brand=brand, css=css, js=js,
            about=about, footer=footer, **inargs)


register_tab(StateTab)


class UrlTab(Tab):
    type = 'link'

    def __init__(self, name, url, **kwargs):
        super(UrlTab, self).__init__(name, **kwargs)
        self.href = url

    @property
    def href(self):
        return self._href

    @href.setter
    def href(self, url):
        self._href = url

    @classmethod
    def from_ini(cls, cp, section, *args, **kwargs):
        kwargs.setdefault('url', cp.get(section, 'url'))
        return super(UrlTab, cls).from_ini(cp, section, *args, **kwargs)

    def write_html(self, **kwargs):
        return


register_tab(UrlTab)
