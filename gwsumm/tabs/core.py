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

"""This module defines the core `Tab` object.

The basic Tab allows for simple embedding of arbitrary text inside a
standardised HTML interface. Most real-world applications will use a
sub-class of `Tab` to create more complex HTML output.

The `Tab` class comes in three flavours:

- 'static', no specific GPS time reference
- 'interval', displaying data in a given GPS [start, stop) interval
- 'event', displaying data around a specific central GPS time

The flavour is dynamically set when each instance is created based on the
`mode` keyword, or the presence of `span`, `start and `end`, or `gpstime`
keyword arguments.
"""

import os
import re
from collections import OrderedDict
from configparser import NoOptionError
from shutil import copyfile

from six import (string_types, add_metaclass)
from six.moves.urllib.parse import urlparse

from MarkupPy import markup

from gwpy.time import (from_gps, to_gps)
from gwpy.segments import Segment

from gwdetchar.io import html as gwhtml

from .. import html
from ..mode import (Mode, get_mode, get_base)
from ..utils import (re_quote, re_cchar)
from .registry import (get_tab, register_tab)
from .._version import get_versions

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['BaseTab', 'Tab', 'TabList']


# -- BaseTab ------------------------------------------------------------------
# this object defines the basic object from which all three flavours inherit

class BaseTab(object):
    """The core `Tab` object, defining basic functionality
    """
    def __init__(self, name, index=None,
                 shortname=None, parent=None, children=list(), group=None,
                 path=os.curdir, mode=None, hidden=False):
        # mode
        self.mode = mode
        # names
        self.name = name
        self.shortname = shortname
        # structure
        self._children = []
        self.parent = parent
        self.children = children
        self.group = group
        # HTML format
        self.path = path
        self.index = index
        self.page = None
        self.hidden = hidden

    # -- properties -----------------------------

    @property
    def name(self):
        """Full name for this `Tab`

        :type: `str`
        """
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    @property
    def shortname(self):
        """Short name for this tab

        This will be displayed in the navigation bar.

        :type: `str`
        """
        return self._shortname or self.name

    @shortname.setter
    def shortname(self, name):
        self._shortname = name

    @property
    def parent(self):
        """Short name of the parent page for this `Tab`

        A given tab can either be a parent for a set of child tabs, or can
        have a parent, it cannot be both. In this system, the `parent`
        attribute defines the heading under which this tab will be linked
        in the HTML navigation bar.

        :type: `str`
        """
        return self._parent

    @parent.setter
    def parent(self, p):
        if p is None:
            del self.parent
        else:
            self.set_parent(p)

    @parent.deleter
    def parent(self):
        self._parent = None

    def set_parent(self, p):
        """Set the parent `Tab` for this tab

        Parameters
        ----------
        p : `Tab`
            the parent tab for this one
        """
        if p and self._children:
            raise ValueError("A tab cannot have both a parent, and a "
                             "set of children.")
        if isinstance(p, BaseTab) and self not in p.children:
            p.add_child(self)
        else:
            p = ParentTab(p, [self], mode=self.mode)
        self._parent = p

    @property
    def children(self):
        """List of child tabs for this `Tab`

        If this tab is given children, it cannot also have a parent, as it
        will define its own dropdown menu in the HTML navigation bar, linking
        to itself and its children.

        :type: `list` of `tabs <Tab>`
        """
        return self._children

    @children.setter
    def children(self, clist):
        if self._parent and clist:
            raise ValueError("A Tab cannot have both a parent, and a "
                             "set of children.")
        self._children = list(clist)

    @property
    def index(self):
        """The HTML path (relative to the `~Tab.path`) for this tab
        """
        if not self._index:
            if self.shortname.lower() == 'summary':
                p = ''
            else:
                p = re_cchar.sub('_', self.shortname.strip('_')).lower()
            tab_ = self
            while tab_.parent:
                p = os.path.join(re_cchar.sub(
                        '_', tab_.parent.shortname.strip('_')).lower(), p)
                tab_ = tab_.parent
            self._index = os.path.normpath(os.path.join(
                self.path, p, 'index.html'))
        return self._index

    @index.setter
    def index(self, p):
        self._index = p

    @index.deleter
    def index(self):
        self._index = None

    @property
    def href(self):
        """HTML href (relative to the `~Tab.path`) for this tab

        This attribute is just a convenience to clean-up the
        `~Tab.index` for a given tab, by removing index.htmls.
        hierarchy.

        :type: `str`
        """
        if os.path.basename(self.index) in ('index.html', 'index.php'):
            return os.path.split(self.index)[0] + os.path.sep
        elif self.index:
            return self.index
        else:
            return ''

    @property
    def title(self):
        """Page title for this tab
        """
        if self.parent:
            title = self.name.strip('_')
            tab_ = self
            while tab_.parent:
                title = '%s/%s' % (tab_.parent.name.strip('_'), title)
                tab_ = tab_.parent
            return title
        else:
            return self.name

    @property
    def shorttitle(self):
        """Page title for this tab
        """
        if self.parent:
            title = self.shortname.strip('_')
            tab_ = self
            while tab_.parent:
                title = '%s/%s' % (tab_.parent.shortname.strip('_'), title)
                tab_ = tab_.parent
            return title
        else:
            return self.shortname.strip('_')

    @property
    def group(self):
        """Dropdown group for this `Tab` in the navigation bar

        :type: `str`
        """
        return self._group

    @group.setter
    def group(self, gp):
        if gp is None:
            self._group = None
        else:
            self._group = str(gp)

    @property
    def mode(self):
        """The date-time mode of this tab.

        :type: `int`

        See Also
        --------
        gwsumm.mode : for details on the modes
        """
        return self._mode

    @mode.setter
    def mode(self, m):
        self._mode = get_mode(m)

    @mode.deleter
    def mode(self):
        self._mode = get_mode(0)

    # -- Tab instance methods -------------------

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

    # -- Tab configuration parser ---------------

    @classmethod
    def from_ini(cls, cp, section, *args, **kwargs):
        """Define a new tab from a `~gwsumm.config.GWConfigParser`

        Parameters
        ----------
        cp : `~gwsumm.config.GWConfigParser`
            customised configuration parser containing given section

        section : `str`
            name of section to parse

        *args, **kwargs
            other positional and keyword arguments to pass to the class
            constructor (`__init__`)

        Returns
        -------
        tab : `Tab`
            a new tab defined from the configuration

        Notes
        -----
        This method parses the following configuration options

        .. autosummary::

           ~Tab.name
           ~Tab.shortname
           ~Tab.parent
           ~Tab.group
           ~Tab.index

        Sub-classes should parse their own configuration values and then pass
        these as ``*args`` and ``**kwargs`` to this method via `super`:

        .. code-block:: python

           class MyTab(Tab):
               [...]
               def from_ini(cls, cp, section)
                   \"\"\"Define a new `MyTab`.
                   \"\"\"
                   foo = cp.get(section, 'foo')
                   bar = cp.get(section, 'bar')
                   return super(MyTab, cls).from_ini(cp, section, foo, bar=bar)
        """
        # get tab name
        try:
            # name given explicitly
            name = re_quote.sub('', cp.get(section, 'name'))
        except NoOptionError:
            # otherwise strip 'tab-' from section name
            name = section[4:]
        try:
            kwargs.setdefault('shortname',
                              re_quote.sub('', cp.get(section, 'shortname')))
        except NoOptionError:
            pass
        # get parent:
        #     if parent is not given, this assumes a top-level tab
        try:
            kwargs.setdefault('parent',
                              re_quote.sub('', cp.get(section, 'parent')))
        except NoOptionError:
            pass
        else:
            if kwargs['parent'] == 'None':
                kwargs['parent'] = None
        # get group
        try:
            kwargs.setdefault('group', cp.get(section, 'group'))
        except NoOptionError:
            pass
        # get HTML file
        try:
            kwargs.setdefault('index', cp.get(section, 'index'))
        except NoOptionError:
            pass
        # get hidden param
        try:
            hidden = cp.get(section, 'hidden')
        except NoOptionError:
            hidden = False
        else:
            if hidden is None:
                hidden = True
            else:
                hidden = bool(hidden.title())
        kwargs.setdefault('hidden', hidden)

        # get mode and times if required
        try:
            kwargs['mode']
        except KeyError:
            try:
                kwargs['mode'] = get_mode(cp.get(section, 'mode'))
            except NoOptionError:
                kwargs['mode'] = get_mode()
        if kwargs['mode'] >= Mode.gps:
            try:
                kwargs['start']
            except KeyError:
                kwargs['start'] = cp.getint(section, 'gps-start-time')
            try:
                kwargs['end']
            except KeyError:
                kwargs['end'] = cp.getint(section, 'gps-end-time')
        elif kwargs['mode'] == Mode.event:
            try:
                kwargs['gpstime']
            except KeyError:
                kwargs['gpstime'] = cp.getfloat(section, 'gpstime')
            try:
                kwargs['duration']
            except KeyError:
                kwargs['duration'] = cp.getfloat(section, 'duration')

        return cls(name, *args, **kwargs)

    # -- HTML operations ------------------------
    # the following related HTML operations are defined here
    #
    # - navbar: create navbar
    # - banner: create header
    # - content: create main content
    #
    # The `Tab.write_html` method pulls all of these things together and
    # is the primary user-facing HTML method

    def html_banner(self, title=None, subtitle=None):
        """Build the HTML headline banner for this tab.

        Parameters
        ----------
        title : `str`
            title for this page

        subtitle : `str`
            sub-title for this page

        Returns
        -------
        banner : `~MarkupPy.markup.page`
            formatter markup page for the banner
        """
        # work title as Parent Name/Tab Name
        if title is None and subtitle is None:
            title = self.title.replace('/', ' : ', 1)
        return html.banner(title, subtitle=subtitle)

    def html_navbar(self, brand=None, tabs=list(), ifo=None, ifomap=dict(),
                    **kwargs):
        """Build the navigation bar for this tab.

        Parameters
        ----------
        brand : `str`, `~MarkupPy.markup.page`
            content to place inside `<div class="navbar-brand"></div>`

        tabs : `list`, optional
            list of parent tabs (each with a list of children) to include
            in the navigation bar.

        ifo : `str`, optional
            prefix for this IFO.

        ifomap : `dict`, optional
            `dict` of (ifo, {base url}) pairs to map to summary pages for
            other IFOs.

        **kwargs
            other keyword arguments to pass to :meth:`gwsumm.html.navbar`

        Returns
        -------
        page : `~MarkupPy.markup.page`
            a markup page containing the navigation bar.
        """
        class_ = 'navbar navbar-fixed-top'
        # build interferometer cross-links
        if ifo is not None:
            brand_ = html.base_map_dropdown(ifo, id_='ifos', bases=ifomap)
            class_ += ' navbar-%s' % ifo.lower()
        else:
            brand_ = markup.page()
        # build HTML brand
        if brand:
            brand_.add(str(brand))
        # combine and return
        return gwhtml.navbar(self._html_navbar_links(tabs), class_=class_,
                             brand=brand_, **kwargs)

    def _html_navbar_links(self, tabs):
        """Construct the ordered list of tabs to write into the navbar

        Parameters
        ----------
        tabs : `list`
            a list of `Tabs <gwsumm.tabs.Tab`, some of which may contain
            `~Tab.children`.

        Returns
        -------
        links : `str`
            a structured list of navigation menu entries. Each element will
            be a (`name`, `entries`) tuple defining the heading and dropdown
            entries for a single dropdown menu. The `entries` list will again
            be a list of (`name`, (`link`|`entries`)) tuples, each defining
            a name and link for a given entry in a single dropdown menu,
            or a set of (`name`, `link`) tuples for a group in a dropdown
            menu.
        """
        # build navbar links
        navlinks = []
        tabs = TabList(tabs).get_hierarchy()
        for tab in tabs:
            if tab.hidden:
                continue
            children = [t for t in tab.children if not t.hidden]
            if len(children):
                navlinks.append([tab.shortname.strip('_'), []])
                links = []
                active = None
                # build groups
                groups = set([t.group for t in children])
                groups = dict((g, [t for t in children if t.group == g])
                              for g in groups)
                nogroup = sorted(
                    groups.pop(None, []),
                    key=(lambda c:
                         c.shortname.lower() in ['summary', 'overview'] and
                         c.shortname.upper() or c.shortname.lower()))
                for child in nogroup:
                    links.append((child.shortname.strip('_'), child.href))
                    if child == self:
                        active = len(links) - 1
                for group in sorted(list(groups)):
                    # sort group by name
                    re_group = re.compile(
                        r'(\A{0}\s|\s{0}\Z)'.format(group.strip('_')),
                        re.I,
                    )
                    names = [re_group.sub('', t.shortname)
                             for t in groups[group]]
                    groups[group] = list(zip(*sorted(
                        list(zip(groups[group], names)),
                        key=(lambda x:
                             x[1].lower() in ['summary', 'overview'] and
                             ' %s' % x[1].upper() or x[1].lower()),
                    )))[0]
                    # build link sets
                    links.append((group.strip('_'), []))
                    for i, child in enumerate(groups[group]):
                        name = re_group.sub('', child.shortname.strip('_'))
                        links[-1][1].append((name, child.href))
                        if child == self:
                            active = [len(links) - 1, i]
                if (children[0].shortname == 'Summary' and
                        not children[0].group and len(children)) > 1:
                    links.insert(1, None)
                    if active and isinstance(active, int) and active > 0:
                        active += 1
                    elif active and isinstance(active, list) and active[0] > 0:
                        active[0] += 1
                navlinks[-1][1].extend(links)
                navlinks[-1].append(active)
            else:
                navlinks.append((tab.shortname.strip('_'), tab.href))
        return navlinks

    @staticmethod
    def html_content(content):
        """Build the #main div for this tab.

        Parameters
        ----------
        content : `str`, `~MarkupPy.markup.page`
            HTML content to be wrapped

        Returns
        -------
        #main : `~MarkupPy.markup.page`
            A new `page` with the input content wrapped as
        """
        page = markup.page()
        page.div(id_='main')
        page.div(str(content), id_='content')
        page.div.close()
        return page

    def write_html(self, maincontent, title=None, subtitle=None, tabs=list(),
                   ifo=None, ifomap=dict(), brand=None, base=None,
                   css=None, js=None, about=None, footer=None, issues=True,
                   **inargs):
        """Write the HTML page for this `Tab`.

        Parameters
        ----------
        maincontent : `str`, `~MarkupPy.markup.page`
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

        brand : `str`, `~MarkupPy.markup.page`, optional
            non-menu content for navigation bar

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

        issues : `bool` or `str`, default: `True`
            print link to github.com issue tracker for this package

        **inargs
            other keyword arguments to pass to the
            :meth:`~Tab.build_inner_html` method
        """
        # setup directories
        outdir = os.path.split(self.index)[0]
        if outdir and not os.path.isdir(outdir):
            os.makedirs(outdir)

        # get default style and scripts
        if css is None:
            css = list(html.get_css().values())
        if js is None:
            js = list(html.get_js().values())

        # find relative base path
        if base is None:
            n = len(self.index.split(os.path.sep)) - 1
            base = os.sep.join([os.pardir] * n)
        if not base.endswith('/'):
            base += '/'

        # set default title
        if title is None:
            title = self.shorttitle.replace('/', ' | ')

        # construct navigation
        if tabs:
            navbar = str(self.html_navbar(ifo=ifo, ifomap=ifomap,
                                          tabs=tabs, brand=brand))

        # initialize page
        self.page = gwhtml.new_bootstrap_page(
            title=title, base=base, css=css, script=js, navbar=navbar)

        # add banner
        self.page.add(str(self.html_banner(
            title=title.replace('|', ':'), subtitle=subtitle)))

        # add #main content
        self.page.add(str(self.html_content(maincontent)))

        # format custom footer
        version = get_versions()['version']
        commit = get_versions()['full-revisionid']
        url = 'https://github.com/gwpy/gwsumm/tree/{}'.format(commit)
        link = markup.oneliner.a(
            'View gwsumm-{} on GitHub'.format(version),
            href=url, target='_blank')
        issues = markup.oneliner.a(
            'Report an issue', href=issues, target='_blank')

        # close page and write
        gwhtml.close_page(self.page, self.index, about=about,
                          link=link, issues=issues, content=footer)
        return


# -- Mixins -------------------------------------------------------------------
#
# All actual `Tab` objects will come in three favours:
#
#     - `StaticTab` - no GPS associations
#     - `IntervalTab` - associated with GPS start and stop time
#     - `EventTab` - associated with central GPS time (and duration)


class StaticTab(BaseTab):
    """Simple `Tab` with no GPS association
    """
    pass


class GpsTab(BaseTab):
    """Stub for GPS-related tabs
    """
    @property
    def span(self):
        """The GPS [start, end) span of this tab.

        :type: `~gwpy.segments.Segment`
        """
        return self._span

    @span.setter
    def span(self, seg):
        if seg:
            self._span = Segment(*map(to_gps, seg))
        else:
            self._span = None

    @property
    def start(self):
        """The GPS start time of this tab.

        :type: `float`
        """
        try:
            return self.span[0]
        except TypeError:
            return None

    @property
    def end(self):
        """The GPS end time of this tab.

        :type: `float`
        """
        try:
            return self.span[1]
        except TypeError:
            return None


class IntervalTab(GpsTab):
    """`Tab` defined within a GPS [start, end) interval
    """
    def __init__(self, *args, **kwargs):
        try:
            span = kwargs.pop('span')
        except KeyError:
            try:
                start = kwargs.pop('start')
                end = kwargs.pop('end')
            except KeyError:
                mode = get_mode(kwargs.get('mode')).name
                raise TypeError("%s() in %r mode needs keyword argument 'span'"
                                " or both 'start' and 'end'"
                                % (type(self).__name__, mode))
            else:
                span = (start, end)
        self.span = span
        super(IntervalTab, self).__init__(*args, **kwargs)

    def html_calendar(self):
        """Build the datepicker calendar for this tab.

        Notes
        -----
        The datetime for the calendar is taken from this tab's `~GpsTab.span`
        """
        date = from_gps(self.start)
        # double-check path matches what is required from custom datepicker
        try:
            requiredpath = get_base(date, mode=self.mode)
        except ValueError:
            return markup.oneliner.div('%d-%d' % (self.start, self.end),
                                       class_='navbar-brand')
        if requiredpath not in self.path:
            raise RuntimeError("Tab path %r inconsistent with required "
                               "format including %r for archive calendar"
                               % (self.path, requiredpath))
        # format calendar
        return html.calendar(date, mode=self.mode)

    def html_navbar(self, brand=None, calendar=True, **kwargs):
        """Build the navigation bar for this `Tab`.

        The navigation bar will consist of a switch for this page linked
        to other interferometer servers, followed by the navbar brand,
        then the full dropdown-based navigation menus configured for the
        given ``tabs`` and their descendents.

        Parameters
        ----------
        brand : `str`, `~MarkupPy.markup.page`
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
        # build interferometer cross-links
        brand_ = markup.page()
        # add calendar
        if calendar:
            brand_.add(str(self.html_calendar()))
        # build HTML brand
        if brand:
            brand_.add(str(brand))
        # combine and return
        return super(IntervalTab, self).html_navbar(brand=brand_, **kwargs)


class EventTab(GpsTab):
    """`Tab` defined around a central GPS time
    """
    @property
    def gpstime(self):
        """Central GPS time of this tab

        :type: `~gwpy.time.LIGOTimeGPS`
        """
        return self._gpstime

    @gpstime.setter
    def gpstime(self, t):
        self._gpstime = to_gps(t)

    @property
    def datetime(self):
        return from_gps(self.gpstime)

    @property
    def duration(self):
        """Time duration of this tab, centred on the `gpstime`

        :type: `float`
        """
        return abs(self.span)

    @duration.setter
    def duration(self, d):
        d2 = d/2.
        self.span = (self.gpstime - d2, self.gpstime + d2)

    # -- EventTab methods -----------------------

    def __init__(self, *args, **kwargs):
        # parse gpstime and duration
        try:
            gpstime = kwargs.pop('gpstime')
        except KeyError:
            mode = get_mode(kwargs['mode']).name
            raise TypeError("%s() in %r mode needs keyword argument 'gpstime'"
                            % (type(self).__name__, mode))
        duration = kwargs.pop('duration', 200)
        self.gpstime = gpstime
        self.duration = duration
        # create tab and assign properties
        super(EventTab, self).__init__(*args, **kwargs)

    def html_navbar(self, brand=None, **kwargs):
        if brand is None:
            brand = str(self.gpstime)
        super(EventTab, self).html_navbar(brand=brand, **kwargs)
    html_navbar.__doc__ = GpsTab.html_navbar.__doc__


class _MetaTab(type):
    """Metaclass for creating tabs of the right flavour

    The `MetaTab.__call__` method will get executed whenever a subclass of
    `Tab` is created. It works out the 'mode' of the new tab and dynamically
    sets the parent class for `Tab` accordingly.
    """
    def __call__(cls, *args, **kwargs):
        """Parse the `mode` kwarg for the Tab and add the right flavour
        """
        # parse default mode based on other kwargs
        try:
            mode = kwargs['mode']
        except KeyError:
            if 'gpstime' in kwargs:
                mode = 'EVENT'
            else:
                mode = None
        mode = get_mode(mode)
        # parse regular Tab (don't add a mixin)
        if mode == Mode.static:
            kwargs.pop('mode', None)
            base = StaticTab
        # parse event Tab
        elif mode == Mode.event:
            kwargs['mode'] = mode
            base = EventTab
        # parse interval Tab
        else:
            kwargs['mode'] = mode
            base = IntervalTab
        # set bases and create Tab
        Tab.__bases__ = (base,)
        return super(_MetaTab, cls).__call__(*args, **kwargs)


# -- Tab ----------------------------------------------------------------------
# this is the first actual Tab object, all of the functionality is defined
# in the `BaseTab` object

@add_metaclass(_MetaTab)
class Tab(BaseTab):
    """A Simple HTML tab.

    This `class` provides a mechanism to generate a full-formatted
    HTML page including banner, navigation-bar, content, and a footer,
    without the user worrying too much about the details.

    For example::

    >>> # import Tab and make a new one with a given title and HTML file
    >>> from gwsumm.tabs import Tab
    >>> tab = Tab('My new tab', 'mytab.html')
    >>> # write the Tab to disk with some simple content
    >>> tab.write_html('This is my content', brand='Brand name')

    Parameters
    ----------
    name : `str`
        name of this tab (required)

    index : `str`
        HTML file in which to write. By default each tab is written to
        an index.html file in its own directory. Use `~Tab.index`
        to find out the default index, if not given.

    shortname : `str`
        shorter name for this tab to use in the navigation bar. By
        default the regular name is used

    parent : `~gwsumm.tabs.Tab`
        parent of this tab. This is used to position this tab in the
        navigation bar.

    children : `list`
        list of child `Tabs <~gwsumm.tabs.Tab>` of this one. This
        is used to position this tab in the navigation bar.

    group : `str`
        name of containing group for this tab in the navigation bar
        dropdown menu. This is only relevant if this tab has a parent.

    path : `str`
        base output directory for this tab (should be the same directory
        for all tabs in this run)

    Notes
    -----
    A `Tab` cannot have both a `~Tab.parent` and `~tab.Children`.
    This is a limitation imposed by the twitter bootstrap navigation bar
    implementation, which does not allow nested dropdown menus. In order
    to collect child tabs in a given place, assign them all the same
    `~Tab.group`.
    """
    type = 'basic'


register_tab(Tab)


class ParentTab(Tab):
    """Dummy `Tab` only for navigation
    """
    def __init__(self, name, children, **kwargs):
        # parse list of children
        if not isinstance(children, list):
            children = [children]
        child = children[0]
        # parse mode and GPS arguments (if required)
        kwargs.setdefault('mode', child.mode)
        if isinstance(self, EventTab):
            kwargs.setdefault('gpstime', child.gpstime)
            kwargs.setdefault('duration', child.duration)
        elif isinstance(self, IntervalTab):
            kwargs.setdefault('span', child.span)
        # create Tab
        super(ParentTab, self).__init__(name, children=children, **kwargs)


# -- TabList -----------------------------------------------------------------

class TabList(list):
    """Custom `list` of `Tab` objects with sorting and parsing
    """

    def __init__(self, entries=[]):
        super(TabList, self).__init__(entries)

    def get_hierarchy(self):
        parents = OrderedDict()
        # 1. Assume all tabs without parents are parents themselves
        for tab in [tab for tab in self if tab.parent is None]:
            parents[tab.name] = tab
        # 2. All remaining tabs without a defined parent define that parent
        # 3. Sort all tabs into their parent sets
        for tab in [tab for tab in self if tab.parent is not None]:
            if tab.parent in parents:
                tab.set_parent(parents[tab.parent])
            elif not isinstance(tab.parent, Tab):
                tab.set_parent(get_tab('default')(
                    tab.parent, mode=tab.mode, span=tab.span))
            parents.setdefault(tab.parent.name, tab.parent)
            if tab not in tab.parent.children:
                tab.parent.add_child(tab)
        return list(parents.values())

    @staticmethod
    def _sortkey(tab):
        # NOTE: we need all return values to be strings for
        #       the sorting to actually work
        if 'Home' in tab.shortname:
            return '1'
        if tab.shortname == 'Summary' and tab.parent is None:
            return '2'
        if tab.shortname == 'Summary':
            return '3'
        if 'ODC' in tab.shortname:
            return '4'
        if tab.shortname.islower():
            return tab.shortname.upper()
        return tab.shortname.lower()

    def sort(self, key=None, reverse=False):
        """Sort this `TabList` in place
        """
        if key is None:
            key = self._sortkey
        hlist = sorted(self.get_hierarchy(), key=key)
        for tab in hlist:
            tab.children.sort(key=key)
        super(TabList, self).sort(key=key, reverse=reverse)

    @classmethod
    def from_ini(cls, config, tag='tab[_-]', match=[],
                 path=os.curdir, plotdir='plots'):
        if isinstance(tag, string_types):
            tag = re.compile(tag)
        tabs = cls()
        parents = {}
        for section in filter(tag.match, config.sections()):
            # if user gave matches, test match and skip
            if match and section[4:] not in match:
                continue
            # otherwise, get type and create instance of class
            try:
                type_ = config.get(section, 'type')
            except NoOptionError:
                type_ = 'default'
            Tab = get_tab(type_)
            if issubclass(Tab, get_tab('data')):
                tab = Tab.from_ini(config, section, plotdir=plotdir, path=path)
            else:
                tab = Tab.from_ini(config, section, path=path)
            tabs.append(tab)
            if tab.parent and tab.parent.name in parents:
                tab.set_parent(parents[tab.parent.name])
            elif tab.parent:
                parents[tab.parent.name] = tab.parent
        tabs.get_hierarchy()  # call this to resolve map parent names to tabs
        return tabs
