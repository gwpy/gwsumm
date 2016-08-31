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

"""This module defines the core :class:`Tab` object.
"""

import os
import re
from urlparse import urlparse
from shutil import copyfile

from gwpy.utils.compat import OrderedDict

from gwpy.time import (from_gps, to_gps)
from gwpy.segments import Segment

from .. import (html, mode)
from ..utils import (re_quote, re_cchar)
from ..config import *
from .registry import (get_tab, register_tab)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


class Tab(object):
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
        for all tabs in this run)

    Notes
    -----
    A `Tab` cannot have both a :attr:`~Tab.parent` and :attr:`~tab.Children`.
    This is a limitation imposed by the twitter bootstrap navigation bar
    implementation, which does not allow nested dropdown menus. In order
    to collect child tabs in a given place, assign them all the same
    :attr:`~Tab.group`.
    """
    type = 'basic'
    """Type identifier for this `Tab`"""

    def __init__(self, name, index=None, shortname=None, parent=None,
                 children=list(), group=None, path=os.curdir, hidden=False):
        """Initialise a new `Tab`.
        """
        # names
        self.name = name
        self.shortname = shortname
        # structure
        self._children = []
        self.set_parent(parent)
        self.children = children
        self.group = group
        # HTML format
        self.path = path
        self.index = index
        self.page = None
        self.hidden = hidden

    # -------------------------------------------
    # Tab properties

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
        self.set_parent(p)

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
        """The HTML path (relative to the :attr:`~Tab.path`) for this tab
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
        """HTML href (relative to the :attr:`~Tab.path`) for this tab

        This attribute is just a convenience to clean-up the
        :attr:`~Tab.index` for a given tab, by removing index.htmls.
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
    def from_ini(cls, cp, section, *args, **kwargs):
        """Define a new tab from a :class:`~gwsumm.config..GWConfigParser`

        Parameters
        ----------
        cp : :class:`~gwsumm.config.GWConfigParser`
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
        these as ``*args`` and ``**kwargs`` to this method via :class:`super`:

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
        return cls(name, *args, **kwargs)

    # -------------------------------------------------------------------------
    # HTML operations

    def initialise_html_page(self, title=None, subtitle=None, css=list(),
                             js=list(), copy=True, **initargs):
        """Initialise a new HTML `~markup.page` with a headline.

        This method opens the HTML page with <html> and <body> tags,
        writes a complete <head>, opens the #wrap <div> and writes the
        headline banner.

        Parameters
        ----------
        title : `str`, optional, default: {parent.name}
            level 1 heading for this `Tab`.
        subtitle : `str`, optional, default: {self.name}
            level 2 heading for this `Tab`.
        css : `list`, optional
            list of resolvable URLs for CSS files
        js : `list`, optional
            list of resolvable URLs for javascript files
        copy : `bool`, optional, default: `True`
            copy CSS or javascript files that exist on disk into this
            tab's output directory

        Returns
        -------
        page : :class:`~gwsumm.html.markup.page`
            initialised HTML markup page

        Notes
        -----
        The :meth:`~Tab.finalize_html_page` method is provided to clean
        up the open tags from this method, so if you want to subclass either,
        you should preserve the same tag structure, or subclass both.
        """
        # find relative base path
        base = initargs.pop('base', None)
        if base is None:
            n = len(self.index.split(os.path.sep)) - 1
            base = os.sep.join([os.pardir] * n)
        if not base.endswith('/'):
            base += '/'

        # move css and js files if needed
        if copy:
            staticdir = os.path.join(self.path, 'static')
            for i, script in enumerate(css + js):
                pr = urlparse(script)
                if not pr.netloc and os.path.isfile(script):
                    localscript = os.path.join(staticdir,
                                               os.path.basename(script))
                    if not os.path.isdir(staticdir):
                        os.makedirs(staticdir)
                    if (not os.path.isfile(localscript) or not
                            os.path.samefile(script, localscript)):
                        copyfile(script, localscript)
                        if i >= len(css):
                            js[i-len(css)] = localscript
                        else:
                            css[i] = localscript
        # initialise page
        if title is None:
            title = self.shorttitle
        initargs.setdefault('doctype', html.DOCTYPE)
        initargs.setdefault('metainfo', html.META)
        self.page = html.markup.page()
        self.page.init(css=css, script=js, title=title, base=base, **initargs)
        return self.page

    def finalize_html_page(self, user=True, issues=True, about=None,
                           content=None):
        """Finalize the HTML content for this tab.

        This method closes the #wrap div (opened by
        :meth:`~Tab.initialise_html_page`), prints a footer, and closes
        the <body> and <html> tags to close the HTML content.

        Parameters
        ----------
        user : `bool`, default: `True`
            print details of user running this job
        issues : `bool`, default: `True`
            print link to github.com issue tracker for this package
        about : `str`
            URL for the 'About this page' HTML page
        content : `str`, `~gwsumm.html.markup.page`
            user-defined content for the footer (placed below everything
            else).

        Returns
        -------
        page : :class:`~gwsumm.html.markup.page`
            a copy of the :attr:`~Tab.page` for this tab.

        Notes
        -----
        This method is twinned with the :meth:`~Tab.initialize_html_page`
        method, and is designed to clean up any open tags known to have been
        created from that method. So, if you want to subclass either method,
        you should preserve the same tag structure, or subclass both.
        """
        if not self.page:
            raise RuntimeError("Cannot finalize HTML page before intialising "
                               "one")
        # add footer
        self.page.add(str(html.footer(user=user, issues=issues, about=about,
                                      content=content)))
        self.page.body.close()
        self.page.html.close()

    def build_html_banner(self, title=None, subtitle=None):
        """Build the HTML headline banner for this tab.

        Parameters
        ----------
        title : `str`
            title for this page
        subtitle : `str`
            sub-title for this page

        Returns
        -------
        banner : :class:`~gwsumm.html.markup.page`
            formatter markup page for the banner
        """
        # work title as Parent Name/Tab Name
        if title is None and subtitle is None:
            title = self.title.replace('/', ' : ', 1)
        return html.banner(title, subtitle=subtitle)

    def build_html_navbar(self, brand=None, tabs=list(),
                          ifo=None, ifomap=dict()):
        """Build the navigation bar for this `Tab`.

        Parameters
        ----------
        brand : `str`, :class:`~gwsumm.html.markup.page`
            content for navbar-brand
        tabs : `list`, optional
            list of parent tabs (each with a list of children) to include
            in the navigation bar.
        ifo : `str`, optional
            prefix for this IFO.
        ifomap : `dict`, optional
            `dict` of (ifo, {base url}) pairs to map to summary pages for
            other IFOs.

        Returns
        -------
        page : `~gwsumm.html.markup.page`
            a markup page containing the navigation bar.
        """
        class_ = 'navbar navbar-fixed-top'
        # build interferometer cross-links
        if ifo is not None:
            brand_ = html.base_map_dropdown(ifo, id_='ifos', **ifomap)
            class_ += ' navbar-%s' % ifo.lower()
        else:
            brand_ = html.markup.page()
        # build HTML brand
        if isinstance(brand, html.markup.page):
            brand_.add(str(brand))
        elif brand:
            brand_.div(str(brand), class_='navbar-brand')
        # combine and return
        return html.navbar(self._build_nav_links(tabs), class_=class_,
                           brand=brand_)

    def _build_nav_links(self, tabs):
        """Construct the ordered list of tabs to write into the navbar

        Parameters
        ----------
        tabs : `list`
            a list of `Tabs <gwsumm.tabs.Tab`, some of which may contain
            :attr:`~Tab.children`.

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
                nogroup = sorted(groups.pop(None, []),
                                 key=lambda c: c.shortname.lower()
                                               in ['summary', 'overview'] and
                                               c.shortname.upper() or
                                               c.shortname.lower())
                for child in nogroup:
                    links.append((child.shortname.strip('_'), child.href))
                    if child == self:
                        active = len(links) - 1
                for group in sorted(groups.keys()):
                    # sort group by name
                    re_group = re.compile('(\A{0}\s|\s{0}\Z)'.format(
                                              group.strip('_')), re.I)
                    names = [re_group.sub('', t.shortname)
                             for t in groups[group]]
                    groups[group] = zip(*sorted(
                        zip(groups[group], names),
                        key=lambda (t, n): n.lower() in ['summary', 'overview']
                                           and ' %s' % n.upper() or
                                           n.lower()))[0]
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
    def build_html_content(content):
        """Build the #main div for this tab.

        Parameters
        ----------
        content : `str`, :class:`~gwsumm.html.markup.page`
            HTML content to be wrapped

        Returns
        -------
        #main : :class:`~gwsumm.html.markup.page`
            A new `page` with the input content wrapped as
        """
        page = html.markup.page()
        page.div(id_='main')
        page.div(str(content), id_='content')
        page.div.close()
        return page

    def write_html(self, maincontent, title=None, subtitle=None, tabs=list(),
                   ifo=None, ifomap=dict(), brand=None, base=os.curdir,
                   css=None, js=None, about=None, footer=None, **inargs):
        """Write the HTML page for this `Tab`.

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
            non-menu content for navigation bar
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
        # setup directories
        outdir = os.path.split(self.index)[0]
        if outdir and not os.path.isdir(outdir):
            os.makedirs(outdir)

        # initialise HTML page
        if css is None:
            css = html.get_css().values()
        if js is None:
            js = html.get_js().values()
        self.initialise_html_page(title=title, subtitle=subtitle, css=css,
                                  base=base, js=js)

        # add navigation
        self.page.add(str(self.build_html_navbar(ifo=ifo, ifomap=ifomap,
                                                 tabs=tabs, brand=brand)))
        # -- open main page container
        self.page.div(class_='container')

        # add banner
        self.page.add(str(self.build_html_banner(title=title,
                                                 subtitle=subtitle)))

        # add #main content
        self.page.add(str(self.build_html_content(maincontent)))

        self.page.div.close()  # container
        # close page and write
        self.finalize_html_page(about=about, content=footer)
        with open(self.index, 'w') as fobj:
            fobj.write(str(self.page))
        return

register_tab(Tab)


class SummaryArchiveMixin(object):
    """A `Tab` designed to exist in a calendar-linked archive.
    """
    @property
    def span(self):
        """The GPS [start, end) span of this tab.
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
        """
        try:
            return self.span[0]
        except TypeError:
            return None

    @property
    def end(self):
        """The GPS end time of this tab.
        """
        try:
            return self.span[1]
        except TypeError:
            return None

    @property
    def mode(self):
        """The date-time mode of this tab.

        See Also
        --------
        gwsumm.mode : for details on the modes
        """
        if self._mode is None:
            return mode.get_mode()
        else:
            return self._mode

    @mode.setter
    def mode(self, m):
        if isinstance(m, str):
            m = mode.MODE_ENUM[m]
        if m is None:
            self._mode = m
        else:
            self._mode = int(m)
        return self._mode

    # -------------------------------------------------------------------------
    # SummaryArchive methods

    @classmethod
    def from_ini(cls, config, section, start=None, end=None, hidden=None,
                 **kwargs):
        config = GWSummConfigParser.from_configparser(config)
        if start is None:
            start = config.getint(section, 'gps-start-time')
        if end is None:
            end = config.getint(section, 'gps-end-time')

        return super(SummaryArchiveMixin, cls).from_ini(config, section, start,
                                                        end, **kwargs)

    def build_html_calendar(self):
        """Build the datepicker calendar for this tab.

        Parameters
        ----------
        mode : `int`
            the mode in which to print the calendar, defaults to the selected
            mode in `gwsumm.globalv`.

        Notes
        -----
        The datetime for the calendar is taken from this tab's
        :attr:`~StateTab.span`.
        """
        date = from_gps(self.start)
        # double-check path matches what is required from custom datepicker
        try:
            requiredpath = mode.get_base(date, mode=self.mode)
        except ValueError:
            return html.markup.oneliner.div('%d-%d' % (self.start, self.end),
                                            class_='navbar-brand')
        if not requiredpath in self.path:
            raise RuntimeError("Tab path %r inconsistent with required "
                               "format including %r for archive calendar"
                               % (self.path, requiredpath))
        # format calendar
        return html.calendar(date, mode=self.mode % 3)

    def build_html_navbar(self, brand=None, calendar=True, **kwargs):
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
        brand_ = html.markup.page()
        # add calendar
        if calendar:
            brand_.add(str(self.build_html_calendar()))
        # build HTML brand
        if isinstance(brand, html.markup.page):
            brand_.add(brand)
        elif brand:
            brand_.div(str(brand), class_='navbar-brand')
        # combine and return
        return super(SummaryArchiveMixin, self).build_html_navbar(
            brand=brand_, **kwargs)


# -- TabList -----------------------------------------------------------------

class TabList(list):
    """Custom `list` of `Tab` objects with sorting and parsing
    """

    def __init__(self, entries=[]):
        super(TabList, self).__init__(entries)

    def get_hierarchy(self):
        parents = OrderedDict()
        # 1. Assume all tabs without parents are parents themselves
        for tab in filter(lambda tab: tab.parent is None, self):
            parents[tab.name] = tab
        # 2. All remaining tabs without a defined parent define that parent
        # 3. Sort all tabs into their parent sets
        for tab in filter(lambda tab: tab.parent is not None, self):
            if tab.parent in parents:
                tab.set_parent(parents[tab.parent])
            elif not isinstance(tab.parent, Tab):
                tab.set_parent(get_tab('default')(tab.parent, *tab.span))
            parents.setdefault(tab.parent.name, tab.parent)
            if tab not in tab.parent.children:
                tab.parent.add_child(tab)
        return parents.values()

    @staticmethod
    def _sortkey(tab):
        if tab.shortname == 'Summary' and tab.parent is None:
            return 1
        elif tab.shortname == 'Summary':
            return 2
        elif 'ODC' in tab.shortname:
            return 3
        elif tab.shortname.islower():
            return tab.shortname.upper()
        else:
            return tab.shortname.lower()

    def sort(self, key=None, reverse=False):
        """Sort this `TabList` in place
        """
        if key is None:
            key = self._sortkey
        hlist = self.get_hierarchy()
        hlist.sort(key=key)
        for tab in hlist:
           tab.children.sort(key=key)
        super(TabList, self).sort(key=key, reverse=reverse)

    @classmethod
    def from_ini(cls, config, tag='tab[_-]', match=[],
                 path=os.curdir, plotdir='plots'):
        if isinstance(tag, str):
            tag = re.compile(tag)
        tabs = cls()
        for section in filter(tag.match, config.sections()):
            # if user gave matches, test match and skip
            if match and section[4:] not in match:
                continue
            # otherwise, get type and create instance of class
            try:
                type = config.get(section, 'type')
            except NoOptionError:
                type = 'default'
            else:
                if not type.startswith('archived-'):
                    type = 'archived-%s' % type
            Tab = get_tab(type)
            if issubclass(Tab, get_tab('archived-data')):
                tab = Tab.from_ini(config, section, plotdir=plotdir, path=path)
            else:
                tab = Tab.from_ini(config, section, path=path)
            tabs.append(tab)
        tabs.get_hierarchy()  # call this to resolve map parent names to tabs
        return tabs
