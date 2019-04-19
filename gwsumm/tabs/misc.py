# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013-2016)
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

"""This module defines some utility `Tab` subclasses, including HTTP
error handlers.
"""

from MarkupPy import markup

from ..config import GWSummConfigParser
from .registry import (get_tab, register_tab)
from gwdetchar.io import html

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['AboutTab', 'Error404Tab', 'HTMLContentTab']

Tab = get_tab('basic')


# -- About --------------------------------------------------------------------

class AboutTab(Tab):
    """Page describing how the containing HTML pages were generated
    """
    type = 'about'

    def __init__(self, name='About', **kwargs):
        super(AboutTab, self).__init__(name, **kwargs)

    def write_html(self, config=list(), **kwargs):
        return super(AboutTab, self).write_html(
            html.about_this_page(config=config), **kwargs)


register_tab(AboutTab)


# -- HTTP errors --------------------------------------------------------------

class Error404Tab(Tab):
    """Custom HTTP 404 error page
    """
    type = '404'

    def __init__(self, name='404', **kwargs):
        super(Error404Tab, self).__init__(name, **kwargs)

    def write_html(self, config=list(), top=None, **kwargs):
        if top is None:
            top = kwargs.get('base', self.path)
        kwargs.setdefault('title', '404: Page not found')
        page = markup.page()
        page.div(class_='alert alert-danger')
        page.p()
        page.strong("The page you are looking for doesn't exist")
        page.p.close()
        page.p("This could be because the times for which you are looking "
               "were never processed (or haven't even happened yet), or "
               "because no page exists for the specific data products you "
               "want. Either way, if you think this is in error, please "
               "contact <a class=\"alert-link\" "
               "href=\"mailto:detchar+code@ligo.org\">the DetChar group</a>.")
        page.p("Otherwise, you might be interested in one of the following:")
        page.div(style="padding-top: 10px;")
        page.a("Take me back", role="button", class_="btn btn-lg btn-info",
               title="Back", href="javascript:history.back()")
        page.a("Take me up one level", role="button",
               class_="btn btn-lg btn-warning", title="Up",
               href="javascript:linkUp()")
        page.a("Take me to the top level", role="button",
               class_="btn btn-lg btn-success", title="Top", href=top)
        page.div.close()
        page.div.close()
        page.script("""
  function linkUp() {
    var url = window.location.href;
    if (url.substr(-1) == '/') url = url.substr(0, url.length - 2);
    url = url.split('/');
    url.pop();
    window.location = url.join('/');
  }""", type="text/javascript")
        return super(Error404Tab, self).write_html(page, **kwargs)


register_tab(Error404Tab)


# -- HTMLContent --------------------------------------------------------------
class HTMLContentTab(Tab):
    """A simple tab to add any generic HTML content into the #main div.

    Parameters
    ----------
    name : `str`
        name of this tab (required)
    content : `~MarkupPy.markup.page`, `str`, optional
        content to include in the #main HTML
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
    type = 'htmlcontent'

    def __init__(self, name, content=None, **kwargs):
        """Initialise a new :class:`HTMLContentTab`.
        """
        super(HTMLContentTab, self).__init__(name, **kwargs)

        self.content = content

    @property
    def content(self):
        """HTML content to be included on the page
        """
        return self._content

    @content.setter
    def content(self, content):
        if isinstance(content, markup.page) or content is None:
            self._content = content
        else:
            self._content = markup.page()
            if not str(content).startswith('<'):
                self._content.p(str(content))
            else:
                self._content.add(str(content))

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
        tab : `HTMLContentTab`
            a new tab defined from the configuration
        """
        cp = GWSummConfigParser.from_configparser(cp)

        kwargs.setdefault('path', '')

        # get content
        try:
            kwargs.setdefault('content', cp.get(section, 'content'))
        except NoOptionError:
            pass

        # build and return tab
        return super(HTMLContentTab, cls).from_ini(cp, section, *args, **kwargs)

    def html_content(self, content):
        page = markup.page()
        if self.content:
            page.add(str(self.content))
        if content:
            page.add(str(content))
        return Tab.html_content(str(page))

    def write_html(self, content=None, **kwargs):
        """Write the HTML page for this tab.

        Parameters
        ----------
        content : `str`, :class:`~MarkupPy.markup.page`, optional
            content to place on the page, defaults to
            :attr:`HTMLContentTab.content`
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
        if content is not None:
            self.content = content
        return super(HTMLContentTab, self).write_html(None, **kwargs)


register_tab(HTMLContentTab)
