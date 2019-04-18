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

from .registry import (get_tab, register_tab)
from gwdetchar.io import html

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__all__ = ['AboutTab', 'Error404Tab']

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
