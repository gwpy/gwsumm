# -*- coding: utf-8 -*-
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

"""Helper functions for twitter-bootstrap HTML constructs.
"""

import sys
import getpass
import datetime
import os.path

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

from .. import markup
from ..utils import highlight_syntax
from ...mode import (Mode, get_mode)
from ...utils import re_cchar
from ..._version import get_versions

# set resources
BOOTSTRAP_VERSION = "3.3.6"
BOOTSTRAP_CDN = '//netdna.bootstrapcdn.com/bootstrap/%s' % BOOTSTRAP_VERSION
BOOTSTRAP_CSS = '%s/css/bootstrap.min.css' % BOOTSTRAP_CDN
BOOTSTRAP_JS = '%s/js/bootstrap.min.js' % BOOTSTRAP_CDN

# set <meta> for bootstrap
META = {'viewport': 'width=device-width, initial-scale=1.0'}

# -----------------------------------------------------------------------------
# global HTML constructs

NAVBAR_TOGGLE = """<button class="navbar-toggle" data-toggle="collapse" type="button" data-target=".navbar-collapse">
  <span class="icon-bar"></span>
  <span class="icon-bar"></span>
  <span class="icon-bar"></span>
</button>"""


# -----------------------------------------------------------------------------
# variable HTML constructs

def banner(title, subtitle=None, titleclass=None, substitleclass=None):
    """Construct a banner heading in bootstrap format

    Parameters
    ----------
    title : `str`
        name of page (<h1>)
    subtitle : `str`, optional
        description of page (<p>)
    titleclass : `str`, optional
        class option for <h1>
    subtitleclass : `str`, optional
        class option for <p>

    Returns
    -------
    banner : `~gwsumm.html.markup.page`
        markup.py `page` instance
    """
    page = markup.page()
    page.div(class_='banner')
    if titleclass is None:
        page.h1(str(title))
    else:
        page.h1(str(title), class_=titleclass)
    if subtitle is not None and subtitleclass is not None:
        page.p(subtitle, class_=subtitleclass)
    elif subtitle is not None:
        page.p(subtitle)
    page.div.close()

    return page


def navbar(links, class_='navbar navbar-fixed-top', brand=None, collapse=True):
    """Construct a navigation bar in bootstrap format.

    Parameters
    ----------
    links : `list`
        list of either (Link text, linkurl) tuples or
        :class:`~gwsumm.html.markup.page` objects. Tuples will be
        written in ``<a>`` tags and `pages` will be copied directly.
        Both will be enclosed in a <li> element inside the navbar

    Returns
    -------
    page : :class:`~gwsumm.html.markup.page`
        stand-alone navbar HTML `page`
    """
    # set up page
    page = markup.page()
    markup.element('header', parent=page)(class_=class_)
    page.div(class_="container")

    # ---- non-collapable part (<div class="navbar-header">) ----

    page.div(class_="navbar-header")

    # add collapsed menu toggle
    if collapse:
        page.add(NAVBAR_TOGGLE)

    # add branding (generic non-collapsed content)
    if brand:
        page.add(str(brand))

    page.div.close()  # navbar-header

    if collapse:
        page.nav(class_="collapse navbar-collapse")
    else:
        page.nav()

    # ---- collapsable part (<nav>) ----

    if links:
        page.ul(class_='nav navbar-nav')
        for i, link in enumerate(links):
            if (isinstance(link, (list, tuple)) and
                    isinstance(link[1], basestring)):
                page.li()
                text, link = link
                page.a(text, href=link)
            elif (isinstance(link, (list, tuple)) and
                  isinstance(link[1], (list, tuple))):
                page.li(class_='dropdown')
                page.add(str(dropdown(*link)))
            else:
                page.li()
                page.add(str(link))
            page.li.close()
        page.ul.close()

    page.nav.close()
    page.div.close()
    markup.element('header', parent=page).close()
    return page


def dropdown(text, links, active=None, class_='dropdown-toggle'):
    """Construct a dropdown menu in bootstrap format

    Parameters
    ----------
    text : `str`
        dropdown menu header
    links : `list`
        list of (Link text, linkurl) tuples or dict of such tuples for
        grouped dropdowns

    Returns
    -------
    page : :class:`~gwsumm.html.markup.page`
        HTML element with the following grammar::

        .. code:: html

            <a>text</a>
            <ul>
                <li>link</li>
                <li>link</li>
            </ul>
    """
    page = markup.page()
    # dropdown header
    page.a(href='#', class_=class_, **{'data-toggle': 'dropdown'})
    page.add(text)
    page.b('', class_='caret')
    page.a.close()

    # work out columns
    ngroup = sum([isinstance(x, (tuple, list)) and len(x) and
                 isinstance(x[1], (tuple, list)) for x in links])
    if ngroup < 1:
        column = ''
    else:
        ncol = min(ngroup, 4)
        column = 'col-xs-12 col-sm-%d' % (12 // ncol)

    # dropdown elements
    if column:
        page.ul(class_='dropdown-menu dropdown-%d-col row' % ncol)
    else:
        page.ul(class_='dropdown-menu')
    for i, link in enumerate(links):
        if isinstance(active, int) and i == active:
            active_ = True
        elif isinstance(active, (list, tuple)) and i == active[0]:
            active_ = active[1]
        else:
            active_ = False
        dropdown_link(page, link, active=active_, class_=column)
    page.ul.close()
    return page


def dropdown_link(page, link, active=False, class_=''):
    if link is None:
        page.li(class_='divider')
    elif active is True:
        page.li(class_='active')
    else:
        page.li()
    if isinstance(link, (tuple, list)):
        if isinstance(link[1], (tuple, list)):
            page.ul(class_=class_ + ' list-unstyled')
            page.li(link[0], class_='dropdown-header')
            for j, link2 in enumerate(link[1]):
                dropdown_link(page, link2,
                              active=(type(active) is int and active == j))
            page.ul.close()
        else:
            page.a(link[0], href=link[1])

    elif link is not None:
        page.add(str(link))
    page.li.close()


def calendar(date, tag='a', class_='navbar-brand dropdown-toggle',
             id_='calendar', dateformat=None, mode=None):
    """Construct a bootstrap-datepicker calendar.

    Parameters
    ----------
    date : :class:`datetime.datetime`, :class:`datetime.date`
        active date for the calendar
    tag : `str`
        type of enclosing HTML tag, default: ``<a>``

    Returns
    -------
    calendar : `str`
        a onliner string of HTML containing the calendar text and a
        triggering dropdown
    """
    mode = get_mode(mode)
    if mode < Mode.day:
        raise ValueError("Cannot generate calendar for Mode %s" % mode)
    if dateformat is None:
        if mode == Mode.day:
            dateformat = '%B %d %Y'
        elif mode == Mode.week:
            dateformat = 'Week of %B %d %Y'
        elif mode == Mode.month:
            dateformat = '%B %Y'
        elif mode == Mode.year:
            dateformat = '%Y'
        else:
            raise ValueError("Cannot set dateformat for Mode %s" % mode)
    datestring = date.strftime(dateformat).replace(' 0', ' ')
    data_date = date.strftime('%d-%m-%Y')
    page = markup.page()
    page.a('&laquo;', class_='navbar-brand step-back', title='Step back',
           onclick='stepDate(-1)')
    page.a(id_=id_, class_=class_, title='Show/hide calendar',
           **{'data-date': data_date, 'data-date-format': 'dd-mm-yyyy',
              'data-viewmode': '%ss' % mode.name})
    page.add(datestring)
    page.b('', class_='caret')
    page.a.close()
    page.a('&raquo;', class_='navbar-brand step-forward',
           title='Step forwards', onclick='stepDate(1)')
    return page


def wrap_content(page):
    """Utility to wrap some HTML into the relevant <div>s for the main
    body of a page in bootstrap format

    Parameters
    ----------
    page : :class:`~gwsumm.html.markup.page`, `str`
        HTML content to be wrapped
    span : `int`
        column span of content, default: 'full' (``12``)

    Returns
    -------
    wrappedpage : :class:`~gwsumm.html.markup.page`
        A new `page` with the input content wrapped as

        .. code:: html

            <div class="container">
            </div>
    """
    out = markup.page()
    out.div(class_='container', id_='main')
    out.add(str(page))
    out.div.close()
    return out


def footer(user=True, about=None, issues=True, content=None, span=12,
           class_='footer'):
    """Construct a footer with the given HTML content

    Parameters
    ----------
    user : `bool`, optional, default: `True`
        print details of user running this job
    span : `int`
        column span of content, default: 'full' (``12``)

    Returns
    -------
    page : :class:`~gwsumm.html.markup.page`
        HTML footer around content as

        .. code:: html

            <footer>
                <div class="container">
                    <div class="col-md-{span}">
                        content
                    </div>
                </div>
            </footer>
    """
    page = markup.page(twotags=['footer'])
    page.FOOTER(class_=class_)
    page.div(class_='container')
    page.div(class_='row')
    page.div(class_='col-md-%d' % span)
    if user:
        page.p('This page was created by %s at %s.'
               % (getpass.getuser(),
                  datetime.datetime.now().strftime('%H:%m on %B %d %Y')))
    if issues:
        version = get_versions()['version']
        commit = get_versions()['full-revisionid']
        url = 'https://github.com/gwpy/gwsumm/tree/%s' % commit
        page.p()
        page.a('View GWSumm %s on GitHub' % version, href=url, target='_blank')
        page.add('|')
        page.a('Report an issue', href='https://github.com/gwpy/gwsumm/issues',
               target='_blank')
        page.p.close()
    if about is not None:
        page.p(markup.oneliner.a('How was this page generated?', href=about))
    if isinstance(content, markup.page):
        page.add(str(content))
    elif content is not None:
        page.p(str(content))
    page.div.close()
    page.div.close()
    page.div.close()
    page.FOOTER.close()
    return page


def state_switcher(states, default=0):
    """Build a state switch button, including all of the given
    states, with the default selected by index
    """
    current, chref = states[default]
    page = markup.page()
    page.div(class_="btn-group pull-right state-switch")
    page.a(class_='navbar-brand dropdown-toggle', href='#', id_='states',
           title="Show/hide state menu", **{'data-toggle': 'dropdown'})
    page.add(str(current))
    page.b('', class_='caret')
    page.a.close()
    page.ul(class_='dropdown-menu', id_='statemenu')
    page.li("Select an option below to view these data in another state "
            "(different time segments).", class_="dropdown-header")
    page.li('', class_="divider")
    for i, (state, href) in enumerate(states):
        page.li()
        page.a(str(state), class_='state', title=str(state),
               id_='state_%s' % re_cchar.sub('_', str(state)).lower(),
               onclick='$(this).load_state(\'%s\');' % href)
        page.li.close()
    page.ul.close()
    page.div.close()  # btn-group
    return page


def about_this_page(cmdline=True, config=None):
    """Write a blurb about the generation of this page, including the
    command-line arguments used, and the content of the configuration
    files.

    Parameters
    ----------
    cmdline : `bool`, optional, default: `True`
        print the full command-line as recorded by the system.
    config : `str`, `list`, optional
        the absolute path(s) to one or a number of INI files used in this
        process.

    Returns
    -------
    page : :class:`~gwpy.html.markup.page`
        the HTML page to be inserted into the #main <div>.
    """
    page = markup.page()
    page.div(class_='row')
    page.div(class_='col-md-12')
    if cmdline:
        page.h2("On the command-line")
        page.p("This page was generated using %s with the following "
               "command-line call:"
               % (markup.oneliner.a('GWSumm',
                                    href='https://github.com/gwpy/gwsumm')))
        if sys.argv[0].endswith('__main__.py'):
            package = sys.argv[0].rsplit(os.path.sep, 2)[1]
            cmdline = 'python -m %s %s' % (package, ' '.join(sys.argv[1:]))
        elif (os.path.split(sys.argv[0])[0] in
              os.environ.get('PATH', '').split(os.pathsep)):
            cmdline = ' '.join([os.path.basename(sys.argv[0])] + sys.argv[1:])
        else:
            cmdline = ' '.join(sys.argv)
        page.pre(cmdline)
    if config is not None:
        page.h2("Configuration files")
        page.p("GWSumm processes are configured through INI-format files. The "
               "following files were passed on the comand-line and are "
               "reproduced in full below.")
        if isinstance(config, str):
            config = config.split(',')
        page.div(class_='panel-group', id="accordion")
        for i, cpfile in enumerate(config):
            page.div(class_='panel panel-default')
            page.a(href='#file%d' % i, **{'data-toggle': 'collapse',
                                          'data-parent': '#accordion'})
            page.div(class_='panel-heading')
            page.h4(os.path.basename(cpfile), class_='panel-title')
            page.div.close()
            page.a.close()
            page.div(id_='file%d' % i, class_='panel-collapse collapse')
            page.div(class_='panel-body')
            page.pre(highlight_syntax(cpfile, 'ini'))
            page.div.close()
            page.div.close()
            page.div.close()
        page.div.close()

    page.div.close()
    page.div.close()
    return page


def base_map_dropdown(this, class_='btn-group pull-left base-map', id_=None,
                      bases=dict()):
    """Construct a dropdown menu that links to a version of the current
    page on another server, based on a new base.
    """
    # format id
    if id_:
        id_ = dict(id_=id_)
    else:
        id_ = dict()
    # format links
    baselinks = [markup.oneliner.a(key, title=key, **{'data-new-base': val})
                 for (key, val) in bases.iteritems() if key != this]
    # slam it all together
    page = markup.page()
    if baselinks:
        page.div(class_=class_, **id_)
        page.add(str(dropdown(this, baselinks,
                              class_='navbar-brand dropdown-toggle')))
        page.div.close()
    else:
        page.div(str(this), class_='navbar-brand')
    return page
