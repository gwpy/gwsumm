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

"""Helper functions for twitter-bootstrap HTML constructs
"""

import sys
import getpass
import datetime
import os.path

from ... import version
__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

from .. import markup
from ..utils import highlight_syntax


def banner(title, subtitle=None):
    """Construct a banner heading in bootstrap format

    Parameters
    ----------
    """
    page = markup.page()
    page.div(class_='container', id='header')
    page.div(class_='row')
    page.div(class_='col-sm-12')
    # title
    if subtitle is not None:
        page.h1('%s : ' % title, class_="inline-block")
    else:
        page.h1('%s' % title)
    # subtitle
    if subtitle is not None:
        page.h1(subtitle, class_="normal inline-block")
    page.div.close()
    page.div.close()
    page.div.close()

    return page


def navbar(links, brand=None, states=None):
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
    page = markup.page()
    page.div(id_='nav-wrapper')
    page.nav(id_='nav', class_='navbar', role='navigation')
    page.div(class_='container')
    page.div(class_='row')
    page.div(class_='col-md-12')
    if states:
        page.add(str(states))
    page.div(class_='navbar-header')
    if brand is not None:
        page.add(brand)
    page.button(type='button', class_='navbar-toggle',
                **{'data-toggle': 'collapse',
                   'data-target': '.navbar-collapse'})
    page.span('', class_='icon-bar')
    page.span('', class_='icon-bar')
    page.span('', class_='icon-bar')
    page.button.close()
    page.div.close()
    page.div(class_='navbar-collapse collapse pull-left')
    page.ul(class_='nav navbar-nav')
    for i, link in enumerate(links):
        if isinstance(link, (list, tuple)) and isinstance(link[1], basestring):
            page.li()
            text, link = link
            page.a(text, href=link)
        elif (isinstance(link, (list, tuple)) and
              isinstance(link[1], (list, tuple))):
            page.li(class_="dropdown")
            page.add(str(dropdown(*link)))
        else:
            page.li()
            page.add(str(link))
        page.li.close()
    page.ul.close()
    page.div.close()
    page.div.close()
    page.div.close()
    page.div.close()
    page.nav.close()
    page.div.close()

    return page


def dropdown(text, links, active=None):
    """Construct a dropdown menu in bootstrap format

    Parameters
    ----------
    text : `str`
        dropdown menu header
    links : `list`
        list of (Link text, linkurl) tuples

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
    page.a(href='#', class_='dropdown-toggle', **{'data-toggle': 'dropdown'})
    page.add(text)
    page.b('', class_='caret')
    page.a.close()
    # dropdown elements
    page.ul(class_='dropdown-menu')
    for i, link in enumerate(links):
        if link is None:
            page.li(class_='divider')
        elif i == active:
            page.li(class_='active')
        else:
            page.li()
        if isinstance(link, (tuple, list)):
            page.a(link[0], href=link[1])
        elif link is not None:
            page.add(str(link))
        page.li.close()
    page.ul.close()
    return page


def calendar(date, tag='a', class_='dropdown-toggle', id_='calendar',
             dateformat='%B %d %Y', mode='days'):
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
    datestring = date.strftime(dateformat)
    data_date = date.strftime('%d-%m-%Y')
    return markup.element(tag)(datestring, id_=id_, class_=class_,
                               **{'data-date': data_date,
                                  'data-date-format': 'dd-mm-yyyy',
                                  'data-viewmode': mode})


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


def footer(user=True, about=None, issues=True, span=12):
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
    page.FOOTER()
    page.div(class_='container')
    page.div(class_='row')
    page.div(class_='col-md-%d' % span)
    if user:
        page.p('This page was generated by user %s at %s.'
               % (getpass.getuser(),
                  datetime.datetime.now().strftime('%H:%m on %B %d %Y')))
    if issues:
        page.a('Report an issue', href='https://github.com/gwpy/gwsumm/issues',
               target='_blank')
    if about is not None:
        page.a('How was this page generated?', href=about)
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
    page.div(class_="btn-group pull-right")
    page.a(class_='btn dropdown-toggle', href='#', id_='state',
           **{'data-toggle': 'dropdown'})
    page.add(current.name)
    page.b('', class_='caret')
    page.a.close()
    page.ul(class_='dropdown-menu', id_='statemenu')
    for i, (state, href) in enumerate(states):
        page.li()
        page.a(state.name, href=href)
        page.li.close()
    page.ul.close()
    page.div.close()
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
