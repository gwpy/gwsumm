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

from six import string_types

from MarkupPy import markup

from gwdetchar.io import html as gwhtml

from ..mode import (Mode, get_mode)
from ..utils import re_cchar

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

# set <meta> for bootstrap
META = {'viewport': 'width=device-width, initial-scale=1.0'}

# -----------------------------------------------------------------------------
# variable HTML constructs

def banner(title, subtitle=None, titleclass=None, subtitleclass=None):
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
    banner : `~MarkupPy.markup.page`
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
            raise ValueError("Cannot generate calendar for Mode %s" % mode)
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
    page : :class:`~MarkupPy.markup.page`, `str`
        HTML content to be wrapped
    span : `int`
        column span of content, default: 'full' (``12``)

    Returns
    -------
    wrappedpage : :class:`~MarkupPy.markup.page`
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
                 for (key, val) in bases.items() if key != this]
    # slam it all together
    page = markup.page()
    if baselinks:
        page.div(class_=class_, **id_)
        page.add(str(gwhtml.dropdown(this, baselinks,
                                     class_='navbar-brand dropdown-toggle')))
        page.div.close()
    else:
        page.div(str(this), class_='navbar-brand', **id_)
    return page
