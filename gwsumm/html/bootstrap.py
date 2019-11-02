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

from MarkupPy import markup

from gwdetchar.io import html

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


def calendar(date, tag='a', class_='nav-link dropdown-toggle',
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
    calendar : `list`
        a list of three oneliner strings of HTML containing the calendar
        text and a triggering dropdown
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
    # get navigation objects
    backward = markup.oneliner.a(
        '&laquo;', class_='nav-link step-back', title='Step backward')
    cal = markup.oneliner.a(
        datestring, id_=id_, class_=class_, title='Show/hide calendar',
        **{'data-date': data_date, 'data-date-format': 'dd-mm-yyyy',
           'data-viewmode': '%ss' % mode.name})
    forward = markup.oneliner.a(
        '&raquo;', class_='nav-link step-forward', title='Step forward')
    return [backward, cal, forward]


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
    out.div(class_='container-fluid', id_='main')
    out.add(str(page))
    out.div.close()
    return out


def state_switcher(states, default=0):
    """Build a state switch button, including all of the given
    states, with the default selected by index
    """
    current, chref = states[default]
    page = markup.page()
    page.ul(class_='nav navbar-nav')
    page.li(class_='nav-item dropdown')
    page.a(str(current), class_='nav-link dropdown-toggle', href='#',
           id_='states', role='button', title='Show/hide state menu',
           **{'data-toggle': 'dropdown'})
    page.div(
        class_='dropdown-menu dropdown-menu-right state-switch shadow',
        id_='statemenu',
    )
    page.h6('Select below to view this page in another state (different '
            'time segments).', class_='dropdown-header')
    page.div('', class_='dropdown-divider')
    for i, (state, href) in enumerate(states):
        page.a(str(state), class_='dropdown-item state', title=str(state),
               id_='state_%s' % re_cchar.sub('_', str(state)).lower(),
               onclick='jQuery(this).load_state(\'%s\');' % href)
    page.div.close()  # dropdown-menu dropdown-menu-right
    page.li.close()  # nav-item dropdown state-switch
    page.ul.close()  # nav navbar-nav
    return page


def base_map_dropdown(this, id_=None, bases=dict()):
    """Construct a dropdown menu that links to a version of the current
    page on another server, based on a new base.
    """
    # format id
    if id_:
        id_ = dict(id_=id_)
    else:
        id_ = dict()
    # format links
    baselinks = [markup.oneliner.a(
        key, title=key, class_='dropdown-item', **{'data-new-base': val}
    ) for (key, val) in bases.items() if key != this]
    # slam it all together
    page = markup.page()
    if baselinks:
        page.div(class_='dropdown base-map', **id_)
        page.add(str(html.dropdown(
            this,
            baselinks,
            class_='navbar-brand nav-link border border-white '
                   'rounded dropdown-toggle',
        )))
        page.div.close()
    else:
        page.div(
            str(this),
            class_='navbar-brand border border-white rounded',
            **id_,
        )
    return page
