# -*- coding: utf-8 -*-
# Copyright (C) Alex Urban (2019)
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

"""Unit tests for gwsumm.html.bootstrap
"""

__author__ = 'Alex Urban <alexander.urban@ligo.org>'

import pytest
from datetime import datetime

from gwdetchar.utils import parse_html

from .. import bootstrap

# global variables
DATE = datetime.strptime('20140410', '%Y%m%d')
BACKWARD = '<a class="nav-link step-back" title="Step backward">&laquo;</a>'
CAL = ('<a id="calendar" class="nav-link dropdown-toggle" title="Show/hide '
       'calendar" data-date="10-04-2014" data-date-format="dd-mm-yyyy" '
       'data-viewmode="{}">{}</a>')
FORWARD = '<a class="nav-link step-forward" title="Step forward">&raquo;</a>'


# test utilities

def test_banner():
    banner = bootstrap.banner('Test', subtitle='Subtest')
    assert parse_html(str(banner)) == parse_html(
        '<div class="banner">\n<h1>Test</h1>\n<p>Subtest</p>\n</div>')
    # test with classes
    banner_wclass = bootstrap.banner(
        'Test', subtitle='Subtest', titleclass='test', subtitleclass='subtest')
    assert parse_html(str(banner_wclass)) == parse_html(
        '<div class="banner">\n<h1 class=\"test\">Test</h1>\n'
        '<p class=\"subtest\">Subtest</p>\n</div>')


@pytest.mark.parametrize('mode, datefmt', [
    ('day', 'April 10 2014'),
    ('week', 'Week of April 10 2014'),
    ('month', 'April 2014'),
    ('year', '2014'),
])
def test_calendar(mode, datefmt):
    backward, cal, forward = bootstrap.calendar(DATE, mode=mode)
    assert parse_html(str(backward)) == parse_html(str(BACKWARD))
    assert parse_html(str(cal)) == parse_html(
        CAL.format('%ss' % mode, datefmt))
    assert parse_html(str(forward)) == parse_html(str(FORWARD))


def test_calendar_no_mode():
    # test with no Mode
    with pytest.raises(ValueError) as exc:
        bootstrap.calendar(DATE)
    assert str(exc.value).startswith('Cannot generate calendar for Mode')


def test_wrap_content():
    content = bootstrap.wrap_content('test')
    assert parse_html(str(content)) == parse_html(
        '<div class="container-fluid" id="main">\ntest\n</div>')


def test_state_switcher():
    switcher = bootstrap.state_switcher([('Test', '#test')])
    assert parse_html(str(switcher)) == parse_html(
        '<ul class="nav navbar-nav">\n<li class="nav-item dropdown">\n'
        '<a class="nav-link dropdown-toggle" href="#" id="states" role='
        '"button" title="Show/hide state menu" data-toggle="dropdown">Test</a>'
        '\n<div class="dropdown-menu dropdown-menu-right state-switch shadow" '
        'id="statemenu">\n<h6 class="dropdown-header">Select below to view '
        'this page in another state (different time segments).</h6>\n<div '
        'class="dropdown-divider"></div>\n<a class="dropdown-item state" '
        'title="Test" id="state_test" onclick="jQuery(this).load_state'
        '(&quot;#test&quot;);">Test</a>\n</div>\n</li>\n</ul>')


def test_base_map_dropdown():
    menu = bootstrap.base_map_dropdown('test', id_='id')
    assert parse_html(str(menu)) == parse_html(
        '<div class="navbar-brand border border-white '
        'rounded" id="id">test</div>')
    # test with bases
    menu_wbases = bootstrap.base_map_dropdown('test', bases={'key': 'value'})
    assert parse_html(str(menu_wbases)) == parse_html(
        '<div class="dropdown base-map">\n<a href="#" class="navbar-brand '
        'nav-link border border-white rounded dropdown-toggle" role="button" '
        'data-toggle="dropdown">test</a>\n<div class="dropdown-menu '
        'dropdown-1-col shadow">\n<a title="key" class="dropdown-item" '
        'data-new-base="value">key</a>\n</div>\n</div>')
