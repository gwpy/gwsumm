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
CALENDAR = """<a class="navbar-brand step-back" title="Step back" onclick="stepDate(-1)">&laquo;</a>
<a id="calendar" class="navbar-brand dropdown-toggle" title="Show/hide calendar" data-date="10-04-2014" data-date-format="dd-mm-yyyy" data-viewmode="{}">
{}
<b class="caret"></b>
</a>
<a class="navbar-brand step-forward" title="Step forwards" onclick="stepDate(1)">&raquo;</a>"""  # nopep8


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
    cal = bootstrap.calendar(DATE, mode=mode)
    assert parse_html(str(cal)) == parse_html(CALENDAR.format(
        '%ss' % mode, datefmt))


def test_calendar_no_mode():
    # test with no Mode
    with pytest.raises(ValueError) as exc:
        cal = bootstrap.calendar(DATE)
    assert str(exc.value).startswith('Cannot generate calendar for Mode')


def test_wrap_content():
    content = bootstrap.wrap_content('test')
    assert parse_html(str(content)) == parse_html(
        '<div class="container" id="main">\ntest\n</div>')


def test_state_switcher():
    switcher = bootstrap.state_switcher([('Test', '#test')])
    assert parse_html(str(switcher)) == parse_html(
        '<div class="btn-group pull-right state-switch">\n'
        '<a class="navbar-brand dropdown-toggle" href="#" id="states" '
        'title="Show/hide state menu" data-toggle="dropdown">\n'
        'Test\n<b class="caret"></b>\n</a>\n'
        '<ul class="dropdown-menu" id="statemenu">\n'
        '<li class="dropdown-header">Select an option below to view '
        'these data in another state (different time segments).</li>\n'
        '<li class="divider"></li>\n<li>\n'
        '<a class="state" title="Test" id="state_test" '
        'onclick="$(this).load_state(&quot;#test&quot;);">Test</a>\n'
        '</li>\n</ul>\n</div>')


def test_base_map_dropdown():
    menu = bootstrap.base_map_dropdown('test', id_='id')
    assert parse_html(str(menu)) == parse_html(
        '<div class=\"navbar-brand\" id=\"id\">test</div>')
    # test with bases
    menu_wbases = bootstrap.base_map_dropdown('test', bases={'key': 'value'})
    assert parse_html(str(menu_wbases)) == parse_html(
        '<div class="btn-group pull-left base-map">\n'
        '<a href="#" class="navbar-brand dropdown-toggle" '
        'data-toggle="dropdown">\ntest\n<b class="caret"></b>\n</a>\n'
        '<ul class="dropdown-menu">\n<li>\n'
        '<a title="key" data-new-base="value">key</a>\n</li>\n</ul>\n</div>')
