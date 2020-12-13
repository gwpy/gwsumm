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

"""Unit tests for gwsumm.html.static
"""

__author__ = 'Alex Urban <alexander.urban@ligo.org>'

from collections import OrderedDict

from .. import static


# test simple utils

def test_get_css():
    css = static.get_css()
    assert isinstance(css, OrderedDict)
    # test dict keys
    assert list(css.keys()) == [
        'font-awesome',
        'font-awesome-solid',
        'gwbootstrap',
    ]
    # test list of files
    css_files = list(x.split('/')[-1] for x in css.values())
    assert css_files == [
        'fontawesome.min.css',
        'solid.min.css',
        'gwbootstrap.min.css',
    ]


def test_get_js():
    js = static.get_js()
    assert isinstance(js, OrderedDict)
    # test dict keys
    assert list(js.keys()) == [
        'jquery',
        'jquery-ui',
        'moment',
        'bootstrap',
        'fancybox',
        'datepicker',
        'gwbootstrap',
    ]
    # test list of files
    js_files = list(x.split('/')[-1] for x in js.values())
    assert js_files == [
        'jquery-3.5.1.min.js',
        'jquery-ui.min.js',
        'moment.min.js',
        'bootstrap.bundle.min.js',
        'jquery.fancybox.min.js',
        'bootstrap-datepicker.min.js',
        'gwbootstrap-extra.min.js',
    ]
