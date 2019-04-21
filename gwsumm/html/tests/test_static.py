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

from gwdetchar.io.html import (CSS_FILES, JS_FILES)

from .. import static


# test simple utils

def test_get_css():
    css = static.get_css()
    assert isinstance(css, OrderedDict)
    # test dict keys
    keys = list(css.keys())
    assert keys == [
        'bootstrap',
        'fancybox',
        'datepicker',
        'bootstrap-ligo',
        'gwdetchar',
        'gwsumm',
    ]
    # test list of files
    files = list(css.values())
    assert len(files) == len(CSS_FILES) + 2
    assert set(CSS_FILES) <= set(files)
    assert files[2].endswith('bootstrap-datepicker.min.css')
    assert files[5].endswith('gwsumm.min.css')


def test_get_js():
    js = static.get_js()
    assert isinstance(js, OrderedDict)
    # test dict keys
    keys = list(js.keys())
    assert keys == [
        'jquery',
        'moment',
        'bootstrap',
        'fancybox',
        'datepicker',
        'bootstrap-ligo',
        'gwdetchar',
        'gwsumm',
    ]
    # test list of files
    files = list(js.values())
    assert len(files) == len(JS_FILES) + 2
    assert set(JS_FILES) <= set(files)
    assert files[4].endswith('bootstrap-datepicker.min.js')
    assert files[7].endswith('gwsumm.min.js')
