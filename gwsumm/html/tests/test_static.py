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

import os.path

from collections import OrderedDict

from gwdetchar.io.html import (CSS_FILES as GWDETCHAR_CSS_FILES,
                               JS_FILES as GWDETCHAR_JS_FILES)

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
    assert list(css_files.values()) == GWDETCHAR_CSS_FILES


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
        'gwbootstrap-extra',
    ]
    # test list of files
    js_files = list(js.values())
    assert len(js_files) == len(GWDETCHAR_JS_FILES) + 3
    assert not set(js_files).isdisjoint(GWDETCHAR_JS_FILES)
    assert os.path.basename(js_files[1]) == 'jquery-ui.min.js'
    assert os.path.basename(js_files[2]) == 'moment.min.js'
    assert os.path.basename(js_files[5]) == 'bootstrap-datepicker.min.js'
    assert os.path.basename(js_files[6]) == 'gwbootstrap-extra.min.js'
