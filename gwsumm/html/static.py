# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2016)
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

"""HTML <head> helphers

This module mainly declares the resources used by standard on HTML pages
"""

import os.path
from collections import OrderedDict
from gwdetchar.io.html import (CSS_FILES, JS_FILES)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__credits__ = 'Alex Urban <alexander.urban@ligo.org>'


STATICDIR = os.path.join(os.path.dirname(__file__), 'static')

# build list of javascript resources
JS = OrderedDict((
    ('jquery', JS_FILES[0]),
    ('moment', JS_FILES[1]),
    ('bootstrap', JS_FILES[2]),
    ('fancybox', JS_FILES[3]),
    ('datepicker', (
        '//cdnjs.cloudflare.com/ajax/libs/bootstrap-datepicker/'
        '1.6.0/js/bootstrap-datepicker.min.js')),
    ('bootstrap-ligo', JS_FILES[4]),
    ('gwdetchar', JS_FILES[5]),
    ('gwsumm', os.path.join(STATICDIR, 'gwsumm.min.js')),
))

# build list of CSS resources
CSS = OrderedDict((
    ('bootstrap', CSS_FILES[0]),
    ('fancybox', CSS_FILES[1]),
    ('datepicker', (
        '//cdnjs.cloudflare.com/ajax/libs/bootstrap-datepicker/'
        '1.6.0/css/bootstrap-datepicker.min.css')),
    ('bootstrap-ligo', CSS_FILES[2]),
    ('gwdetchar', CSS_FILES[3]),
    ('gwsumm', os.path.join(STATICDIR, 'gwsumm.min.css')),
))


# -- utilities ----------------------------------------------------------------

def get_css():
    """Return a `dict` of CSS files to link in the HTML <head>
    """
    return CSS


def get_js():
    """Return a `dict` of javascript files to link in the HTML <head>
    """
    return JS
