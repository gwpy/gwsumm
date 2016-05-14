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

"""HTML output utilities

HTML output is built upon the
`markup.py module <http://markup.sourceforge.net/>`_ and primarily
formatted to fit the
`twitter bootstrap library <http://getbootstrap.com/>`_.
"""

import os.path

try:
    from collections import OrderedDict
except:
    from ordereddict import OrderedDict

from .html5 import *
from .bootstrap import *
from . import markup

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

# ----------------------------------------------------------------------------
# Collect CSS and JS scripts

STATICDIR = os.path.join(os.path.dirname(__file__), 'static')

JS = OrderedDict()
JS['jquery'] = '//code.jquery.com/jquery-1.12.3.min.js'
JS['moment'] = (
    '//cdnjs.cloudflare.com/ajax/libs/moment.js/2.13.0/moment.min.js')
JS['bootstrap'] = BOOTSTRAP_JS
JS['fancybox'] = (
    '//cdnjs.cloudflare.com/ajax/libs/fancybox/2.1.5/jquery.fancybox.pack.js')
JS['datepicker'] = (
    '//cdnjs.cloudflare.com/ajax/libs/bootstrap-datepicker/'
    '1.6.0/js/bootstrap-datepicker.min.js')
JS['bootstrap-ligo'] = os.path.join(STATICDIR, 'bootstrap-ligo.min.js')
JS['gwsumm'] = os.path.join(STATICDIR, 'gwsumm.min.js')

CSS = OrderedDict()
CSS['bootstrap'] = BOOTSTRAP_CSS
CSS['fancybox'] = (
    '//cdnjs.cloudflare.com/ajax/libs/fancybox/2.1.5/jquery.fancybox.css')
CSS['datepicker'] = (
    '//cdnjs.cloudflare.com/ajax/libs/bootstrap-datepicker/'
    '1.6.0/css/bootstrap-datepicker.min.css')
CSS['bootstrap-ligo'] = os.path.join(STATICDIR, 'bootstrap-ligo.min.css')
CSS['gwsumm'] = os.path.join(STATICDIR, 'gwsumm.min.css')


def get_css():
    """Return a `dict` of CSS files to link in the HTML <head>
    """
    return CSS


def get_js():
    """Return a `dict` of javascript files to link in the HTML <head>
    """
    return JS


# ----------------------------------------------------------------------------
# Write HTML table

def data_table(headers, data, **class_):
    """Write a <table> with a single row of headers, followed by multiple
    rows of data
    """
    # construct arguments
    page = markup.page()
    if 'table' in class_:
        page.table(class_=class_['table'])
    else:
        page.table(class_='table table-condensed table-hover')

    # get row class
    trargs = {}
    if 'tr' in class_:
        trargs['class_'] = class_['tr']

    # write headers
    page.thead()
    page.tr(**trargs)
    thargs = {}
    if 'th' in class_:
        thargs['class_'] = class_['th']
    for th in headers:
        page.th(th, **thargs)
    page.tr.close()
    page.thead.close()

    # write data
    tdargs = {}
    if 'td' in class_:
        tdargs['class_'] = class_['td']
    for row in data:
        page.tr(**trargs)
        for td in row:
            page.td(td, **tdargs)
        page.tr.close()
    page.table.close()
    return page


