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

from .. import version
from .html5 import *
from .bootstrap import *
from . import markup

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

# ----------------------------------------------------------------------------
# Collect CSS and JS scripts

sharedir = os.path.join(os.path.split(__file__)[0], os.path.pardir,
                        os.path.pardir, 'share', 'gwsumm', 'html')

JQUERYJS = ['//code.jquery.com/jquery-1.11.0.min.js',
            '//cdnjs.cloudflare.com/ajax/libs/moment.js/2.4.0/moment.min.js']

BOOTSTRAPCSS = CSS
BOOTSTRAPJS = JS

FANCYBOXCSS = ['//cdnjs.cloudflare.com/ajax/libs/fancybox/2.1.5/'
               'jquery.fancybox.css']
FANCYBOXJS = ['//cdnjs.cloudflare.com/ajax/libs/fancybox/2.1.5/'
              'jquery.fancybox.pack.js']

GWSUMMCSS = list(map(lambda p: os.path.normpath(os.path.join(sharedir, p)),
                     ['datepicker.css', 'gwsummary.css']))
GWSUMMJS = list(map(lambda p: os.path.normpath(os.path.join(sharedir, p)),
                    ['datepicker.js', 'gwsummary.js']))

CSS = BOOTSTRAPCSS + FANCYBOXCSS + GWSUMMCSS
JS = JQUERYJS + BOOTSTRAPJS + FANCYBOXJS + GWSUMMJS


def ifo_css(ifo):
    return os.path.normpath(os.path.join(
        sharedir, 'gwsummary_%s.css' % str(ifo).lower()))

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
        page.table()

    # get row class
    trargs = {}
    if 'tr' in class_:
        trargs['class_'] = class_['tr']

    # write headers
    page.tr(**trargs)
    thargs = {}
    if 'th' in class_:
        thargs['class_'] = class_['th']
    for th in headers:
        page.th(th, **thargs)
    page.tr.close()

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


