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

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


STATICDIR = os.path.join(os.path.dirname(__file__), 'static')

# set bootstrap info
BOOTSTRAP_VERSION = "3.3.6"
BOOTSTRAP_CDN = '//netdna.bootstrapcdn.com/bootstrap/%s' % BOOTSTRAP_VERSION
BOOTSTRAP_CSS = '%s/css/bootstrap.min.css' % BOOTSTRAP_CDN
BOOTSTRAP_JS = '%s/js/bootstrap.min.js' % BOOTSTRAP_CDN

# build list of javascript resources
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

# build list of CSS resources
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
