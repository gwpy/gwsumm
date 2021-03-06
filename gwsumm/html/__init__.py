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

"""HTML helpers

HTML output is built upon the
`markup.py module <http://markup.sourceforge.net/>`_ and primarily
formatted to fit the
`twitter bootstrap library <http://getbootstrap.com/>`_.
"""

from .static import (
    get_css,
    get_js,
)
from .html5 import (
    _expand_path,
    load_state,
    load,
    comments_box,
    ldvw_qscan,
    dialog_box,
    overlay_canvas,
)
from .bootstrap import (
    banner,
    calendar,
    wrap_content,
    state_switcher,
    base_map_dropdown,
)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
