# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013-2016)
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

"""This module defines the `Tab` API, and all of the built-in tab objects
"""

# core
from .registry import *
from .core import *
from .builtin import *
from .misc import *

# data
from .data import *

# application-specific extras
from .sei import *
from .guardian import *
from .stamp import *
from .management import *
from .etg import *
from .fscan import *
from .gracedb import *

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
