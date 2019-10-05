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
from .registry import (
    register_tab,
    get_tab,
)
from .core import (
    _MetaTab,
    BaseTab,
    StaticTab,
    GpsTab,
    IntervalTab,
    EventTab,
    Tab,
    ParentTab,
    TabList,
)
from .builtin import (
    ExternalTab,
    PlotTab,
    StateTab,
    UrlTab,
)
from .misc import (
    AboutTab,
    Error404Tab,
)

# data
from .data import (ProcessedTab, DataTab)

# application-specific extras
from .sei import SEIWatchDogTab
from .guardian import GuardianTab
from .stamp import StampPEMTab
from .management import AccountingTab
from .etg import EventTriggerTab
from .fscan import FscanTab
from .gracedb import GraceDbTab

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
