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

"""A `Tab` is a single, configurable page of output, containing some data.
Each `Tab` is written in its own HTML page, and can be written to contain
any set of data, with any format.

Abstract Tab classes
--------------------

The :mod:`gwsumm.tabs` module provides a few base classes from which users
can write their own tab `classes <class>`.
These include:

.. autosummary::
   :nosignatures:
   :toctree: api

   Tab
   SummaryArchiveMixin
   DataTabBase


Available Tabs
--------------

Along with the base classes, GWSumm provides a set of `Tab` subclasses covering
some common data displays.
These include:

.. autosummary::
   :nosignatures:
   :toctree: api

   ExternalTab
   PlotTab
   StateTab

Each of these tab classes is also available with the `SummaryArchiveMixin`
applied, allowing users to archive their output with a calendar, and crosslinks:

.. autosummary::
   :nosignatures:
   :toctree: api

   ArchivedExternalTab
   ArchivedPlotTab
   ArchivedStateTab


Data generation
---------------

All of the above tabs simply format data that has already been generated.
GWSumm also provides a `DataTab`, allowing users to configure plots and have
them generated on-the-fly by reading in the raw data, processing, and saving
plots to disk.

This is achieved via the following class

.. autosummary::
   :toctree: api

   DataTab
   EventTriggerTab

Process-specific tabs
---------------------

The following `Tab` classes are provided to interface to a specific GW
detector characterization or analysis group or process:

.. autosummary::
   :nosignatures:
   :toctree: api

   AccountingTab
   GuardianTab
   HvetoTab
   SEIWatchDogTab
   StampPEMTab

"""

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

from .registry import *
from .core import *
from .builtin import *
from .data import *
from .ihope import *
from .hveto import *
from .sei import *
from .guardian import *
from .stamp import *
from .management import *
from .etg import *
from .fscan import *
from .gracedb import *
