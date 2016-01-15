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

"""A `SummaryState` defines a sub-set of time over which a `~gwsumm.tabs.Tab`
should be processed.
Each `SummaryState` is normally tied to one or more data-quality flags marking
times during which each of the LIGO instruments was operating in a certain
configuration, or was subject to a known noise interference.

==================
The state registry
==================

GWSumm defines a state 'registry', simply a record of all `SummaryState`
objects that have been defined (and registered) so far in a given program.
The registry just makes remembering states in complicated programs a little
easier.

Any `SummaryState` can be registered with an arbitrary name as follows::

    >>> from gwsumm.state.registry import register_state
    >>> register_state(mystate, 'my state')

and can be recovered later::

    >>> from gwsumm.state.registry import get_state
    >>> mystate = get_state('my state')

=============
API reference
=============

.. autosummary::
   :toctree: api

   SummaryState
   get_state
   get_states
   register_state

"""

from .core import SummaryState
from .registry import (get_state, get_states, register_state)
from .all import (ALLSTATE, generate_all_state)

__all__ = ['ALLSTATE', 'SummaryState', 'get_state', 'get_states',
           'register_state', 'generate_all_state']
