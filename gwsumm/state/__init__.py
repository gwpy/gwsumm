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

import re

from gwpy.segments import (Segment, SegmentList, DataQualityFlag)

from .. import globalv
from ..config import (GWSummConfigParser, NoSectionError)
from ..utils import re_cchar
from ..segments import get_segments

from .registry import (get_state, get_states, register_state)

from .core import SummaryState
from .registry import (get_state, get_states, register_state)
from .all import (ALLSTATE, generate_all_state)

__all__ = ['ALLSTATE', 'SummaryState', 'get_state', 'get_states',
           'register_state', 'load_states', 'generate_all_state']


def load_states(config, section='states'):
    """Read and format a list of `SummaryState` definitions from the
    given :class:`~configparser.ConfigParser`
    """
    config = GWSummConfigParser.from_configparser(config)
    # get [start, stop) job interval
    start = config.getint(section, 'gps-start-time')
    end = config.getint(section, 'gps-end-time')
    # parse the [states] section into individual state definitions
    try:
        states = dict(config.nditems(section))
    except NoSectionError:
        states = {}
    for state in states:
        if not (config.has_section('state-%s' % state) or
                config.has_section('state %s' % state)):
            section = 'state-%s' % state
            config.add_section(section)
            config.set(section, 'name', state)
            config.set(section, 'definition', states[state])

    # parse each state section into a new state
    states = []
    for section in config.sections():
        if re.match('state[-\s]', section):
            states.append(register_state(
                SummaryState.from_ini(config, section)))

    # register All state
    try:
        all_ = generate_all_state(start, end)
    except ValueError:
        all_ = get_state(ALLSTATE)
        pass
    else:
        states.insert(0, all_)

    return states
