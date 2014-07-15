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

ALLSTATE = 'All'

__all__ = ['SummaryState', 'get_state', 'get_states', 'register_state',
           'load_states']


class SummaryState(DataQualityFlag):
    """An operating state over which to process a `~gwsumm.tabs.DataTab`.

    Parameters
    ----------
    name : `str`
        name for this state
    valid : `~gwpy.segments.SegmentList`, optional
        list of valid segments
    active : `~gwpy.segments.SegmentList`, optional
        list of active segments
    description : `str`, optional
        text describing what this state means
    definition : `str`, optional
        logical combination of flags that define valid and active segments
        for this state (see :attr:`documentation <SummaryState.definition>`
        for details)
    key : `str`, optional
        registry key for this state, defaults to :attr:`~SummaryState.name`
    """
    def __init__(self, name, valid=SegmentList(), active=SegmentList(),
                 description=None, definition=None, key=None):
        """Initialise a new `SummaryState`
        """
        # allow users to specify valid as (start, end)
        if (isinstance(valid, Segment) or
                (isinstance(valid, tuple) and len(valid) == 2 and
                 not isinstance(valid[0], tuple))):
            valid = [valid]
        super(SummaryState, self).__init__(name=name, valid=valid,
                                           active=active)
        self.description = description
        self.definition = definition
        self.key = key
        if valid and active:
            self.ready = True
        else:
            self.ready = False

    @property
    def name(self):
        """Name of this state
        """
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    @property
    def tag(self):
        """File tag for images generated using this state
        """
        return re_cchar.sub("_", self.name).upper()

    @property
    def start(self):
        """GPS start time of this state's validity.
        """
        return self.extent[0]

    @property
    def end(self):
        """GPS end time of this state's validity
        """
        return self.extent[1]

    @property
    def definition(self):
        """The combination of data-quality flags that define this `SummaryState`

        For example::

           >>> state = SummaryState(definition='L1:DMT-SCIENCE:1')

        would define a `SummaryState` based on the validity and activity of
        the single flag ``'L1:DMT-SCIENCE:1'`` (the science-mode flag for the
        LIGO Livingston Observatory interferometer).
        Similarly::

           >>> state = SummaryState(
                   definition='L1:DMT-SCIENCE:1&!L1:DMT-LIGHTDIP_10_PERCENT:1')

        would define a `SummaryState` as active when the ``'L1:DMT-SCIENCE:1'``
        flag was active and the ``'L1:DMT-LIGHTDIP_10_PERCENT:1'`` flag was
        not active.

        The following logical identifiers are acceptable:

        ==  ===================================================================
        &   Union (i.e. flag1 **and** flag2 must be active)
        \|  Intersection (i.e. flag1 **or** flag2 must be active)
        &!  One-sided difference (i.e. flag1 is active and flag2 **is not**\
            active)
        !=  Two-sided difference (i.e. flag1 is active and flag2 is not **OR**\
            flag2 is active and flag2 is not)
        ==  ===================================================================

        :type: str
        """
        return self._definition

    @definition.setter
    def definition(self, d):
        self._definition = str(d)

    @property
    def key(self):
        """The registry key for this `SummaryState`.

        :type: `str`
        """
        if self._key is not None:
            return self._key
        else:
            return self._name

    @key.setter
    def key(self, k):
        self._key = k

    @classmethod
    def from_ini(cls, config, section):
        """Create a new `SummaryState` from a section in a `ConfigParser`.

        Parameters
        ----------
        config : :class:`~gwsumm.config.GWConfigParser`
            customised configuration parser containing given section
        section : `str`
            name of section to parse

        Returns
        -------
        `SummaryState`
            a new state, with attributes set from the options in the
            configuration
        """
        config = GWSummConfigParser.from_configparser(config)
        # get span times
        start = config.getint(section, 'gps-start-time')
        end = min(globalv.NOW, config.getint(section, 'gps-end-time'))
        # get parameters
        params = dict(config.nditems(section))
        # parse name
        name = params.pop('name', section)
        if re.match('state[-\s]', name):
            name = section[6:]
        return cls(name, valid=[(start, end)], **params)

    def fetch(self, config=GWSummConfigParser(), **kwargs):
        """Finalise this state by fetching its defining segments,
        either from global memory, or from the segment database
        """
        # check we haven't done this before
        if self.ready:
            return self
        # if state == 'ALL': return full span as 'active'
        if self.name == ALLSTATE:
            now = min(self.valid[0][1], globalv.NOW)
            self.active = [(self.valid[0][0], now)]
            self.ready = True
            return self
        # otherwise find out which flags we need
        segs = get_segments(self.definition, self.valid, config=config,
                            **kwargs)
        self.valid = segs.valid
        self.active = segs.active
        self.ready = True
        return self

    def copy(self):
        new = super(SummaryState, self).copy()
        new.description = self.description
        new.definition = self.definition
        new.ready = self.ready
        return new
    copy.__doc__ = DataQualityFlag.copy.__doc__

    def __str__(self):
        return self.name


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
        all_ = register_state(
                   SummaryState(ALLSTATE,
                                valid=SegmentList([Segment(start, end)])))
    except ValueError:
        pass
    else:
        states.insert(0, all_)

    return states
