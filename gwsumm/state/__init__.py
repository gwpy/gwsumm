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

"""Segment-state definitions for GWSumm
"""

import re
try:
    import configparser
except ImportError:
    import ConfigParser as configparser

from gwpy.segments import (Segment, SegmentList, DataQualityFlag)

from .. import globalv
from ..utils import re_cchar
from ..segments import get_segments
from .registry import (get_state, register_state)

ALLSTATE = 'All'


class SummaryState(DataQualityFlag):
    """Tab run state - defining an interferometer operating mode and
    associated segments over which to generate summary information
    """
    def __init__(self, name, valid=SegmentList(), active=SegmentList(),
                 description=None, definition=None):
        """Initialise a new `SummaryState`
        """
        # allow users to specify valid as (start, end)
        if (isinstance(valid, Segment) or
                (isinstance(valid, tuple) and len(valid) == 2 and
                 not isinstance(valid[0], (tuple)))):
            valid = [valid]
        super(SummaryState, self).__init__(name=name, valid=valid,
                                           active=active)
        self.description = description
        self.definition = definition
        self.ready = False

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, n):
        self._name = n

    @property
    def tag(self):
        return re_cchar.sub("_", self.name).upper()

    @property
    def start(self):
        return self.extent[0]

    @property
    def end(self):
        return self.extent[1]

    @classmethod
    def from_ini(cls, section, config=configparser.ConfigParser()):
        """Create a new `SummaryState` object from its definition in
        a :class:`configparser.ConfigParser` object.
        """
        # get span times
        start = config.getint('general', 'gps-start-time')
        end = config.getint('general', 'gps-end-time')
        # get parameters
        try:
            params = dict(config.nditems(section))
        except AttributeError:
            params = dict(config.items(section))
        # parse name
        name = params.pop('name', section)
        if re.match('state[-\s]', name):
            name = section[6:]
        return cls(name, valid=[(start, end)], **params)

    def fetch(self, config=configparser.ConfigParser(), **kwargs):
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
        new.include = self.include
        new.exclude = self.exclude
        new.ready = self.ready
        return new
    copy.__doc__ = DataQualityFlag.copy.__doc__

    def __str__(self):
        return self.name


def load_states(config, section='states'):
    """Read and format a list of `SummaryState` definitions from the
    given :class:`~configparser.ConfigParser`
    """
    # get [start, stop) job interval
    start = config.getint('general', 'gps-start-time')
    end = config.getint('general', 'gps-end-time')

    # parse the [states] section into individual state definitions
    try:
        states = dict(config.nditems(section))
    except configparser.NoSectionError:
        states = {}
    for state in states:
        if not (config.has_section('state-%s' % state) or
                config.has_section('state %s' % state)):
            section = 'state-%s' % state
            config.add_section(section)
            config.set(section, 'name', state)
            config.set(section, 'definition', states[state])

    # parse each state section into a new state
    for section in config.sections():
        if re.match('state[-\s]', section):
            register_state(SummaryState.from_ini(section, config))

    # register All state
    try:
        register_state(SummaryState(ALLSTATE, valid=[(start, end)]))
    except ValueError:
        pass

    return globalv.STATES
