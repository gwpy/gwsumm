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

from gwpy.segments import DataQualityFlag

from .. import globalv
from ..utils import re_cchar
from ..segments import get_segments

ALLSTATE = 'All'


class SummaryState(DataQualityFlag):
    """Tab run state - defining an interferometer operating mode and
    associated segments over which to generate summary information
    """
    def __init__(self, name, start, stop, description=None):
        """Initialise a new `SummaryState`
        """
        super(SummaryState, self).__init__()
        self.name = name
        self.valid = [(start, stop)]
        self.active = [(start, stop)]
        self.description = description
        self.goodflags = []
        self.badflags = []
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

    def fetch(self, config=configparser.ConfigParser()):
        """Finalise this state by fetching its defining segments,
        either from global memory, or from the segment database
        """
        # check we haven't done this before
        if self.ready:
            return self
        # if state == 'ALL': return full span as 'active'
        if self.name == ALLSTATE:
            self.active = self.valid
            self.ready = True
            return self
        # otherwise find out which flags we need
        goodflags, badflags = self.parse_description()
        for flag in badflags:
            self -= get_segments(flag, self.valid, config=config)
        for flag in goodflags:
            self &= get_segments(flag, self.valid, config=config)
        self.ready = True
        return self

    def parse_description(self):
        """Parse the configuration for this state.

        Returns
        -------
        flags : `tuple`
            a 2-tuple containing lists of flags defining required ON
            and OFF segments respectively for this state
        """
        # find flags
        div = re.compile("[&!]")
        divs = div.findall(self.description)
        keys = div.split(self.description)
        # load flags and vetoes
        flags = []
        vetoes = []
        for i,key in enumerate(keys):
            # get veto bit
            if key.startswith("-") or (i != 0 and divs[i-1] == "!"):
                vetoes.append(key)
            else:
                flags.append(key)
        return flags, vetoes

    def copy(self):
        new = super(SummaryState, self).copy()
        new.description = self.description
        new.goodflags = self.goodflags
        new.badflags = self.badflags
        new.ready = self.ready
        return new
    copy.__doc__ = DataQualityFlag.copy.__doc__


def load_states(config, section='states'):
    """Read and format a list of `SummaryState` definitions from the
    given :class:`~configparser.ConfigParser`
    """
    # get [start, stop) job interval
    start = config.getint('general', 'gps-start-time')
    end = config.getint('general', 'gps-end-time')

    try:
        items = list(config.nditems('states'))
    except configparser.NoSectionError:
        items = []
    if ALLSTATE not in map(str.lower, zip(*items)[0]):
        items.insert(0, (ALLSTATE, None))

    for name, def_ in items:
        if not def_:
            try:
                def_ = config.get('state-%s' % name.lower(), 'description')
            except (configparser.NoSectionError, configparser.NoOptionError):
                def_ = None
        globalv.STATES[name] = SummaryState(name, start, end,
                                            description=def_)
