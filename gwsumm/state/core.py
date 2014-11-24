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

"""Definition of the `SummaryState` class.
"""

import datetime
import re

from astropy.time import Time

from gwpy.detector import get_timezone_offset
from gwpy.segments import (Segment, SegmentList, DataQualityFlag)
from gwpy.time import (to_gps, from_gps)

from .. import globalv
from ..config import (GWSummConfigParser, NoOptionError, DEFAULTSECT)
from ..utils import re_cchar
from ..segments import get_segments


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
                 description=None, definition=None, hours=None, key=None,
                 url=None):
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
        if definition:
            self.definition = re.sub('(\s|\n)', '', definition)
        else:
            self.definition = None
        self.key = key
        self.hours = hours
        self.url = url
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
        if d:
            self._definition = str(d)
        else:
            self._definition = None

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
        from .all import ALLSTATE, generate_all_state
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
        # get hours
        hours = params.pop('hours', None)
        if hours is not None:
            segs = re.split('(,|, )', hours)[::2]
            hours = []
            offset = 0
            for seg in segs:
                try:
                    # parse hour segment
                    hours.append(map(float, seg.split('-', 1)))
                except ValueError:
                    # parse time-zone
                    if seg == segs[-1]:
                        if seg.lower() == 'utc':
                            offset = 0
                        elif seg.lower() == 'local':
                            try:
                                ifo = config.get(DEFAULTSECT, 'ifo')
                            except NoOptionError:
                                raise ValueError("The relevant IFO must be "
                                                 "given either from the --ifo "
                                                 "command-line option, or the "
                                                 "[DEFAULT] section of any "
                                                 "INI file")
                            offset = get_timezone_offset(ifo, from_gps(start))
                        else:
                            offset = get_timezone_offset(seg, from_gps(start))
                    else:
                        raise
            # apply time-zone
            for i, (h0, h1) in enumerate(hours):
                hours[i] = (h0 - offset / 3600., h1 - offset / 3600.)
        # generate state
        if name == ALLSTATE:
            return generate_all_state(start, end, register=False, **params)
        else:
            return cls(name, valid=[(start, end)], hours=hours, **params)

    def fetch(self, config=GWSummConfigParser(), **kwargs):
        """Finalise this state by fetching its defining segments,
        either from global memory, or from the segment database
        """
        # check we haven't done this before
        if self.ready:
            return self
        # otherwise find out which flags we need
        if self.definition:
            segs = get_segments(self.definition, self.valid, config=config,
                                **kwargs)
            self.valid = segs.valid
            self.active = segs.active
        else:
            start = config.getfloat(DEFAULTSECT, 'gps-start-time')
            end = config.getfloat(DEFAULTSECT, 'gps-end-time')
            self.valid = [(start, end)]
            self.active = self.valid
        # restrict to given hours
        if self.hours:
            segs_ = SegmentList()
            # get start day
            d = Time(float(self.start), format='gps', scale='utc').datetime
            d.replace(hour=0, minute=0, second=0, microsecond=0)
            end_ = Time(float(self.end), format='gps', scale='utc').datetime
            while d < end_:
                # get GPS of day
                t = to_gps(d)
                # for each [start, end) hour pair, build a segment
                for h0, h1 in self.hours:
                    segs_.append(Segment(t + h0 * 3600, t + h1*3600))
                # increment and return
                d += datetime.timedelta(1)
            self.valid &= segs_
            self.active &= segs_
        # FIXME
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
