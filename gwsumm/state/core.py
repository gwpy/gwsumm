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
import operator
from configparser import (NoOptionError, DEFAULTSECT)

from astropy.time import Time

from gwpy.detector import get_timezone_offset
from gwpy.segments import (Segment, SegmentList, DataQualityFlag)
from gwpy.time import (to_gps, from_gps)

from .. import globalv
from ..config import (GWSummConfigParser)
from ..utils import re_cchar
from ..segments import get_segments
from ..data import get_timeseries

MATHOPS = {
    '<': operator.lt,
    '<=': operator.le,
    '=': operator.eq,
    '>=': operator.ge,
    '>': operator.gt,
    '==': operator.is_,
    '!=': operator.is_not,
}


class SummaryState(DataQualityFlag):
    """An operating state over which to process a `~gwsumm.tabs.DataTab`.

    Parameters
    ----------
    name : `str`
        name for this state
    known : `~gwpy.segments.SegmentList`, optional
        list of known segments
    active : `~gwpy.segments.SegmentList`, optional
        list of active segments
    description : `str`, optional
        text describing what this state means
    definition : `str`, optional
        logical combination of flags that define known and active segments
        for this state (see :attr:`documentation <SummaryState.definition>`
        for details)
    key : `str`, optional
        registry key for this state, defaults to :attr:`~SummaryState.name`
    """
    MATH_DEFINITION = re.compile(r'(%s)' % '|'.join(MATHOPS.keys()))

    def __init__(self, name, known=SegmentList(), active=SegmentList(),
                 description=None, definition=None, hours=None, key=None,
                 filename=None, url=None):
        """Initialise a new `SummaryState`
        """
        # allow users to specify known as (start, end)
        if (isinstance(known, Segment) or
                (isinstance(known, tuple) and len(known) == 2 and
                 not isinstance(known[0], tuple))):
            known = [known]
        super(SummaryState, self).__init__(name=name, known=known,
                                           active=active)
        self.description = description
        if definition:
            self.definition = re.sub(r'(\s|\n)', '', definition)
        else:
            self.definition = None
        self.key = key
        self.hours = hours
        self.url = url
        if known and active:
            self.ready = True
        else:
            self.ready = False
        self.filename = filename

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
        r"""The combination of data-quality flags that define this `SummaryState`

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
        config = GWSummConfigParser.from_configparser(config)
        # get span times
        start = config.getint(section, 'gps-start-time')
        end = min(globalv.NOW, config.getint(section, 'gps-end-time'))
        # get parameters
        params = dict(config.nditems(section))
        # parse name
        name = params.pop('name', section)
        if re.match(r'state[-\s]', name):
            name = section[6:]
        # get hours
        hours = params.pop('hours', None)
        if hours is not None:
            segs = re.split(r'(,|, )', hours)[::2]
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
        return cls(name, known=[(start, end)], hours=hours, **params)

    def _fetch_segments(self, config=GWSummConfigParser(), **kwargs):
        kwargs.setdefault('url', self.url)
        segs = get_segments([self.definition], self.known, config=config,
                            **kwargs)[self.definition].round(contract=True)
        self.known = segs.known
        self.active = segs.active
        return self

    def _fetch_data(self, channel, thresh, op, config=GWSummConfigParser(),
                    **kwargs):
        # if 0 is out-of-state, allowing padding gaps with 0.
        if (
                (op == '<' and thresh <= 0.) or
                (op == '<=' and thresh < 0.) or
                (op == '>' and thresh >= 0.) or
                (op == '>=' and thresh > 0.) or
                (op in ('=', '==') and thresh != 0.) or
                (op == '!=' and thresh == 0.)
        ):
            kwargs.setdefault('pad', 0.)
        data = get_timeseries(channel, self.known, config=config, **kwargs)
        for ts in data:
            if isinstance(thresh, (float, int)) and ts.unit is not None:
                thresh *= ts.unit
            segs = MATHOPS[op](ts, thresh).to_dqflag()
            try:
                globalv.SEGMENTS[self.definition] += segs
            except KeyError:
                globalv.SEGMENTS[self.definition] = segs
        return self._fetch_segments(query=False)

    def _read_segments(self, filename):
        segs = DataQualityFlag.read(filename, self.definition)
        # XXX HACK around malformed segment files with no segment_summary table
        if segs.active and not segs.known:
            segs.known = type(segs.active)(segs.active)
        if self.known:
            self.known = self.known & segs.known
            self.active = self.known & segs.active
        else:
            self.known = segs.known
            self.active = segs.active
        return self

    def fetch(self, config=GWSummConfigParser(),
              segmentcache=None, segdb_error='raise',
              datacache=None, datafind_error='raise', nproc=1, nds=None,
              **kwargs):
        """Finalise this state by fetching its defining segments,
        either from global memory, or from the segment database
        """
        # check we haven't done this before
        if self.ready:
            return self
        # fetch data
        match = self.MATH_DEFINITION.search(str(self.definition))
        if self.filename:
            self._read_segments(self.filename)
        elif self.definition and match is not None:
            channel, thresh = self.definition.split(match.groups()[0])
            channel = channel.rstrip()
            thresh = float(thresh.strip())
            self._fetch_data(channel, thresh, match.groups()[0], config=config,
                             cache=datacache, nproc=nproc, nds=nds,
                             datafind_error=datafind_error, **kwargs)
        # fetch segments
        elif self.definition:
            self._fetch_segments(config=config, cache=segmentcache,
                                 segdb_error=segdb_error, **kwargs)
        # fetch null
        else:
            start = config.getfloat(DEFAULTSECT, 'gps-start-time')
            end = config.getfloat(DEFAULTSECT, 'gps-end-time')
            self.known = [(start, end)]
            self.active = self.known
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
            self.known &= segs_
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
