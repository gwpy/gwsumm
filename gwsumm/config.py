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

"""Thin wrapper around configparser
"""

import os.path
import re
import sys
from StringIO import StringIO

if sys.version_info[0] >= 3:
    from configparser import *
    from configparser import InterpolationMissingOptionError
    from configparser import __all__ as _cp__all__
else:
    from ConfigParser import *
    from ConfigParser import InterpolationMissingOptionError
    from ConfigParser import __all__ as _cp__all__

try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

from gwpy.time import tconvert

__all__ = _cp__all__ + ['InterpolationMissingOptionError', 'GWSummConfigParser']


class GWSummConfigParser(ConfigParser):
    # preserve case in options
    optionxform = str
    # disable colon separator
    OPTCRE = re.compile(
        r'(?P<option>[^=\s][^=]*)\s*(?P<vi>[=])\s*(?P<value>.*)$')
    # use OrderedDict for dict type
    _dict = OrderedDict
    # set interpolation match
    _interpvar_re = re.compile(r"%\(([^)]+)\)s")

    def ndoptions(self, section, **kwargs):
        try:
            options = ConfigParser.options(self, section, **kwargs)
        except ValueError as e:
            e.args = ('[%s]: %s' % (section, str(e)),)
            raise
        return [o for o in options if o not in self._defaults]

    def nditems(self, section, **kwargs):
        try:
            items = ConfigParser.items(self, section, **kwargs)
        except ValueError as e:
            e.args = ('[%s]: %s' % (section, str(e)),)
            raise
        return [i for i in items if i[0] not in self._defaults]

    def read(self, filenames):
        readok = ConfigParser.read(self, filenames)
        if isinstance(filenames, (unicode, str)):
            filenames = filenames.split(',')
        for fp in filenames:
            if fp not in readok:
                raise IOError("Cannot read file: %s" % fp)
        self.files = map(os.path.abspath, readok)
        return readok

    @classmethod
    def from_configparser(cls, cp):
        """Copy an existing :class:`~ConfigParser.ConfigParser`.
        """
        # set up temporary buffer
        buf = StringIO()
        # write to buffer
        cp.write(buf)
        buf.seek(0)
        # read new GWSummConfigParser
        new = cls()
        new.readfp(buf)
        return new

    def __repr__(self):
        return '<GWSummConfigParser()>'

    def interpolate_section_names(self, **kwargs):
        """Interpolate a specific key in a section name using val
        """
        for section in self.sections():
            s = section
            for key in self._interpvar_re.findall(section):
                try:
                    val = kwargs[key]
                except KeyError:
                    raise InterpolationMissingOptionError(
                        '[%s]' % section, section, key, '%%(%s)s' % key)
                s = section % {key: val}
            self._sections[s] = self._sections.pop(section)

    def set_date_options(self, start, end, section=DEFAULTSECT):
        """Set datetime options in [DEFAULT] based on the given times

        The following options are set

        - `gps-start-time` - the integer GPS start time of this job
        - `gps-end-time` - the integer GPS end time of this job
        - `yyyy` - the four-digit year of the start date
        - `mm` - the two-digit month of the start date
        - `dd` - the two-digit day-of-month of the start date
        - `yyyymm` - the six-digit year and month of the start date
        - `yyyymmdd` - the eight-digit year-month-day of the start date
        - `duration` - the duration of the job (seconds)

        Additionally, if LAL is available, the following extra options
        are also set

        - `leap-seconds` - the number of leap seconds for the start date
        - `gps-start-time-noleap` - the leap-corrected integer GPS start time
        - `gps-end-time-noleap` - the leap-corrected integer GPS end time

        """
        utc = tconvert(start)
        self.set(section, 'gps-start-time', str(int(start)))
        self.set(section, 'gps-end-time', str(int(end)))
        self.set(section, 'yyyy', utc.strftime('%Y'))
        self.set(section, 'yy', utc.strftime('%y'))
        self.set(section, 'mm', utc.strftime('%m'))
        self.set(section, 'dd', utc.strftime('%d'))
        self.set(section, 'yyyymm', utc.strftime('%Y%m'))
        self.set(section, 'yyyymmdd', utc.strftime('%Y%m%d'))
        self.set(section, 'duration', str(int(end - start)))
        try:
            from lal import GPSLeapSeconds as leap_seconds
        except ImportError:
            pass
        else:
            nleap = leap_seconds(int(start))
            self.set(section, 'leap-seconds', str(nleap))
            self.set(section, 'gps-start-time-noleap',
                       str(int(start) - nleap))
            self.set(section, 'gps-end-time-noleap', str(int(end) - nleap))
