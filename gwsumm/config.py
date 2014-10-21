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

__all__ = _cp__all__ + ['InterpolationMissingOptionError', 'GWSummConfigParser']


class GWSummConfigParser(ConfigParser):
    optionxform = str

    def nditems(self, section, **kwargs):
        try:
            items = ConfigParser.items(self, section, **kwargs)
        except ValueError as e:
            msg = str(e)
            raise ValueError('[%s]: %s' % (section, msg))
        return [i for i in items if i[0] not in self._defaults]

    def read(self, filenames):
        readok = ConfigParser.read(self, filenames)
        if isinstance(filenames, (unicode, str)):
            filenames = filenames.split(',')
        for fp in filenames:
            if fp not in readok:
                raise IOError("Cannot read file: %s" % fp)
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
