# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2017)
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

"""Input/output utilities
"""

import os.path
import re

from astropy.io.registry import IORegistryError

from gwpy.frequencyseries import FrequencySeries

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

HDF5_FILENAME = re.compile(r'(?P<ext>(.hdf5|.hdf|.h5))\/')


def read_frequencyseries(filename):
    """Read a `~gwpy.frequencyseries.FrequencySeries` from a file

    IF using HDF5, the filename can be given as a combined filename/path, i.e.
    ``test.hdf5/path/to/dataset``.

    Parameters
    ----------
    filename : `str`
        path of file to read

    Returns
    -------
    series : `~gwpy.frequencyseries.FrequencySeries`
        the data as read

    Raises
    ------
    astropy.io.registry.IORegistryError
        if the input format cannot be identified or is not registered
    """
    # try and parse path in HDF5 file if given
    try:
        ext = HDF5_FILENAME.search(filename).groupdict()['ext']
    except AttributeError:  # no match
        kwargs = {}
    else:
        kwargs = {'path': filename.rsplit(ext, 1)[1]}
    # read file
    try:
        return FrequencySeries.read(filename, **kwargs)
    except IORegistryError:
        if filename.endswith('.gz'):
            fmt = os.path.splitext(filename[:-3])[-1]
        else:
            fmt = os.path.splitext(filename)[-1]
        return FrequencySeries.read(filename, format=fmt.lstrip('.'), **kwargs)
