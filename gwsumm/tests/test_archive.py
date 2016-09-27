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

"""Tests for `gwsumm.archive`

"""

import os
import tempfile
from functools import wraps

from numpy import testing as nptest

from gwpy.timeseries import TimeSeries

from common import unittest
from gwsumm import (archive, data, globalv, channels)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

TEST_DATA = TimeSeries([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], epoch=100,
                       unit='meter', sample_rate=1, channel='X1:TEST-CHANNEL',
                       name='TEST DATA')
TEST_DATA.channel = channels.get_channel(TEST_DATA.channel)

def empty_globalv_DATA(f):
    @wraps(f)
    def wrapped_f(*args, **kwargs):
        _data = globalv.DATA
        globalv.DATA = type(globalv.DATA)()
        try:
            return f(*args, **kwargs)
        finally:
            globalv.DATA = _data
    return wrapped_f


class ArchiveTests(unittest.TestCase):
    """`TestCase` for the `gwsumm.archive` module
    """
    def test_write_archive(self, delete=True):
        data.add_timeseries(TEST_DATA)
        fname = tempfile.mktemp(suffix='.hdf', prefix='gwsumm-tests-')
        try:
            archive.write_data_archive(fname)
        finally:
            if delete and os.path.isfile(fname):
                os.remove(fname)
        return fname

    @empty_globalv_DATA
    def test_read_archive(self):
        fname = self.test_write_archive(delete=False)
        try:
            archive.read_data_archive(fname)
        finally:
            os.remove(fname)
        ts = data.get_timeseries('X1:TEST-CHANNEL',
                                 [(100, 110)], query=False).join()
        nptest.assert_array_equal(ts.value, TEST_DATA.value)
        for attr in ['epoch', 'unit', 'sample_rate', 'channel', 'name']:
            self.assertEqual(getattr(ts, attr), getattr(TEST_DATA, attr))
