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

"""Test suite for the mode module

"""

import datetime

from common import unittest
from gwsumm import (globalv, mode)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


class ModeTests(unittest.TestCase):
    """`TestCase` for the channels.py module
    """
    @classmethod
    def setUpClass(cls):
        cls._defaultmode = mode.get_mode()

    @classmethod
    def tearDownClass(cls):
        mode.set_mode(cls._defaultmode)

    def test_get_mode(self):
        m = mode.get_mode()
        self.assertEqual(m, globalv.MODE)

    def test_set_mode(self):
        mode.set_mode(0)
        self.assertEqual(globalv.MODE, mode.Mode(0))
        self.assertEqual(globalv.MODE, mode.Mode.static)
        mode.set_mode('GPS')
        self.assertEqual(globalv.MODE, mode.Mode.gps)

    def test_get_base(self):
        date = datetime.date(2015, 9, 22)
        for m, basestr in zip(
                ['DAY', 'WEEK', 'MONTH', 'YEAR'],
                ['day/20150922', 'week/20150922',
                 'month/201509', 'year/2015']):
            mode.set_mode(m)
            base = mode.get_base(date)
            self.assertEqual(base, basestr)
