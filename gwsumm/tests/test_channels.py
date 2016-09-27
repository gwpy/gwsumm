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

"""Test suite

"""

from common import (unittest, empty_globalv_CHANNELS)
from gwsumm import (globalv, channels)
from gwsumm.mode import (get_mode, set_mode)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

TEST_NAME = 'H1:GDS-CALIB_STRAIN'
TREND_NAME = 'H1:TEST-TREND_CHANNEL.rms,m-trend'
TREND_NAME2 = 'H1:TEST-TREND_CHANNEL.mean,m-trend'
TREND_NAME3 = 'H1:TEST-TREND_CHANNEL_3.mean'
TREND_NAME4 = 'H1:TEST-TREND_CHANNEL_4.mean'


class ChannelTests(unittest.TestCase):
    """`TestCase` for the channels.py module
    """
    @classmethod
    def setUpClass(cls):
        cls._defaultmode = get_mode()

    @classmethod
    def tearDownClass(cls):
        set_mode(cls._defaultmode)

    @empty_globalv_CHANNELS
    def test_get_channel(self):
        """Test get_channel for normal channels
        """
        nchan = len(globalv.CHANNELS)
        chan = channels.get_channel(TEST_NAME)
        self.assertEqual(len(globalv.CHANNELS), nchan+1)
        self.assertEqual(chan.name, TEST_NAME)
        chan2 = channels.get_channel(TEST_NAME)
        self.assertEqual(len(globalv.CHANNELS), nchan+1)
        self.assertIs(chan, chan2)

    @empty_globalv_CHANNELS
    def test_get_channel_trend(self):
        """Test get_channel for trends

        `get_channel` should query for the trend and the underlying
        raw channel
        """
        # test simple query
        nchan = len(globalv.CHANNELS)
        chan = channels.get_channel(TREND_NAME)
        self.assertEqual(len(globalv.CHANNELS), nchan+2)
        # test that raw doesn't get built again
        chan = channels.get_channel(TREND_NAME2)
        self.assertEqual(len(globalv.CHANNELS), nchan+3)
        # test that raw matches trend
        raw = channels.get_channel(TREND_NAME.split('.')[0])
        self.assertEqual(len(globalv.CHANNELS), nchan+3)
        self.assertEqual(raw.name, TREND_NAME.split('.')[0])
        # test default trend type
        set_mode('DAY')
        chan = channels.get_channel(TREND_NAME3)
        self.assertEqual(chan.type, 'm-trend')
        set_mode('GPS')
        chan = channels.get_channel(TREND_NAME4)
        self.assertEqual(chan.type, 's-trend')

    @empty_globalv_CHANNELS
    def test_get_channels(self):
        names = [TEST_NAME, TREND_NAME, TREND_NAME2]
        nchan = len(globalv.CHANNELS)
        chans = channels.get_channels(names)
        self.assertEqual(len(globalv.CHANNELS), nchan+3)
        for name, chan in zip(names, chans):
            self.assertEqual(name, chan.ndsname)
