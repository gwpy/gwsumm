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

import pytest

from astropy import units

from gwpy.detector import ChannelList

from gwsumm import (globalv, channels)
from gwsumm.mode import (get_mode, set_mode)

from gwpy.detector import Channel

from .common import empty_globalv_CHANNELS

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

TEST_NAME = 'H1:GDS-CALIB_STRAIN'
TREND_NAME = 'H1:TEST-TREND_CHANNEL.rms,m-trend'
TREND_NAME2 = 'H1:TEST-TREND_CHANNEL.mean,m-trend'
TREND_NAME3 = 'H1:TEST-TREND_CHANNEL_3.mean'
TREND_NAME4 = 'H1:TEST-TREND_CHANNEL_4.mean'

DEFAULT_MODE = get_mode()


def teardown_module():
    """Undo any set_mode() operations from this module
    """
    set_mode(DEFAULT_MODE)


@empty_globalv_CHANNELS
def test_get_channel():
    """Test :func:`gwsumm.channels.get_channel`
    """
    nchan = len(globalv.CHANNELS)

    # test simple query
    chan = channels.get_channel(TEST_NAME)
    assert len(globalv.CHANNELS) == nchan + 1
    assert chan.name == TEST_NAME

    # make sure that querying again returns the same object
    chan2 = channels.get_channel(TEST_NAME)
    assert len(globalv.CHANNELS) == nchan + 1
    assert chan2 is chan


@empty_globalv_CHANNELS
def test_get_channel_trend():
    """Test get_channel for trends

    `get_channel` should query for the trend and the underlying
    raw channel
    """
    # test simple query
    nchan = len(globalv.CHANNELS)
    chan = channels.get_channel(TREND_NAME)
    assert len(globalv.CHANNELS) == nchan + 2

    # test that raw doesn't get built again
    chan = channels.get_channel(TREND_NAME2)
    assert len(globalv.CHANNELS) == nchan + 3

    # test that raw matches trend
    raw = channels.get_channel(TREND_NAME.split('.')[0])
    assert len(globalv.CHANNELS) == nchan + 3
    assert raw.name == TREND_NAME.split('.')[0]

    # test default trend type
    chan = channels.get_channel(TREND_NAME3)
    assert chan.type is None


@empty_globalv_CHANNELS
def test_get_channels():
    names = [TEST_NAME, TREND_NAME, TREND_NAME2]
    nchan = len(globalv.CHANNELS)
    chans = channels.get_channels(names)
    # trend channels auto create entry for the upstream channel so '+ 1'
    assert len(globalv.CHANNELS) == nchan + 3 + 1
    for name, chan in zip(names, chans):
        assert name == chan.ndsname


@empty_globalv_CHANNELS
def test_update_missing_channel_params():
    # define empty channel
    chan = channels.get_channel('X1:TEST:1')
    assert chan.unit is None

    # update using kwargs
    channels.update_missing_channel_params('X1:TEST:1', unit='meter')
    assert chan.unit == units.meter
    chan.unit = None

    # update from another channel
    c2 = Channel('X1:TEST:1', unit='V')
    channels.update_missing_channel_params(c2)
    assert chan.unit == units.volt


@pytest.mark.parametrize('cstr, clist', [
    ('X1:TEST,Y1:TEST,\nZ1:TEST', ['X1:TEST', 'Y1:TEST', 'Z1:TEST']),
    ('X1:TEST.rms,m-trend,Y1:TEST.mean,s-trend',
     ['X1:TEST.rms,m-trend', 'Y1:TEST.mean,s-trend']),
])
def test_split(cstr, clist):
    assert channels.split(cstr) == clist


@pytest.mark.parametrize('cstr, clist', [
    ('X1:TEST + Y1:TEST', ['X1:TEST', 'Y1:TEST']),
    ('X1:TEST.mean,m-trend * Y1:TEST + 1',
     ['X1:TEST.mean,m-trend', 'Y1:TEST']),
    ('G1:DER_DATA_H-rms - 4 + G1:DER_DATA_BLAH',
     ['G1:DER_DATA_H-rms', 'G1:DER_DATA_BLAH']),
])
def test_split_combination(cstr, clist):
    split = channels.split_combination(cstr)
    assert isinstance(split, ChannelList)
    assert list(map(str, split)), clist
