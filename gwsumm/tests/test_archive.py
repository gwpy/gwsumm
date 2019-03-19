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

import pytest

import h5py

from numpy import (random, testing as nptest)

from gwpy.table import EventTable
from gwpy.timeseries import (TimeSeries, StateVector)
from gwpy.spectrogram import Spectrogram
from gwpy.segments import (Segment, SegmentList)

from gwsumm import (archive, data, globalv, channels, triggers)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

TEST_DATA = TimeSeries([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], epoch=100,
                       unit='meter', sample_rate=1, channel='X1:TEST-CHANNEL',
                       name='TEST DATA')
TEST_DATA.channel = channels.get_channel(TEST_DATA.channel)


# -- utilities ----------------------------------------------------------------

def empty_globalv():
    globalv.DATA = type(globalv.DATA)()
    globalv.SPECTROGRAMS = type(globalv.SPECTROGRAMS)()
    globalv.SEGMENTS = type(globalv.SEGMENTS)()
    globalv.TRIGGERS = type(globalv.TRIGGERS)()


def create(data, **metadata):
    SeriesClass = metadata.pop('series_class', TimeSeries)
    d = SeriesClass(data, **metadata)
    d.channel = channels.get_channel(d.channel)
    if not d.name:
        d.name = d.channel.texname
    return d


# -- tests --------------------------------------------------------------------

def test_write_archive(delete=True):
    empty_globalv()
    data.add_timeseries(TEST_DATA)
    data.add_timeseries(create([1, 2, 3, 4, 5],
                               dt=60., channel='X1:TEST-TREND.mean'))
    data.add_timeseries(create([1, 2, 3, 2, 1],
                               series_class=StateVector,
                               channel='X1:TEST-STATE_VECTOR'))
    data.add_spectrogram(create([[1, 2, 3], [3, 2, 1], [1, 2, 3]],
                                series_class=Spectrogram,
                                channel='X1:TEST-SPECTROGRAM'))
    t = EventTable(random.random((100, 5)), names=['time', 'a', 'b', 'c', 'd'])
    t.meta['segments'] = SegmentList([Segment(0, 100)])
    triggers.add_triggers(t, 'X1:TEST-TABLE,testing')
    fname = tempfile.mktemp(suffix='.h5', prefix='gwsumm-tests-')
    try:
        archive.write_data_archive(fname)
        archive.write_data_archive(fname)  # test again to validate backups
    finally:
        if delete and os.path.isfile(fname):
            os.remove(fname)
    return fname


def test_read_archive():
    fname = test_write_archive(delete=False)
    empty_globalv()
    try:
        archive.read_data_archive(fname)
    finally:
        os.remove(fname)
    # check timeseries
    ts = data.get_timeseries('X1:TEST-CHANNEL',
                             [(100, 110)], query=False).join()
    nptest.assert_array_equal(ts.value, TEST_DATA.value)
    for attr in ['epoch', 'unit', 'sample_rate', 'channel', 'name']:
        assert getattr(ts, attr) == getattr(TEST_DATA, attr)
    # check trend series
    ts = data.get_timeseries('X1:TEST-TREND.mean,m-trend', [(0, 300)],
                             query=False).join()
    assert ts.channel.type == 'm-trend'
    assert ts.span == (0, 300)
    # check triggers
    t = triggers.get_triggers('X1:TEST-TABLE', 'testing', [(0, 100)])
    assert len(t) == 100


def test_archive_load_table():
    t = EventTable(random.random((100, 5)),
                   names=['a', 'b', 'c', 'd', 'e'])
    empty = EventTable(names=['a', 'b'])
    try:
        fname = tempfile.mktemp(suffix='.h5', prefix='gwsumm-tests-')
        h5file = h5py.File(fname)
        # check table gets archived and read transparently
        archive.archive_table(t, 'test-table', h5file)
        t2 = archive.load_table(h5file['test-table'])
        nptest.assert_array_equal(t.as_array(), t2.as_array())
        assert t.dtype == t2.dtype
        # check empty table does not get archived, with warning
        with pytest.warns(UserWarning):
            n = archive.archive_table(empty, 'test-empty', h5file)
        assert n is None
        assert 'test-empty' not in h5file
    finally:
        if os.path.exists(fname):
            os.remove(fname)
