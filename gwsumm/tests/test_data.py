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

"""Tests for `gwsumm.data`

"""

import os.path
import operator
import tempfile
import shutil

try:  # py3
    from urllib.request import urlopen
except ImportError:  # py2
    from urllib2 import urlopen

from glue.lal import (Cache, CacheEntry)

from gwpy.timeseries import TimeSeries
from gwpy.detector import Channel
from gwpy.segments import (Segment, SegmentList)

from common import (unittest, empty_globalv_CHANNELS)
from gwsumm import (data, globalv)
from gwsumm.data import (utils, mathutils)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

LOSC_DATA = {
    'H1:LOSC-STRAIN': ['https://losc.ligo.org/s/events/GW150914/'
                       'H-H1_LOSC_4_V1-1126259446-32.gwf'],
    'L1:LOSC-STRAIN': ['https://losc.ligo.org/s/events/GW150914/'
                       'L-L1_LOSC_4_V1-1126259446-32.gwf'],
}
LOSC_SEGMENTS = SegmentList([Segment(1126259446, 1126259478)])


def download(remote, target=None):
    """Download a file
    """
    if target is None:
        suffix = os.path.splitext(remote)[1]
        _, target = tempfile.mkstemp(suffix=suffix, prefix='gwsumm-tests-')
    response = urlopen(remote)
    with open(target, 'w') as f:
        f.write(response.read())
    return target


class DataTests(unittest.TestCase):
    """`TestCase` for the `gwsumm.data` module
    """
    @classmethod
    def setUpClass(cls):
        cls.FRAMES = {}
        cls._tempdir = tempfile.mkdtemp(prefix='gwsumm-test-data-')
        # get data
        for channel in LOSC_DATA:
            cls.FRAMES[channel] = Cache()
            for gwf in LOSC_DATA[channel]:
                target = os.path.join(cls._tempdir, os.path.basename(gwf))
                download(gwf, target)
                cls.FRAMES[channel].append(CacheEntry.from_T050017(target))

    @classmethod
    def tearDownClass(cls):
        # remove the temporary data
        shutil.rmtree(cls._tempdir)

    # -- test utilities -------------------------

    def test_find_frame_type(self):
        channel = Channel('L1:TEST-CHANNEL')
        self.assertEqual(data.find_frame_type(channel), 'L1_R')
        channel = Channel('C1:TEST-CHANNEL')
        self.assertEqual(data.find_frame_type(channel), 'R')
        channel = Channel('H1:TEST-CHANNEL.rms,s-trend')
        self.assertEqual(data.find_frame_type(channel), 'H1_T')
        channel = Channel('H1:TEST-CHANNEL.rms,m-trend')
        self.assertEqual(data.find_frame_type(channel), 'H1_M')
        channel = Channel('H1:TEST-CHANNEL.rms,reduced')
        self.assertEqual(data.find_frame_type(channel), 'H1_LDAS_C02_L2')
        channel = Channel('H1:TEST-CHANNEL.rms,online')
        self.assertEqual(data.find_frame_type(channel), 'H1_lldetchar')

    def test_get_channel_type(self):
        self.assertEqual(data.get_channel_type('L1:TEST-CHANNEL'), 'adc')
        self.assertEqual(data.get_channel_type('G1:DER_DATA_HL'), 'proc')
        self.assertEqual(data.get_channel_type('H1:GDS-CALIB_STRAIN'), 'proc')
        self.assertEqual(data.get_channel_type('V1:GDS-CALIB_STRAIN'), 'adc')

    @empty_globalv_CHANNELS
    def test_make_globalv_key(self):
        fftparams = utils.get_fftparams('L1:TEST-CHANNEL',
            stride=123.456, window='test-window', method='test-method')
        key = utils.make_globalv_key('L1:TEST-CHANNEL', fftparams)
        self.assertEqual(
            key, 'L1:TEST-CHANNEL;test-method;;;test-window;123.456')

    def test_get_fftparams(self):
        fftparams = utils.get_fftparams('L1:TEST-CHANNEL')
        self.assertIsInstance(fftparams, utils.FftParams)
        for key in utils.FFT_PARAMS.keys():
            self.assertIsNone(getattr(fftparams, key))
        fftparams = utils.get_fftparams('L1:TEST-CHANNEL', window='hanning',
                                        overlap=0)
        self.assertEqual(fftparams.window, 'hanning')
        self.assertEqual(fftparams.overlap, 0)
        self.assertRaises(ZeroDivisionError, utils.get_fftparams,
                          None, stride=0)

    def test_parse_math_definition(self):
        chans, operators = mathutils.parse_math_definition(
            "L1:TEST*2 + L1:TEST2^5")
        self.assertEqual(len(chans), 2)
        self.assertListEqual(operators, [operator.add])
        self.assertIsInstance(chans[0], tuple)
        self.assertEqual(chans[0][0], 'L1:TEST')
        self.assertIsInstance(chans[0][1], list)
        self.assertEqual(len(chans[0][1]), 1)
        self.assertTupleEqual(chans[0][1][0], (operator.mul, 2))
        self.assertIsInstance(chans[1], tuple)
        self.assertEqual(chans[1][0], 'L1:TEST2')
        self.assertTupleEqual(chans[1][1][0], (operator.pow, 5))

    # -- test add/get methods -------------------

    def test_add_timeseries(self):
        a = TimeSeries([1, 2, 3, 4, 5], name='test name', epoch=0,
                       sample_rate=1)
        # test simple add using 'name'
        data.add_timeseries(a)
        self.assertIn('test name', globalv.DATA)
        self.assertEqual(globalv.DATA['test name'], [a])
        # test add using key kwarg
        data.add_timeseries(a, key='test key')
        self.assertIn('test key', globalv.DATA)
        self.assertEqual(globalv.DATA['test key'], [a])
        # test add to existing key with coalesce
        b = TimeSeries([6, 7, 8, 9, 10], name='test name 2', epoch=5,
                       sample_rate=1)
        data.add_timeseries(b, key='test key', coalesce=True)
        self.assertEqual(globalv.DATA['test key'],
                         [a.append(b, inplace=False)])

    def test_get_timeseries(self):
        # test simple get after add
        a = TimeSeries([1, 2, 3, 4, 5], name='test name', epoch=0,
                       sample_rate=1)
        data.add_timeseries(a)
        b = data.get_timeseries('test name', [(0, 5)])
        self.assertEqual(a, b)
        # test more complicated add with a cache
        a = data.get_timeseries('H1:LOSC-STRAIN', LOSC_SEGMENTS,
                                cache=self.FRAMES['H1:LOSC-STRAIN'])
        b = data.get_timeseries('H1:LOSC-STRAIN', LOSC_SEGMENTS)
        self.assertEqual(a, b)

    @empty_globalv_CHANNELS
    def test_get_spectrogram(self):
        self.assertRaises(TypeError, data.get_spectrogram, 'H1:LOSC-STRAIN',
                          LOSC_SEGMENTS, cache=self.FRAMES['H1:LOSC-STRAIN'])
        a = data.get_spectrogram('H1:LOSC-STRAIN', LOSC_SEGMENTS,
                                 cache=self.FRAMES['H1:LOSC-STRAIN'],
                                 stride=4, fftlength=2, overlap=1)

    def test_get_spectrum(self):
        a = data.get_spectrum('H1:LOSC-STRAIN', LOSC_SEGMENTS,
                              cache=self.FRAMES['H1:LOSC-STRAIN'])
        a = data.get_spectrum('H1:LOSC-STRAIN', LOSC_SEGMENTS, format='asd',
                              cache=self.FRAMES['H1:LOSC-STRAIN'])

    def test_get_coherence_spectrogram(self):
        cache = Cache([e for c in self.FRAMES for e in self.FRAMES[c]])
        a = data.get_coherence_spectrogram(
            ('H1:LOSC-STRAIN', 'L1:LOSC-STRAIN'), LOSC_SEGMENTS, cache=cache,
            stride=4, fftlength=2, overlap=1)

    def test_get_coherence_spectrum(self):
        cache = Cache([e for c in self.FRAMES for e in self.FRAMES[c]])
        a = data.get_coherence_spectrogram(
            ('H1:LOSC-STRAIN', 'L1:LOSC-STRAIN'), LOSC_SEGMENTS, cache=cache,
            stride=4, fftlength=2, overlap=1)
