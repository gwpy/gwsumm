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
from collections import OrderedDict

from six.moves.urllib.request import urlopen

import pytest

from numpy import (arange, testing as nptest)

from lal.utils import CacheEntry

from glue.lal import Cache

from gwpy.timeseries import TimeSeries
from gwpy.detector import Channel
from gwpy.segments import (Segment, SegmentList)

from gwsumm import (data, globalv)
from gwsumm.data import (utils, mathutils)

from .common import empty_globalv_CHANNELS

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
    with open(target, 'wb') as f:
        f.write(response.read())
    return target


class TestData(object):
    """Tests for :mod:`gwsumm.data`:
    """
    @classmethod
    def setup_class(cls):
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
    def teardown_class(cls):
        # remove the temporary data
        shutil.rmtree(cls._tempdir)

    # -- test utilities -------------------------

    def test_find_frame_type(self):
        channel = Channel('L1:TEST-CHANNEL')
        assert data.find_frame_type(channel) == 'L1_R'

        channel = Channel('C1:TEST-CHANNEL')
        assert data.find_frame_type(channel) == 'R'

        channel = Channel('H1:TEST-CHANNEL.rms,s-trend')
        assert data.find_frame_type(channel) == 'H1_T'

        channel = Channel('H1:TEST-CHANNEL.rms,m-trend')
        assert data.find_frame_type(channel) == 'H1_M'

        channel = Channel('H1:TEST-CHANNEL.rms,reduced')
        assert data.find_frame_type(channel) == 'H1_LDAS_C02_L2'

        channel = Channel('H1:TEST-CHANNEL.rms,online')
        assert data.find_frame_type(channel) == 'H1_lldetchar'

    def test_get_channel_type(self):
        assert data.get_channel_type('L1:TEST-CHANNEL') == 'adc'
        assert data.get_channel_type('G1:DER_DATA_HL') == 'proc'
        assert data.get_channel_type('H1:GDS-CALIB_STRAIN') == 'proc'
        assert data.get_channel_type('V1:GDS-CALIB_STRAIN') == 'adc'

    @empty_globalv_CHANNELS
    def test_make_globalv_key(self):
        fftparams = utils.get_fftparams(
            'L1:TEST-CHANNEL',
            stride=123.456,
            window='test-window',
            method='scipy-welch',
        )
        key = utils.make_globalv_key('L1:TEST-CHANNEL', fftparams)
        assert key == ';'.join([
            'L1:TEST-CHANNEL',  # channel
            'scipy-welch',  # method
            '',  # fftlength
            '',  # overlap
            'test-window',  # window
            '123.456',  # stride
            '',  # FFT scheme
        ])

    def test_get_fftparams(self):
        fftparams = utils.get_fftparams('L1:TEST-CHANNEL')
        assert isinstance(fftparams, utils.FftParams)

        for key in utils.FFT_PARAMS:
            assert (getattr(fftparams, key) is
                    utils.DEFAULT_FFT_PARAMS.get(key, None))

        fftparams = utils.get_fftparams('L1:TEST-CHANNEL', window='hanning',
                                        overlap=0)
        assert fftparams.window == 'hanning'
        assert fftparams.overlap == 0

        with pytest.raises(ZeroDivisionError):
            utils.get_fftparams(None, stride=0)

    @pytest.mark.parametrize('definition, math', [
        (
             'L1:TEST + L1:TEST2',
             ([('L1:TEST', 'L1:TEST2'), ([], [])], [operator.add]),
        ),
        (
             'L1:TEST + L1:TEST2 * 2',
             ([('L1:TEST', 'L1:TEST2'), ([], [(operator.mul, 2)])],
              [operator.add]),
        ),
        (
             'L1:TEST * 2 + L1:TEST2 ^ 5',
             ([('L1:TEST', 'L1:TEST2'),
               ([(operator.mul, 2)], [(operator.pow, 5)])],
              [operator.add]),
        ),
    ])
    def test_parse_math_definition(self, definition, math):
        chans, operators = mathutils.parse_math_definition(definition)
        assert chans == OrderedDict(list(zip(*math[0])))
        assert operators == math[1]

    # -- test add/get methods -------------------

    def test_add_timeseries(self):
        a = TimeSeries([1, 2, 3, 4, 5], name='test name', epoch=0,
                       sample_rate=1)

        # test simple add using 'name'
        data.add_timeseries(a)
        assert 'test name' in globalv.DATA
        assert len(globalv.DATA['test name']) == 1
        assert globalv.DATA['test name'][0] is a

        # test add using key kwarg
        data.add_timeseries(a, key='test key')
        assert globalv.DATA['test key'][0] is a

        # test add to existing key with coalesce
        b = TimeSeries([6, 7, 8, 9, 10], name='test name 2', epoch=5,
                       sample_rate=1)
        data.add_timeseries(b, key='test key', coalesce=True)
        assert len(globalv.DATA['test key']) == 1
        nptest.assert_array_equal(globalv.DATA['test key'][0].value,
                                  arange(1, 11))

    def test_get_timeseries(self):
        # empty globalv.DATA
        globalv.DATA = type(globalv.DATA)()

        # test simple get after add
        a = TimeSeries([1, 2, 3, 4, 5], name='test name', epoch=0,
                       sample_rate=1)
        data.add_timeseries(a)
        b, = data.get_timeseries('test name', [(0, 5)], nproc=1)
        nptest.assert_array_equal(a.value, b.value)
        assert a.sample_rate.value == b.sample_rate.value

        # test more complicated add with a cache
        a, = data.get_timeseries('H1:LOSC-STRAIN', LOSC_SEGMENTS,
                                 cache=self.FRAMES['H1:LOSC-STRAIN'],
                                 nproc=1)
        b, = data.get_timeseries('H1:LOSC-STRAIN', LOSC_SEGMENTS,
                                 nproc=1)
        nptest.assert_array_equal(a.value, b.value)

    @empty_globalv_CHANNELS
    def test_get_spectrogram(self):
        with pytest.raises(TypeError):
            data.get_spectrogram('H1:LOSC-STRAIN', LOSC_SEGMENTS,
                                 cache=self.FRAMES['H1:LOSC-STRAIN'],
                                 nproc=1)
        data.get_spectrogram('H1:LOSC-STRAIN', LOSC_SEGMENTS,
                             cache=self.FRAMES['H1:LOSC-STRAIN'],
                             stride=4, fftlength=2, overlap=1,
                             nproc=1)

    def test_get_spectrum(self):
        a, _, _ = data.get_spectrum('H1:LOSC-STRAIN', LOSC_SEGMENTS,
                                    cache=self.FRAMES['H1:LOSC-STRAIN'],
                                    nproc=1)
        b, _, _ = data.get_spectrum('H1:LOSC-STRAIN', LOSC_SEGMENTS,
                                    format='asd',
                                    cache=self.FRAMES['H1:LOSC-STRAIN'],
                                    nproc=1)
        nptest.assert_array_equal(a.value ** (1/2.), b.value)

    def test_get_coherence_spectrogram(self):
        cache = Cache([e for c in self.FRAMES for e in self.FRAMES[c]])
        data.get_coherence_spectrogram(
            ('H1:LOSC-STRAIN', 'L1:LOSC-STRAIN'), LOSC_SEGMENTS, cache=cache,
            stride=4, fftlength=2, overlap=1, nproc=1,
        )

    def test_get_coherence_spectrum(self):
        cache = Cache([e for c in self.FRAMES for e in self.FRAMES[c]])
        data.get_coherence_spectrogram(
            ('H1:LOSC-STRAIN', 'L1:LOSC-STRAIN'), LOSC_SEGMENTS, cache=cache,
            stride=4, fftlength=2, overlap=1, nproc=1,
        )
