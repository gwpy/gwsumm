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

from gwpy.detector import Channel

from .compat import unittest
from .. import data

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


class DataTests(unittest.TestCase):
    """`TestCase` for the `gwsumm.data` module
    """
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

    def test_make_key(self):
        key = data._make_key('L1:TEST-CHANNEL',
                             {'stride': 123.456, 'window': 'test-window'},
                             method='test-method', sampling=654.321)
        self.assertEqual(
            key, 'L1:TEST-CHANNEL;test-method;test-window;123.456;'
                 'None;None;654.321')

    def test_clean_fftparams(self):
        fftparams = data._clean_fftparams({}, 'L1:TEST-CHANNEL')
        self.assertDictEqual(
            fftparams,
            {'window': None, 'fftlength': 1, 'overlap': 0.5, 'stride': 2.0})
        fftparams = data._clean_fftparams(
            {'window': 'median-mean', 'overlap': 0}, 'L1:TEST-CHANNEL')
        self.assertDictEqual(
            fftparams,
            {'window': 'median-mean', 'fftlength': 1.0, 'overlap': 0.5,
             'stride': 2.0})
        self.assertRaises(ZeroDivisionError, data._clean_fftparams,
                          {'stride': 0}, None)
