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

"""Tests for `gwsumm.utils`

"""

import os.path
import time
import re
import sys
import shutil
from math import pi

from common import unittest
from gwsumm import (utils, globalv)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


class UtilsTestCase(unittest.TestCase):
    def test_elapsed_time(self):
        e = time.time() - globalv.START
        self.assertTrue(utils.elapsed_time() - e < 1)

    def test_get_odc_bitmask(self):
        # test fast ODC
        self.assertEqual(
            utils.get_odc_bitmask('L1:TEST-ODC_CHANNEL_OUT_DQ'),
            'L1:TEST-ODC_CHANNEL_BITMASK')
        # test slow ODC
        self.assertEqual(
            utils.get_odc_bitmask('L1:TEST-ODC_CHANNEL_OUTMON'),
            'L1:TEST-ODC_CHANNEL_BITMASK')
        # test LATCH
        self.assertEqual(
            utils.get_odc_bitmask('L1:TEST-ODC_CHANNEL_LATCH'),
            'L1:TEST-ODC_CHANNEL_BITMASK')
        # test other (noop)
        self.assertEqual(
            utils.get_odc_bitmask('L1:TEST-CHANNEL_NAME'),
            'L1:TEST-CHANNEL_NAME')

    def test_get_default_ifo(self):
        self.assertEqual(
            utils.get_default_ifo('host.ligo-la.caltech.edu'), 'L1')
        self.assertEqual(
            utils.get_default_ifo('host.ligo-wa.caltech.edu'), 'H1')
        self.assertEqual(
            utils.get_default_ifo('host.atlas.aei.uni-hannover.de'), 'G1')
        self.assertEqual(
            utils.get_default_ifo('host.virgo.ego.it'), 'V1')
        self.assertRaises(ValueError, utils.get_default_ifo,
                          'host.ligo.caltech.edu')

    def test_safe_eval(self):
        # test normal string
        self.assertEqual(utils.safe_eval('test'), 'test')
        self.assertEqual(utils.safe_eval('my random content'),
                         'my random content')
        # test float/int
        t = utils.safe_eval('1')
        self.assertEqual(t, 1)
        self.assertIsInstance(t, int)
        t = utils.safe_eval('1.')
        self.assertEqual(t, 1.)
        self.assertIsInstance(t, float)
        # test tuple/list
        t = utils.safe_eval('1,')
        self.assertEqual(t, (1,))
        self.assertIsInstance(t, tuple)
        t = utils.safe_eval('1,2,\'test\',4')
        self.assertEqual(t, (1, 2, 'test', 4))
        # test with math
        t = utils.safe_eval("[], [0], 1/(2*pi)")
        self.assertIsInstance(t, tuple)
        self.assertEqual(t, ([], [0], 1/(2*pi)))
        # test lambda
        t = utils.safe_eval("lambda x: x**2")
        self.assertTrue(callable(t))
        self.assertEqual(t(4), 16)
        # test unsafe
        self.assertRaises(ValueError, utils.safe_eval,
                          "os.remove('file-that-doesnt-exist')")
        self.assertRaises(ValueError, utils.safe_eval,
                          "lambda x: shutil.remove('file-that-doesnt-exist')")
        # test locals or globals
        t = utils.safe_eval('test', globals_={'test': 4})
        self.assertEqual(t, 4)
        t = utils.safe_eval('type(self)', locals_={'self': self})
        self.assertEqual(t, type(self))

    def test_mkdir(self):
        d = 'test-dir/test-dir2'
        try:
            utils.mkdir(d)
        finally:
            if os.path.isdir(d):
                shutil.rmtree(d)

    def test_nat_sorted(self):
        # sorted strings numerically
        self.assertListEqual(
            utils.nat_sorted(['1', '10', '2']), ['1', '2', '10'])


class TestVprint(object):
    # this is separate from the other tests so we can use pytest's
    # capsys fixture to check stdout for content when testing vprint()

    def test_vprint(self, capsys):
        # test non-verbose
        globalv.VERBOSE = False
        utils.vprint('anything', stream=sys.stdout)
        out, err = capsys.readouterr()
        assert out == ''
        # test verbose
        globalv.VERBOSE = True
        utils.vprint('anything', stream=sys.stdout)
        out, err = capsys.readouterr()
        assert out == 'anything'
        # test profiled
        globalv.PROFILE = True
        utils.vprint('anything\n', stream=sys.stdout)
        out, err = capsys.readouterr()
        assert re.match('\Aanything \(\d+\.\d\d\)\n\Z', out) is not None
