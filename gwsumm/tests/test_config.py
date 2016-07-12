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

import tempfile

from matplotlib import rcParams

from astropy import units

from common import unittest
from gwsumm import (state, config, html)
from gwsumm.channels import get_channel

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


class ConfigTestCase(unittest.TestCase):

    def new(self):
        return config.GWSummConfigParser()

    def test_init(self):
        cp = self.new()
        self.assertIs(cp.optionxform, str)
        self.assertIs(cp._dict, config.OrderedDict)

    def test_read(self):
        cp = self.new()
        self.assertRaises(IOError, cp.read, 'does-not-exist.ini')

    def test_set_date_options(self):
        cp = self.new()
        cp.set_date_options(0, 100)

    def test_get_css(self):
        cp = self.new()
        css = cp.get_css()
        self.assertEqual(css, html.get_css().values())

    def test_get_javascript(self):
        cp = self.new()
        js = cp.get_javascript()
        self.assertEqual(js, html.get_js().values())

    def test_load_rcParams(self):
        cp = self.new()
        cp.add_section('rcParams')
        cp.set('rcParams', 'axes.labelsize', '100')
        cp.load_rcParams()
        self.assertEqual(rcParams['axes.labelsize'], 100)

    def test_load_states(self):
        cp = self.new()
        cp.set_date_options(0, 100)
        cp.add_section('states')
        cp.set('states', 'locked', 'X1:TEST-STATE:1')
        cp.load_states()
        states = state.get_states()
        self.assertEqual(len(states), 2)
        self.assertIn('locked', states)
        self.assertEqual(states['locked'].definition, 'X1:TEST-STATE:1')
        self.assertIn(state.ALLSTATE, states)

    def test_finalize(self):
        cp = self.new()
        cp.set_date_options(0, 100)
        cp.finalize()

    def test_load_plugins(self):
        cp = self.new()
        cp.add_section('plugins')
        cp.set('plugins', 'tempfile', '')
        plugins = cp.load_plugins()
        self.assertListEqual(plugins, [tempfile])

    def test_load_units(self):
        cp = self.new()
        cp.add_section('units')
        cp.set('units', 'myunit', 'meter')
        newunits = cp.load_units()
        self.assertListEqual(newunits, [units.meter])

    def test_load_channels(self):
        cp = self.new()
        cp.add_section('X1:TEST-CHANNEL')
        cp.set('X1:TEST-CHANNEL', 'frametype', 'X1_TEST')
        cp.load_channels()
        c = get_channel('X1:TEST-CHANNEL')
        self.assertEqual(c.frametype, 'X1_TEST')
        # test with interpolation
        cp.set(config.DEFAULTSECT, 'ifo', 'X1')
        cp.add_section('%(ifo)s:TEST-CHANNEL_2')
        cp.set('%(ifo)s:TEST-CHANNEL_2', 'resample', '128')
        cp.interpolate_section_names(ifo='X1')
        cp.load_channels()
        c = get_channel('X1:TEST-CHANNEL_2')
        self.assertEqual(c.resample, 128)
