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

"""Tests for :mod:`gwsumm.config`
"""

import os.path
import tempfile
from collections import OrderedDict
from configparser import (DEFAULTSECT, ConfigParser)

from six.moves import StringIO

import pytest

from matplotlib import rcParams

from astropy import units

from gwsumm import (state, config, html)
from gwsumm.channels import get_channel

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

TEST_CONFIG = StringIO("""
[DEFAULT]
defaultoption = defaultvalue

[section]
option1 = value1
option2 = True
option3 = 4

[plugins]
tempfile = ''

[units]
myunit = meter
cochrane = dimensionless

[%(IFO)s]
""")


def assert_configparser_equal(a, b):
    for sect in set([DEFAULTSECT] + list(a.sections()) +
                    list(b.sections())):
        assert list(a.items(sect)) == list(b.items(sect))


class TestGWSummConfigParser(object):
    PARSER = config.GWSummConfigParser

    @classmethod
    def new(cls):
        TEST_CONFIG.seek(0)
        cp = cls.PARSER()
        cp.read_file(TEST_CONFIG)
        TEST_CONFIG.seek(0)
        return cp

    @classmethod
    @pytest.fixture()
    def cnfg(cls):
        return cls.new()

    # -- test creation --------------------------

    def test_init(self):
        cp = self.new()
        assert cp.optionxform is str
        assert cp._dict is OrderedDict

    # -- test methods ---------------------------

    def test_ndoptions(self, cnfg):
        ndopts = cnfg.ndoptions('section')
        assert isinstance(ndopts, list)
        assert 'defaultoption' not in ndopts

    def test_nditems(self, cnfg):
        nditms = cnfg.nditems('section')
        assert isinstance(nditms, list)
        assert ('defaultoption', 'defaultvalue') not in nditms

    def test_read(self):
        cp = self.new()

        # read config from file
        with tempfile.NamedTemporaryFile(mode='w') as f:
            f.write(TEST_CONFIG.read())
            TEST_CONFIG.seek(0)  # rewind for other users
            read_ = cp.read(f.name)
            assert read_ == [f.name]
            assert cp.files == [os.path.abspath(f.name)]

        # check error gets raised when file isn't read
        with pytest.raises(IOError):
            cp.read('does-not-exist.ini')

    def test_from_configparser(self, cnfg):
        # check that GWSummConfigParser gets returned as is
        copy = self.PARSER.from_configparser(cnfg)
        assert copy is cnfg

        # check that base ConfigParser gets converted to GWSummConfigParser
        cp = ConfigParser()
        try:
            cp.read_file(TEST_CONFIG)
        except AttributeError:
            cp.readfp(TEST_CONFIG)
        TEST_CONFIG.seek(0)
        copy = self.PARSER.from_configparser(cp)
        assert isinstance(copy, self.PARSER)
        print(list(copy.sections()))
        print(list(cnfg.sections()))
        assert_configparser_equal(copy, cnfg)

    def test_interpolate_section_names(self, cnfg):
        assert 'X1' not in cnfg.sections()
        assert '%(IFO)s' in cnfg.sections()
        cnfg.interpolate_section_names(IFO='X1')
        assert 'X1' in cnfg.sections()
        assert '%(IFO)s' not in cnfg.sections()

    @pytest.mark.parametrize('ifo, obs, exp', [
        ('L1', None, 'LIGO Livingston'),
        ('X1', 'Einstein Telescope', 'Einstein Telescope'),
    ])
    def test_set_ifo_options(self, ifo, obs, exp):
        cp = self.new()
        cp.set_ifo_options(ifo, observatory=obs)
        assert cp.get(DEFAULTSECT, 'IFO') == ifo.upper()
        assert cp.get(DEFAULTSECT, 'ifo') == ifo.lower()
        assert cp.get(DEFAULTSECT, 'SITE') == ifo[0].upper()
        assert cp.get(DEFAULTSECT, 'site') == ifo[0].lower()
        assert cp.get(DEFAULTSECT, 'observatory') == exp

    def test_set_date_options(self):
        cp = self.new()
        cp.set_date_options(0, 100)
        assert cp.get(DEFAULTSECT, 'gps-start-time') == '0'
        assert cp.get(DEFAULTSECT, 'gps-end-time') == '100'
        assert cp.get(DEFAULTSECT, 'yyyy') == '1980'
        assert cp.get(DEFAULTSECT, 'duration') == '100'

    def test_load_rcParams(self):
        # check empty config doesn't cause havoc
        cp = self.PARSER()
        assert cp.load_rcParams() == {}

        cp = self.new()
        cp.add_section('rcParams')
        cp.set('rcParams', 'axes.labelsize', '100')
        new = cp.load_rcParams()
        assert new == {'axes.labelsize': 100}
        assert rcParams['axes.labelsize'] == 100

    def test_load_states(self):
        cp = self.new()
        cp.set_date_options(0, 100)
        cp.add_section('states')
        cp.set('states', 'locked', 'X1:TEST-STATE:1')
        cp.load_states()
        states = state.get_states()
        assert len(states) == 2
        assert 'locked' in states
        assert states['locked'].definition == 'X1:TEST-STATE:1'
        assert state.ALLSTATE in states

    def test_load_plugins(self, cnfg):
        # check that empty section doesn't cause havoc
        cp = self.PARSER()
        assert cp.load_plugins() == []

        # check plugins get laoded
        plugins = cnfg.load_plugins()
        assert plugins == [tempfile]

    def test_load_units(self, cnfg):
        # check that empty section doesn't cause havoc
        cp = self.PARSER()
        assert cp.load_units() == []

        newunits = cnfg.load_units()
        assert newunits == [units.meter, units.dimensionless_unscaled]

    def test_load_channels(self):
        # test simple channel section
        cp = self.PARSER()
        cp.add_section('X1:TEST-CHANNEL')
        cp.set('X1:TEST-CHANNEL', 'frametype', 'X1_TEST')
        cp.load_channels()
        c = get_channel('X1:TEST-CHANNEL')
        assert c.frametype == 'X1_TEST'

        # test with interpolation
        cp.set(DEFAULTSECT, 'IFO', 'X1')
        cp.add_section('%(IFO)s:TEST-CHANNEL_2')
        cp.set('%(IFO)s:TEST-CHANNEL_2', 'resample', '128')
        cp.interpolate_section_names(IFO='X1')
        cp.load_channels()
        c = get_channel('X1:TEST-CHANNEL_2')
        assert c.resample == 128

        # test bit parsing
        cp.set('X1:TEST-CHANNEL', '0', 'Bit 0')
        cp.set('X1:TEST-CHANNEL', '1', r'r"A\_B"')
        cp.load_channels()
        c = get_channel('X1:TEST-CHANNEL')
        assert c.bits == ['Bit 0', r'A\_B']

        # test channels section
        cp.add_section('channels-test')
        cp.set('channels-test', 'channels',
               'X1:TEST-CHANNEL,X1:TEST-CHANNEL_2')
        cp.set('channels-test', 'unit', 'urad')
        cp.load_channels()
        assert c.unit == units.microradian

    def test_finalize(self):
        cp = self.new()
        cp.set_date_options(0, 100)
        cp.finalize()

    def test_get_css(self):
        # check empty result returns defaults
        cp = self.PARSER()
        css = cp.get_css()
        assert css == list(html.get_css().values())

        # check overrides
        cp.add_section('html')
        cp.set('html', 'bootstrap-css', 'test.css')
        css = cp.get_css()
        print(css)
        assert 'test.css' in css
        assert html.get_css()['bootstrap'] not in css

        # check custom files
        cp.set('html', 'extra-css', '"extra.css","/static/extra2.css"')
        css = cp.get_css()
        assert 'test.css' in css  # still finds overrides
        assert 'extra.css' in css and '/static/extra2.css' in css

    def test_get_javascript(self):
        # check empty result returns defaults
        cp = self.PARSER()
        js = cp.get_javascript()
        assert js == list(html.get_js().values())

        # check overrides
        cp.add_section('html')
        cp.set('html', 'bootstrap-js', 'test.js')
        js = cp.get_javascript()
        assert 'test.js' in js
        assert html.get_js()['bootstrap'] not in js

        # check custom files
        cp.set('html', 'extra-js', '"extra.js","/static/extra2.js"')
        js = cp.get_javascript()
        assert 'test.js' in js  # still finds overrides
        assert 'extra.js' in js and '/static/extra2.js' in js
