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

from six import string_types

import pytest

try:
    from unittest import mock
except ImportError:
    import mock

from gwsumm import (utils, globalv)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


def test_elapsed_time():
    e = time.time() - globalv.START
    assert utils.elapsed_time() - e < .1


def test_vprint(capsys):
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
    assert re.match(r'\Aanything \(\d+\.\d\d\)\n\Z', out) is not None


def test_mkdir():
    d = 'test-dir/test-dir2'
    try:
        utils.mkdir(d)
        assert os.path.isdir(d)
    finally:
        if os.path.isdir(d):
            shutil.rmtree(d)


def test_nat_sorted():
    # sorted strings numerically
    assert utils.nat_sorted(['1', '10', '2', 'a', 'B']) == [
        '1', '2', '10', 'B', 'a']


def test_which():
    environ = os.environ.copy()
    os.environ['PATH'] = '/usr/bin:/bin'
    try:
        with mock.patch('os.path.isfile', return_value=True), \
                mock.patch('os.access', return_value=True):
            assert utils.which('/bin/bash') == '/bin/bash'
            assert utils.which('bash') == '/usr/bin/bash'
        with mock.patch('os.path.isfile', return_value=False):
            assert utils.which('blah') is None
    finally:
        os.environ = environ


@pytest.mark.parametrize('chan, mask', [
    ('L1:TEST-ODC_CHANNEL_OUT_DQ', 'L1:TEST-ODC_CHANNEL_BITMASK'),
    ('L1:TEST-ODC_CHANNEL_OUTMON', 'L1:TEST-ODC_CHANNEL_BITMASK'),
    ('L1:TEST-ODC_CHANNEL_LATCH', 'L1:TEST-ODC_CHANNEL_BITMASK'),
    ('L1:TEST-CHANNEL', 'L1:TEST-CHANNEL')
])
def test_get_odc_bitmask(chan, mask):
    assert utils.get_odc_bitmask(chan) == mask


@pytest.mark.parametrize('value, out', [
    ('my random content', 'my random content'),
    ('1', 1),
    ('1.', 1.),
    ('1,', (1,)),
    ('1,2,\'test\',4', (1, 2, 'test', 4)),
    ('[], [0], 1/(2*pi)', ([], [0], 1/(2*pi))),
    ('lambda x: x**2', lambda x: x ** 2),
    (pytest, pytest),
])
def test_safe_eval(value, out):
    evalue = utils.safe_eval(value)
    assert type(evalue) is type(out)

    if not isinstance(value, string_types):
        assert evalue is out
    elif callable(out):
        assert evalue(4) == out(4)
    else:
        assert evalue == out


def test_safe_eval_2():
    # test unsafe
    with pytest.raises(ValueError) as exc:
        utils.safe_eval("os.remove('file-that-doesnt-exist')")
    assert str(exc.value).startswith('Will not evaluate string containing')

    with pytest.raises(ValueError):
        utils.safe_eval("lambda x: shutil.remove('file-that-doesnt-exist')")

    # test locals or globals
    assert utils.safe_eval('test', globals_={'test': 4}) == 4
    assert utils.safe_eval('type(self)',
                           locals_={'self': globalv}) == type(globalv)


@pytest.mark.parametrize('ifo, host', [
    ('G1', 'host.atlas.aei.uni-hannover.de'),
    ('H1', 'host.ligo-wa.caltech.edu'),
    ('L1', 'host.ligo-la.caltech.edu'),
    ('V1', 'host.virgo.ego.it'),
])
def test_get_default_ifo(ifo, host):
    assert utils.get_default_ifo(host) == ifo

    with pytest.raises(ValueError):
        utils.get_default_ifo('host.ligo.caltech.edu')
