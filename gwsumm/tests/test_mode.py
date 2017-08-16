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

"""Test suite for the mode module

"""

import datetime

import pytest

from gwsumm import (globalv, mode)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

DEFAULT_MODE = mode.get_mode()


def teardown_module():
    """Undo any set_mode() operations from this module
    """
    mode.set_mode(DEFAULT_MODE)


def test_get_mode():
    assert mode.get_mode().value == globalv.MODE
    assert mode.get_mode(10) == mode.Mode.day
    assert mode.get_mode('week') == mode.Mode.week
    with pytest.raises(ValueError):
        mode.get_mode('invalid mode')


def test_set_mode():
    mode.set_mode(0)
    assert globalv.MODE == mode.Mode(0).value

    mode.set_mode('GPS')
    assert globalv.MODE == mode.Mode.gps.value

    with pytest.raises(ValueError):
        mode.set_mode('invalid mode')


@pytest.mark.parametrize('m, basestr', [
    ('day', 'day/20150914'),
    ('week', 'week/20150914'),
    ('month', 'month/201509'),
    ('year', 'year/2015'),
])
def test_get_base(m, basestr):
    date = datetime.date(2015, 9, 14)
    mode.set_mode(m)
    assert mode.get_base(date) == basestr
