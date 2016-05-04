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

"""Tests for `gwsumm.tabs`

"""

from gwpy.detector import Channel

from .compat import unittest
from .. import tabs

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


class TabTests(unittest.TestCase):
    """`TestCase` for the `gwsumm.tabs` module
    """
    def test_registry(self):
        # test get
        tab = tabs.get_tab('basic')
        self.assertIs(tab, tabs.Tab)
        tab = tabs.get_tab('\'plots\'')
        self.assertIs(tab, tabs.PlotTab)
        self.assertRaises(ValueError, tabs.get_tab, 'dfskaewf')
        # test register
        class TestTab(object):
            type = 'test'
            pass
        tabs.register_tab(TestTab)
        self.assertRaises(ValueError, tabs.register_tab, TestTab)
        tabs.register_tab(TestTab, force=True)
        self.assertIs(tabs.get_tab('test'), TestTab)
        tabs.register_tab(TestTab, name='test-with-name')
        self.assertIs(tabs.get_tab('test-with-name'), TestTab)
