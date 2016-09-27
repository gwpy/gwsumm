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

import os.path

import pytest

from gwsumm import tabs
from gwsumm.plot import SummaryPlot

from common import unittest

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


class TabRegistryTestCase(unittest.TestCase):
    """`TestCase` for the `gwsumm.tabs` registry
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

# -- test tab classes ---------------------------------------------------------


class TabTestCase(unittest.TestCase):
    TYPE = 'basic'
    DEFAULT_ARGS = ['Test']

    @classmethod
    def setUpClass(cls):
        cls.TAB = tabs.get_tab(cls.TYPE)

    def _create(self, *args, **kwargs):
        args = list(args)
        while len(args) < len(self.DEFAULT_ARGS):
            args.append(self.DEFAULT_ARGS[len(args)])
        print(args)
        return self.TAB(*args, **kwargs)

    def test_init(self):
        self._test_init('Test')

    def _test_init(self, *args, **kwargs):
        if len(args) == 0:
           args = self.DEFAULT_ARGS
        # test basic creation and defaults
        tab = self._create(*args, **kwargs)
        self.assertEqual(tab.type, self.TYPE)
        self.assertEqual(tab.name, args[0])
        self.assertEqual(tab.shortname, kwargs.pop('shortname', tab.name))
        self.assertEqual(tab.children, kwargs.pop('children', []))
        self.assertEqual(tab.parent, kwargs.pop('parent', None))
        self.assertEqual(tab.group, kwargs.pop('group', None))
        self.assertEqual(tab.path, kwargs.pop('path', os.path.curdir))
        self.assertEqual(tab.hidden, kwargs.pop('hidden', False))
        return tab

    def test_shortname(self):
        tab = self._create()
        self.assertEqual(tab.shortname, tab.name)
        tab = self._create('Test', shortname='ShortTest')
        self.assertEqual(tab.shortname, 'ShortTest')

    def test_index(self):
        tab = self._create()
        self.assertEqual(tab.index, os.path.join('test', 'index.html'))
        tab2 = self._create('Parent')
        del tab.index
        tab.set_parent(tab2)
        self.assertEqual(tab.index,
                         os.path.join('parent', 'test', 'index.html'))


# -- external tab

class ExternalTabTestCase(TabTestCase):
    TYPE = 'external'
    DEFAULT_ARGS = ['Test', '//test.com']

    def test_init(self):
        tab = self._test_init()
        self.assertEqual(tab.url, '//test.com')


# -- plot tab

class PlotTabTestCase(TabTestCase):
    TYPE = 'plots'

    def test_init(self):
        plots = ['test.png']
        tab = self._test_init('Test', plots=plots)
        self.assertListEqual(tab.plots, map(SummaryPlot, plots))
        self.assertIsNone(tab.layout)

    def test_add_plot(self):
        tab = self._create()
        before = tab.plots[:]
        plot = 'test.png'
        tab.add_plot(plot)
        self.assertListEqual(tab.plots, before + [SummaryPlot(href=plot)])

    def test_layout(self):
        tab = self._create()
        tab.set_layout(1)
        self.assertListEqual(tab.layout, [1])
        tab.set_layout((1, 2))
        self.assertListEqual(tab.layout, [1, 2])
        tab.set_layout((1, (1, 2)))
        self.assertListEqual(tab.layout, [1, (1, 2)])
        self.assertRaises(ValueError, tab.set_layout, 'test')
        self.assertRaises(ValueError, tab.set_layout, [1, (1, 2, 1)])
        with pytest.warns(DeprecationWarning):
            tab.layout = [1]
