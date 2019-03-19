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

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


# -- gwpy.tabs.registry -------------------------------------------------------

@pytest.mark.parametrize('name, tab', [
    ('basic', tabs.Tab),
    ('\'plots\'', tabs.PlotTab),
])
def test_get_tab(name, tab):
    assert tabs.get_tab(name) is tab
    with pytest.raises(ValueError):
        tabs.get_tab('fasmklwea')


def test_register_tab():
    class TestTab(object):
        type = 'test'
        pass

    tabs.register_tab(TestTab)
    with pytest.raises(ValueError):
        tabs.register_tab(TestTab)
    tabs.register_tab(TestTab, force=True)

    assert tabs.get_tab('test') is TestTab

    tabs.register_tab(TestTab, name='test-with-name')
    assert tabs.get_tab('test-with-name') is TestTab


# -- test tab classes ---------------------------------------------------------

class TestTab(object):
    TYPE = 'basic'
    DEFAULT_ARGS = ['Test']

    @classmethod
    def setup_class(cls):
        cls.TAB = tabs.get_tab(cls.TYPE)

    def create(self, *args, **kwargs):
        args = list(args)
        while len(args) < len(self.DEFAULT_ARGS):
            args.append(self.DEFAULT_ARGS[len(args)])
        return self.TAB(*args, **kwargs)

    def test_init(self):
        self._test_init('Test')

    def _test_init(self, *args, **kwargs):
        if len(args) == 0:
            args = self.DEFAULT_ARGS
        # test basic creation and defaults
        tab = self.create(*args, **kwargs)
        assert tab.type == self.TYPE
        assert tab.name == args[0]
        assert tab.shortname == kwargs.pop('shortname', tab.name)
        assert tab.children == kwargs.pop('children', [])
        assert tab.parent == kwargs.pop('parent', None)
        assert tab.group == kwargs.pop('group', None)
        assert tab.path == kwargs.pop('path', os.path.curdir)
        assert tab.hidden == kwargs.pop('hidden', False)
        return tab

    def test_shortname(self):
        tab = self.create()
        assert tab.shortname == tab.name
        tab = self.create('Test', shortname='ShortTest')
        assert tab.shortname == 'ShortTest'

    def test_index(self):
        tab = self.create()
        assert tab.index == os.path.join('test', 'index.html')
        tab2 = self.create('Parent')
        del tab.index
        tab.set_parent(tab2)
        assert tab.index == os.path.join('parent', 'test', 'index.html')


# -- external tab

class TestExternalTab(TestTab):
    TYPE = 'external'
    DEFAULT_ARGS = ['Test', '//test.com']

    def test_init(self):
        tab = self._test_init()
        assert tab.url == '//test.com'


# -- plot tab

class TestPlotTab(TestTab):
    TYPE = 'plots'

    def test_init(self):
        plots = ['test.png']
        tab = self._test_init('Test', plots=plots)
        assert tab.plots == list(map(SummaryPlot, plots))
        assert tab.layout is None

    def test_add_plot(self):
        tab = self.create()
        before = tab.plots[:]
        plot = 'test.png'
        tab.add_plot(plot)
        assert tab.plots == before + [SummaryPlot(href=plot)]

    def test_layout(self):
        tab = self.create()
        tab.set_layout(1)
        assert tab.layout == [1]

        tab.set_layout((1, 2))
        assert tab.layout == [1, 2]

        tab.set_layout((1, (1, 2)))
        assert tab.layout == [1, (1, 2)]

        with pytest.raises(ValueError):
            tab.set_layout('test')
        with pytest.raises(ValueError):
            tab.set_layout([1, (1, 2, 1)])
        with pytest.warns(DeprecationWarning):
            tab.layout = [1]
