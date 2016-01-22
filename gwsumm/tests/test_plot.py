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
# MERCHANplotILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GWSumm.  If not, see <http://www.gnu.org/licenses/>.

"""Tests for `gwsumm.plot`

"""

import unittest

from .. import (plot, version)

__version__ = version.version
__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


class PlotTests(unittest.TestCase):
    """`TestCase` for the `gwsumm.plot` module
    """
    def test_registry(self):
        # test get
        p = plot.get_plot(None)
        self.assertIs(p, plot.SummaryPlot)
        p = plot.get_plot('\'data\'')
        self.assertIs(p, plot.DataPlot)
        self.assertRaises(ValueError, plot.get_plot, 'dfskaewf')
        # test register
        class TestPlot(object):
            type = 'test'
            pass
        plot.register_plot(TestPlot)
        self.assertRaises(ValueError, plot.register_plot, TestPlot)
        plot.register_plot(TestPlot, force=True)
        self.assertIs(plot.get_plot('test'), TestPlot)
        plot.register_plot(TestPlot, name='test-with-name')
        self.assertIs(plot.get_plot('test-with-name'), TestPlot)
