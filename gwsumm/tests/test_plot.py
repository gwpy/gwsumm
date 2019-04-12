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

import os

from configparser import ConfigParser

from matplotlib import use
use('agg')  # noqa

from matplotlib import (rcParams, rc_context)

import pytest

from gwpy.detector import ChannelList
from gwpy.plot import Plot
from gwpy.plot.tex import HAS_TEX
from gwpy.segments import Segment

from gwsumm import plot as gwsumm_plot
from gwsumm.channels import get_channel

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

rcParams.update({
    'text.usetex': False,
})


# -- gwsumm.plot.registry -----------------------------------------------------

@pytest.mark.parametrize('name, plot', [
    (None, gwsumm_plot.SummaryPlot),
    ('\'data\'', gwsumm_plot.DataPlot),
])
def test_get_plot(name, plot):
    assert gwsumm_plot.get_plot(name) is plot
    with pytest.raises(ValueError):
        gwsumm_plot.get_plot('dfskaewf')


def test_registry_plot():
    class TestPlot(object):
        type = 'test'
        pass

    gwsumm_plot.register_plot(TestPlot)
    with pytest.raises(ValueError):
        gwsumm_plot.register_plot(TestPlot)
    gwsumm_plot.register_plot(TestPlot, force=True)

    assert gwsumm_plot.get_plot('test') is TestPlot

    gwsumm_plot.register_plot(TestPlot, name='test-with-name')
    assert gwsumm_plot.get_plot('test-with-name') is TestPlot


# -- gwsumm.plot.utils --------------------------------------------------------

@pytest.mark.parametrize('column, label', [
    ('test', 'Test'),
    ('rho', r'$\rho$'),
    ('frequency', 'Frequency [Hz]'),
    ('mchirp', r'Chirp mass [M$_\odot$]'),
])
def test_get_column_label(column, label):
    assert gwsumm_plot.get_column_label(column) == label


# -- gwsumm.plot.core ---------------------------------------------------------

class TestSummaryPlot(object):
    TYPE = None
    DEFAULT_ARGS = []
    DEFAULT_KWARGS = {}

    @classmethod
    def setup_class(cls):
        cls.PLOT = gwsumm_plot.get_plot(cls.TYPE)

    @classmethod
    def create(cls, *args, **kwargs):
        if args or kwargs:
            return cls.PLOT(*args, **kwargs)
        return cls.PLOT(*cls.DEFAULT_ARGS, **cls.DEFAULT_KWARGS)

    @classmethod
    @pytest.fixture()
    def plot(cls):
        return cls.create()

    def test_init(self):
        plot = self.PLOT(href='test.png')
        assert plot.href == 'test.png'
        assert plot.src == plot.href
        assert plot.new is True
        assert plot.caption == ''

        plot = self.PLOT(href='test.png', src='other.png')
        assert plot.src == 'other.png'

    @pytest.mark.parametrize('url, href', [
        (None, None),
        ('https://test.com/test.png', 'https://test.com/test.png'),
        ('test//test.png', os.path.normpath('test/test.png')),
    ])
    def test_href(self, plot, url, href):
        plot.href = url
        assert plot.href == href

    @pytest.mark.parametrize('isnew, new', [
        (True, True),
        (False, False),
        (100, True),
    ])
    def test_new(self, plot, isnew, new):
        plot.new = isnew
        assert plot.new is new

    def test_src(self, plot):
        plot.href = 'test.html'
        assert plot.src == 'test.html'

        plot.src = 'test.png'
        assert plot.src == 'test.png'

    def test_eq(self, plot):
        other = self.create()
        assert plot == other

        plot.href = 'test.png'
        other.href = 'test2.png'
        assert plot != other

        assert plot != 1

    def test_repr(self, plot):
        plot.href = 'test.png'
        assert repr(plot) == '<{0}(test.png)>'.format(self.PLOT.__name__)

    def test_str(self, plot):
        plot.href = 'test.png'
        assert str(plot) == 'test.png'


class TestDataPlot(TestSummaryPlot):
    TYPE = 'data'
    DEFAULT_ARGS = [['X1:TEST-CHANNEL', 'Y1:TEST-CHANNEL2'], 0, 100]

    def test_init(self):
        plot = self.PLOT('X1:TEST', 0, 100, blah=4)
        assert plot._channels == ['X1:TEST']
        assert plot.span == (0, 100)
        assert plot.pargs == {'blah': 4}

    # -- properties -----------------------------

    def test_span(self, plot):
        assert isinstance(plot.span, Segment)
        assert plot.span == (0, 100)

    def test_start(self, plot):
        assert plot.start == 0

    def test_end(self, plot):
        assert plot.end == 100

    def test_state(self, plot):
        assert plot.state is None

    def test_channels(self, plot):
        assert isinstance(plot.channels, ChannelList)
        assert plot.channels[0] is get_channel(self.DEFAULT_ARGS[0][0])

    def test_allchannels(self, plot):
        assert plot.allchannels == plot.channels
        plot.channels = ['X1:TEST', 'Y1:TEST', 'X1:TEST']
        assert len(plot.channels) == 3
        assert len(plot.allchannels) == 2

    def test_ifos(self, plot):
        assert plot.ifos == {'X1', 'Y1'}

    def test_tag(self, plot):
        assert plot.tag == 'MULTI_3602E9_DATA'

        plot.tag = 'TEST'
        assert plot.tag == 'TEST'

        plot.tag = None
        assert plot.tag == 'MULTI_3602E9_DATA'

    def test_outputfile(self, plot):
        assert plot.outputfile == './X1Y1-MULTI_3602E9_DATA-0-100.png'

    def test_href(self, plot):
        assert plot.href == plot.outputfile
        plot.href = 'test.png'
        assert plot.href == 'test.png'

    # -- methods --------------------------------

    def test_add_channel(self):
        plot = self.PLOT([], 0, 1)
        assert len(plot.channels) == 0
        plot.add_channel('Z1:TEST-CHANNEL')
        assert len(plot.channels) == 1
        assert plot.channels[-1] is get_channel('Z1:TEST-CHANNEL')

    def test_get_channel_groups(self):
        plot = self.PLOT(['X1:TEST.mean', 'X1:TEST.max', 'X1:TEST.min',
                          'Y1:TEST', 'G1:TEST-av', 'G1:TEST-min'], 0, 1)
        assert plot.get_channel_groups() == [
            ('X1:TEST', [get_channel('X1:TEST.mean'),
                         get_channel('X1:TEST.min'),
                         get_channel('X1:TEST.max')]),
            ('Y1:TEST', [get_channel('Y1:TEST')]),
            ('G1:TEST', [get_channel('G1:TEST-av'),
                         get_channel('G1:TEST-min')]),
        ]

    def test_from_ini(self):
        cp = ConfigParser()
        cp.add_section('plot')
        cp.set('plot', 'type', self.PLOT.type)
        cp.set('plot', 'channels', 'X1:TEST,Y1:TEST')
        cp.set('plot', 'marker', 'o')
        cp.set('plot', 'figure.figsize', '2, 1')
        plot = self.PLOT.from_ini(cp, 'plot', 0, 100)
        assert plot.span == (0, 100)
        assert plot.channels == [get_channel('X1:TEST'),
                                 get_channel('Y1:TEST')]
        assert plot.pargs == {'marker': 'o'}
        assert plot.rcParams == {'figure.figsize': (2, 1)}

        cp.set('plot', 'type', 'blah')
        with pytest.warns(UserWarning):
            self.PLOT.from_ini(cp, 'plot', 0, 100)

    def test_parse_legend_kwargs(self, plot):
        plot.pargs = {'a': 1, 'b': 1, 'legend-c': 'test'}
        assert plot.parse_legend_kwargs() == {'c': 'test'}

    def test_parse_plot_kwargs(self, plot):
        plot.pargs = {
            'a': 1,  # ignored (not a valid Axes.plot() param)
            'alpha': .8,  # mapped to each channel
            'colors': '\'red\', \'green\'',  # evaluated as tuple
            'linestyle': ['-', '--'],  # as is
            'labels': '[t.ifo for t in plot.channels]',  # eval'd with Plot
            'labels': '[t.ifo for t in plot.channels]',  # eval'd with Plot
        }
        assert plot.parse_plot_kwargs() == [
            {'alpha': .8, 'color': 'red', 'linestyle': '-', 'label': 'X1'},
            {'alpha': .8, 'color': 'green', 'linestyle': '--', 'label': 'Y1'},
        ]

    @pytest.mark.parametrize('usetex, result', [
        (False,
         [{'label': r'TEST_WITH_UNDERSCORE'},
          {'label': r'TEST_QUOTED'}]),
        pytest.param(
            True,
            [{'label': r'TEST\_WITH\_UNDERSCORE'}, {'label': r'TEST\_QUOTED'}],
            marks=pytest.mark.skipif(not HAS_TEX,
                                     reason='TeX is not available'),
        ),
    ])
    def test_parse_plot_kwargs_labels(self, plot, usetex, result):
        with rc_context(rc={'text.usetex': usetex}):
            plot.pargs = {'labels': ['TEST_WITH_UNDERSCORE', '"TEST_QUOTED"']}
            assert plot.parse_plot_kwargs() == result

    def test_parse_rcParams(self, plot):
        pargs = {
            'a': 1,
            'xtick.labelsize': '4',  # should get eval'd to int(4)
            'ytick.labelsize': 10,
        }
        rcp = plot.parse_rcParams(pargs)
        assert rcp == {'xtick.labelsize': 4, 'ytick.labelsize': 10}
        assert plot.rcParams is rcp
        assert pargs == {'a': 1}

    def test_finalize(self, plot):
        plot.pargs = {'ylim': (0, 10)}
        fig = plot.init_plot()
        ax = fig.gca()
        ax.set_title('TEST', y=1.0)
        assert ax.get_ylim() == (0, 10)

        # check fileformat PDF creates png as well
        plot.fileformat = 'pdf'
        try:
            plot.finalize(close=False)
            assert os.path.isfile(plot.outputfile)
            assert os.path.isfile(plot.outputfile.replace('.pdf', '.png'))
        finally:
            for p in [plot.outputfile,
                      plot.outputfile.replace('.pdf', '.png')]:
                if os.path.isfile(p):
                    os.remove(p)

        plot.fileformat = 'png'

        # now check normal things
        try:
            plot.finalize(close=False)
            assert ax.title.get_position()[1] == 1.
            assert ax.get_autoscalex_on()
            assert not ax.get_autoscaley_on()
            assert os.path.isfile(plot.outputfile)

            plot.finalize(close=True)
            assert ax.get_title() == ''  # Axes.cla() called
        finally:
            if os.path.isfile(plot.outputfile):
                os.remove(plot.outputfile)

    def test_apply_parameters(self, plot):
        fig = Plot()
        ax = fig.gca()
        plot.apply_parameters(ax, **{
            'xlim': (10, 20),
            'no-blah': 'anything',
            'grid': False,
        })
        assert ax.get_xlim() == (10, 20)
