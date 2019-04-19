# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2017)
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

"""Extensions to the spectrum plot for noise budgets
"""

from __future__ import division

import re

import numpy

from gwpy.segments import SegmentList

from ..data import get_spectrum
from .registry import (get_plot, register_plot)
from .utils import usetex_tex

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


class NoiseBudgetPlot(get_plot('spectrum')):
    """Plot of a noise budget ASD
    """
    type = 'noise-budget'
    data = 'spectrum'
    defaults = get_plot('spectrum').defaults.copy()
    defaults.update({
        'xscale': 'log',
        'yscale': 'log',
        'format': 'asd',
        'sum-label': 'Sum of noises',
        'sum-linestyle': '--',
        'sum-color': 'black',
        'residual-label': 'Residual',
        'residual-linestyle': ':',
        'residual-color': 'grey',
    })

    def _parse_extra_params(self, prefix, **defaults):
        """Parse parameters for an extra plot element
        """
        re_prefix = re.compile(r'\A%s[-_]' % prefix.rstrip('-_'))
        extras = defaults.copy()
        for key in list(self.pargs):
            m = re_prefix.match(key)
            if m:
                extras[key[m.span()[1]:]] = self.pargs.pop(key)
        return extras

    def parse_sum_params(self, **defaults):
        return self._parse_extra_params('sum', **defaults)

    def parse_residual_params(self, **defaults):
        return self._parse_extra_params('residual', **defaults)

    def _draw(self):
        """Load all data, and generate this `SpectrumDataPlot`
        """
        plot = self.init_plot()
        ax = plot.gca()
        ax.grid(b=True, axis='both', which='both')

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        usetex_tex(str(self.state))))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.993, va='top')

        # get spectrum format: 'amplitude' or 'power'
        sdform = self.pargs.pop('format')

        # parse plotting arguments
        plotargs = self.parse_plot_kwargs()
        legendargs = self.parse_legend_kwargs()

        # add data
        sumdata = []
        for i, (channel, pargs) in enumerate(zip(self.channels, plotargs)):
            if self.state and not self.all_data:
                valid = self.state
            else:
                valid = SegmentList([self.span])

            data = get_spectrum(str(channel), valid, query=False,
                                format=sdform, method=None)[0]
            if i and data.size:
                sumdata.append(data)
            else:
                darmdata = data

            # anticipate log problems
            if self.logx:
                data = data[1:]
            if self.logy:
                data.value[data.value == 0] = 1e-100

            pargs.setdefault('zorder', -i)
            ax.plot(data, **pargs)

        # assert all noise terms have the same resolution
        if any([x.dx != sumdata[0].dx for x in sumdata]):
            raise RuntimeError("Noise components have different resolutions, "
                               "cannot construct sum of noises")
        # reshape noises if required
        n = max(x.size for x in sumdata)
        for i, d in enumerate(sumdata):
            if d.size < n:
                sumdata[i] = numpy.resize(
                    numpy.require(d, requirements=['O']), (n,))

        # plot sum of noises
        sumargs = self.parse_sum_params()
        sum_ = sumdata[0] ** 2
        for d in sumdata[1:]:
            sum_ += d ** 2
        ax.plot(sum_ ** (1/2.), zorder=1, **sumargs)
        ax.lines.insert(1, ax.lines.pop(-1))

        # plot residual of noises
        if not self.pargs.pop('no-residual', False):
            resargs = self.parse_residual_params()
            try:
                residual = (darmdata ** 2 - sum_) ** (1/2.)
            except ValueError:
                if not darmdata.size:  # if no data, just copy nothing
                    residual = darmdata
                else:  # other error
                    raise
            ax.plot(residual, zorder=-1000, **resargs)
            ax.lines.insert(1, ax.lines.pop(-1))

        # finalize
        self.apply_parameters(ax, **self.pargs)
        ax.legend(**legendargs)

        return self.finalize()


register_plot(NoiseBudgetPlot)


class RelativeNoiseBudgetPlot(get_plot('spectrum')):
    """Spectrum plot for a `SummaryTab`
    """
    type = 'noise-budget-ratio'
    data = 'spectrum'
    defaults = get_plot('spectrum').defaults.copy()
    defaults.update({
        'xscale': 'log',
        'yscale': 'log',
        'format': 'asd',
    })

    def _draw(self):
        """Load all data, and generate this `SpectrumDataPlot`
        """
        plot = self.init_plot()
        ax = plot.gca()
        ax.grid(b=True, axis='both', which='both')

        if self.state:
            self.pargs.setdefault(
                'suptitle',
                '[%s-%s, state: %s]' % (self.span[0], self.span[1],
                                        usetex_tex(str(self.state))))
        suptitle = self.pargs.pop('suptitle', None)
        if suptitle:
            plot.suptitle(suptitle, y=0.993, va='top')

        # get spectrum format: 'amplitude' or 'power'
        sdform = self.pargs.pop('format')

        # parse plotting arguments
        plotargs = self.parse_plot_kwargs()[0]

        # add data
        sumdata = []
        for i, channel in enumerate(self.channels):
            if self.state and not self.all_data:
                valid = self.state
            else:
                valid = SegmentList([self.span])

            data = get_spectrum(str(channel), valid, query=False,
                                format=sdform, method=None)[0]
            if i and data.size:
                sumdata.append(data)
            else:
                target = data

        if target.size:
            # assert all noise terms have the same resolution
            if any([x.dx != target.dx for x in sumdata]):
                raise RuntimeError("Noise components have different "
                                   "resolutions, cannot construct "
                                   "sum of noises")
            # reshape noises if required
            n = target.size
            for i, d in enumerate(sumdata):
                if d.size < n:
                    sumdata[i] = numpy.resize(
                        numpy.require(d, requirements=['O']), (n,))

            # calculate sum of noises
            sum_ = sumdata[0] ** 2
            for d in sumdata[1:]:
                sum_ += d ** 2
            sum_ **= (1/2.)

            relative = sum_ / target
        else:  # no data, so just use anything as a proxy
            relative = target

        # plot ratio of h(t) to sum of noises
        ax.plot(relative, **plotargs)

        # finalize plot
        self.apply_parameters(ax, **self.pargs)

        return self.finalize()


register_plot(RelativeNoiseBudgetPlot)
