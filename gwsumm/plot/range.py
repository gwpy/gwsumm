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

"""Definitions for range plots
"""

from __future__ import division

from .. import (globalv, version)
from .registry import (get_plot, register_plot)
from ..data import get_range

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


class RangePlotMixin(object):
    data = 'spectrogram'
    defaults = {
        'snr': 8.0,
        'stride': 60.,
        'fftlength': 8,
        'overlap': 4,
        'fmin': 10,
    }

    def __init__(self, *args, **kwargs):
        super(RangePlotMixin, self).__init__(*args, **kwargs)
        self.rangeparams = {}
        for key in ['mass1', 'mass2', 'snr', 'energy', 'stride', 'fftlength',
                    'overlap', 'fmin', 'fmax']:
            try:
                value = self.pargs.pop(key)
            except KeyError:
                continue
            if not isinstance(value, (tuple, list)):
                value = [value]*len(self.channels)
            self.rangeparams[key] = value

    def process(self):
        """Read in all necessary data, and generate the figure.
        """
        # generate data
        keys = []
        for i, channel in enumerate(self.channels):
            kwargs = dict((key, self.rangeparams[key][i]) for
                          key in self.rangeparams if
                          self.rangeparams[key][i] is not None)
            if self.state and not self.all_data:
                valid = self.state.active
            else:
                valid = SegmentList([self.span])
            rlist = get_range(channel, valid, **kwargs)
            try:
                keys.append(rlist[0].channel)
            except IndexError:
                keys.append(channel)

        # reset channel lists and generate time-series plot
        channels = self.channels
        outputfile = self.outputfile
        self.channels = keys
        out = super(RangePlotMixin, self).process(outputfile=outputfile)
        self.channels = channels
        return out


class RangeDataPlot(RangePlotMixin, get_plot('timeseries')):
    type = 'range'
    _threadsafe = False
    defaults = get_plot('timeseries').defaults.copy()
    defaults.update(RangePlotMixin.defaults.copy())
    defaults.update({
        'ylabel': 'Sensitive distance [Mpc]',
    })

register_plot(RangeDataPlot)


class RangeDataHistogramPlot(RangePlotMixin, get_plot('histogram')):
    type = 'range-histogram'
    _threadsafe = False
    defaults = get_plot('histogram').defaults.copy()
    defaults.update(RangePlotMixin.defaults.copy())
    defaults.update({
        'xlabel': 'Sensitive distance [Mpc]',
    })

register_plot(RangeDataHistogramPlot)
