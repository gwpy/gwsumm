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

"""Get range data
"""

import re

import numpy

from gwpy import astro
from gwpy.timeseries import (TimeSeries, TimeSeriesList)
from gwpy.frequencyseries import FrequencySeries

from .. import globalv
from ..utils import re_cchar
from ..channels import get_channel
from .utils import (use_segmentlist, make_globalv_key)
from .timeseries import (add_timeseries, get_timeseries)
from .spectral import get_spectrogram

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


def get_range_channel(channel, **rangekwargs):
    """Return the meta-channel name used to store range data
    """
    if not rangekwargs:
        rangekwargs = {'mass1': 1.4, 'mass2': 1.4}
    re_float = re.compile(r'[.-]')
    rkey = '_'.join(['%s_%s' % (re_cchar.sub('_', key),
                                re_float.sub('_', str(val))) for key, val in
                    rangekwargs.items()])
    channel = get_channel(channel)
    return '%s_%s' % (channel.ndsname, rkey)


@use_segmentlist
def get_range(channel, segments, config=None, cache=None,
              query=True, nds=None, return_=True, nproc=1,
              datafind_error='raise', frametype=None,
              stride=None, fftlength=None, overlap=None,
              method=None, **rangekwargs):
    """Calculate the sensitive distance for a given strain channel
    """
    if not rangekwargs:
        rangekwargs = {'mass1': 1.4, 'mass2': 1.4}
    if 'energy' in rangekwargs:
        range_func = astro.burst_range
    else:
        range_func = astro.inspiral_range
    channel = get_channel(channel)
    key = make_globalv_key(get_range_channel(channel, **rangekwargs))
    # get old segments
    havesegs = globalv.DATA.get(key, TimeSeriesList()).segments
    new = segments - havesegs
    query &= abs(new) != 0
    # calculate new range
    if query:
        # get spectrograms
        spectrograms = get_spectrogram(channel, new, config=config,
                                       cache=cache, nproc=nproc,
                                       frametype=frametype, format='psd',
                                       datafind_error=datafind_error, nds=nds,
                                       stride=stride, fftlength=fftlength,
                                       overlap=overlap, method=method)
        # calculate range for each PSD in each spectrogram
        for sg in spectrograms:
            ts = TimeSeries(numpy.zeros(sg.shape[0],), unit='Mpc',
                            epoch=sg.epoch, dx=sg.dx, channel=key)
            for i in range(sg.shape[0]):
                psd = sg[i]
                psd = FrequencySeries(psd.value, f0=psd.x0, df=psd.dx)
                ts[i] = range_func(psd, **rangekwargs)
            add_timeseries(ts, key=key)

    if return_:
        return get_timeseries(key, segments, query=False)
