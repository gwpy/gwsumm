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

from gwpy import astro
from gwpy.frequencyseries import FrequencySeries
from gwpy.spectrogram import SpectrogramList
from gwpy.timeseries import TimeSeriesList

from .. import globalv
from ..utils import re_cchar
from ..channels import get_channel
from .utils import (use_segmentlist, make_globalv_key)
from .timeseries import (add_timeseries, get_timeseries)
from .spectral import (add_spectrogram, get_spectrogram)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__credits__ = 'Alex Urban <alexander.urban@ligo.org>'


def _metadata(channel, **rangekwargs):
    """Return a common set of metadata for range calculations
    """
    channel = get_channel(channel)
    key = make_globalv_key(get_range_channel(channel, **rangekwargs))
    return (channel, key)


def _segments_diff(segments, havesegs, query=True):
    """Return a diff of two `SegmentList`, and determine whether to query for
    new data
    """
    new = segments - havesegs
    query &= abs(new) != 0
    return (new, query)


def get_range_channel(channel, **rangekwargs):
    """Return the meta-channel name used to store range data
    """
    channel = get_channel(channel)
    if rangekwargs:
        re_float = re.compile(r'[.-]')
        rkey = '_'.join(['%s_%s' % (re_cchar.sub('_', key),
                                    re_float.sub('_', str(val))) for key, val
                        in rangekwargs.items()])
        return '%s_%s' % (channel.ndsname, rkey)
    return channel.ndsname


@use_segmentlist
def get_range(channel, segments, config=None, cache=None,
              query=True, nds=None, return_=True, nproc=1,
              datafind_error='raise', frametype=None,
              stride=None, fftlength=None, overlap=None,
              method=None, **rangekwargs):
    """Calculate the sensitive distance for a given strain channel
    """
    channel, key = _metadata(channel, **rangekwargs)
    # get new segments
    havesegs = globalv.DATA.get(key, TimeSeriesList()).segments
    new, query = _segments_diff(segments, havesegs, query)
    if query:  # calculate new range
        spectrograms = get_spectrogram(
            channel, new, config=config, cache=cache, nds=nds, format='psd',
            frametype=frametype, nproc=nproc, datafind_error=datafind_error,
            stride=stride, fftlength=fftlength, overlap=overlap, method=method)
        for sg in spectrograms:  # calculate range for each spectrogram
            ts = astro.range_timeseries(sg, **rangekwargs)
            ts.channel = key
            add_timeseries(ts, key=key)

    if return_:
        return get_timeseries(key, segments, query=False)


@use_segmentlist
def get_range_spectrogram(channel, segments, config=None, cache=None,
                          query=True, nds=None, return_=True, nproc=1,
                          datafind_error='raise', frametype=None, stride=60,
                          fftlength=None, overlap=None, method=None,
                          **rangekwargs):
    """Estimate the spectral contribution to sensitive distance for a given
    strain channel
    """
    channel, key = _metadata(channel, **rangekwargs)
    # get new segments
    havesegs = globalv.SPECTROGRAMS.get(key, SpectrogramList()).segments
    new, query = _segments_diff(segments, havesegs, query)
    if query:  # calculate new data
        spectrograms = get_spectrogram(
            channel, new, config=config, cache=cache, nds=nds, format='psd',
            frametype=frametype, nproc=nproc, datafind_error=datafind_error,
            stride=stride, fftlength=fftlength, overlap=overlap, method=method)
        for sg in spectrograms:  # get contribution from each spectrogram
            outspec = astro.range_spectrogram(sg, **rangekwargs)
            outspec.channel = key
            add_spectrogram(outspec if 'energy' in rangekwargs else
                            outspec**(1/2.), key=key)

    if return_:
        return globalv.SPECTROGRAMS.get(key, SpectrogramList())


@use_segmentlist
def get_range_spectrum(channel, segments, config=None, cache=None, query=True,
                       nds=None, return_=True, nproc=1, datafind_error='raise',
                       frametype=None, stride=60, fftlength=None, overlap=None,
                       method=None, which='all', **rangekwargs):
    """Compute percentile spectra of the range integrand from a set of
    spectrograms
    """
    name = str(channel)
    cmin = '%s.min' % name
    cmax = '%s.max' % name

    if name not in globalv.SPECTRUM:
        speclist = get_range_spectrogram(
            channel, segments, config=config, cache=cache, query=query,
            nds=nds, return_=return_, nproc=nproc, frametype=frametype,
            datafind_error=datafind_error, method=method, stride=stride,
            fftlength=fftlength, overlap=overlap, **rangekwargs)
        specgram = speclist.join(gap='ignore')
        try:  # store median spectrum
            globalv.SPECTRUM[name] = specgram.percentile(50)
        except (ValueError, IndexError):
            unit = 'Mpc' if 'energy' in rangekwargs else 'Mpc^2 / Hz'
            globalv.SPECTRUM[name] = FrequencySeries(
                [], channel=channel, f0=0, df=1, unit=unit)
            globalv.SPECTRUM[cmin] = globalv.SPECTRUM[name]
            globalv.SPECTRUM[cmax] = globalv.SPECTRUM[name]
        else:  # store percentiles
            globalv.SPECTRUM[cmin] = specgram.percentile(5)
            globalv.SPECTRUM[cmax] = specgram.percentile(95)

    if not return_:
        return

    if which == 'all':
        return (globalv.SPECTRUM[name], globalv.SPECTRUM[cmin],
                globalv.SPECTRUM[cmax])
    if which == 'mean':
        return globalv.SPECTRUM[name]
    if which == 'min':
        return globalv.SPECTRUM[cmin]
    if which == 'max':
        return globalv.SPECTRUM[cmax]
    raise ValueError("Unrecognised value for `which`: %r" % which)
