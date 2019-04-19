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

"""Get spectrograms and spectra
"""

from __future__ import division

import os.path
import operator
import warnings
from collections import OrderedDict

from six import string_types
from six.moves import reduce

# imports for filter
from math import pi  # noqa: F401

import numpy

from scipy import interpolate

from astropy import units

from gwpy.segments import DataQualityFlag
from gwpy.frequencyseries import FrequencySeries
from gwpy.spectrogram import SpectrogramList

from .. import (globalv, io)
from ..utils import (vprint, safe_eval)
from ..channels import (
    get_channel,
    split_combination as split_channel_combination,
)
from .utils import (use_segmentlist, make_globalv_key, get_fftparams)
from .mathutils import (get_with_math, parse_math_definition)
from .timeseries import (get_timeseries, get_timeseries_dict)

OPERATOR = {
    '*': operator.mul,
    '-': operator.sub,
    '+': operator.add,
    '/': operator.truediv,
}

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


# -- spectrogram --------------------------------------------------------------

@use_segmentlist
def get_spectrogram(channel, segments, config=None, cache=None,
                    query=True, nds=None, format='power', return_=True,
                    frametype=None, nproc=1, datafind_error='raise',
                    **fftparams):
    """Retrieve the time-series and generate a spectrogram of the given
    channel
    """
    channel = get_channel(channel)

    # read data for all sub-channels
    specs = []
    channels = split_channel_combination(channel)
    for c in channels:
        specs.append(_get_spectrogram(c, segments, config=config, cache=cache,
                                      query=query, nds=nds, format=format,
                                      return_=return_, frametype=frametype,
                                      nproc=nproc,
                                      datafind_error=datafind_error,
                                      **fftparams))
    if return_ and len(channels) == 1:
        return specs[0]
    elif return_:
        return get_with_math(
            channel, segments, _get_spectrogram, _get_spectrogram,
            config=config, query=False, format=format, return_=True,
            **fftparams)


@use_segmentlist
def _get_spectrogram(channel, segments, config=None, cache=None,
                     query=True, nds=None, format='power', return_=True,
                     frametype=None, nproc=1,
                     datafind_error='raise', **fftparams):
    channel = get_channel(channel)

    # if we aren't given a method, check to see whether data have already
    # been processed, if so, choose that one
    if fftparams.get('method', None) is None:
        methods = set([key.split(';')[1] for key in globalv.SPECTROGRAMS
                       if key.startswith('%s;' % channel.ndsname)])
        try:
            fftparams['method'] = list(methods)[0]
        except IndexError:
            fftparams['method'] = 'welch'

    # clean fftparams dict using channel default values
    fftparams = get_fftparams(channel, **fftparams)
    # override special-case methods
    if format in ['rayleigh']:
        fftparams.method = format

    # key used to store the coherence spectrogram in globalv
    key = make_globalv_key(channel, fftparams)

    # keep FftParams as a dict for convenience
    fftparams = fftparams.dict()

    # read segments from global memory
    havesegs = globalv.SPECTROGRAMS.get(key, SpectrogramList()).segments
    new = segments - havesegs
    query &= abs(new) != 0

    globalv.SPECTROGRAMS.setdefault(key, SpectrogramList())

    if query:
        # extract spectrogram stride from dict
        try:
            stride = float(fftparams.pop('stride'))
        except (TypeError, KeyError) as e:
            msg = ('cannot parse a spectrogram stride from the kwargs '
                   'given, please give some or all of fftlength, overlap, '
                   'stride')
            if isinstance(e, TypeError):
                e.args = (msg,)
                raise
            raise TypeError(msg)

        # read channel information
        try:
            filter_ = channel.frequency_response
        except AttributeError:
            filter_ = None
        else:
            if isinstance(filter_, string_types) and os.path.isfile(filter_):
                filter_ = io.read_frequencyseries(filter_)
            elif isinstance(filter_, string_types):
                filter_ = safe_eval(filter_, strict=True)

        # get time-series data
        timeserieslist = get_timeseries(channel, new, config=config,
                                        cache=cache, frametype=frametype,
                                        nproc=nproc, query=query,
                                        datafind_error=datafind_error, nds=nds)
        # calculate spectrograms
        if len(timeserieslist):
            vprint("    Calculating (%s) spectrograms for %s"
                   % (fftparams['method'], str(channel)))
        for ts in timeserieslist:
            # if too short for a single segment, continue
            if abs(ts.span) < (stride + fftparams.get('overlap', 0)):
                continue
            # truncate timeseries to integer number of strides
            d = size_for_spectrogram(ts.duration.to('s').value, stride,
                                     fftparams['fftlength'],
                                     fftparams.get('overlap', 0))
            ts = ts.crop(ts.span[0], ts.span[0] + d, copy=False)
            # calculate spectrogram
            try:
                # rayleigh spectrogram has its own instance method
                if fftparams.get('method', None) == 'rayleigh':
                    spec_kw = fftparams.copy()
                    for fftkey in ('method', 'scheme',):  # remove ASD keys
                        spec_kw.pop(fftkey, None)
                    spec_func = ts.rayleigh_spectrogram
                else:
                    spec_kw = fftparams
                    spec_func = ts.spectrogram
                specgram = spec_func(stride, nproc=nproc, **spec_kw)
            except ZeroDivisionError:
                if stride == 0:
                    raise ZeroDivisionError("Spectrogram stride is 0")
                elif fftparams['fftlength'] == 0:
                    raise ZeroDivisionError("FFT length is 0")
                else:
                    raise
            except ValueError as e:
                if 'has no unit' in str(e):
                    unit = ts.unit
                    ts._unit = units.Unit('count')
                    specgram = ts.spectrogram(stride, nproc=nproc, **fftparams)
                    specgram._unit = unit ** 2 / units.Hertz
                else:
                    raise
            if isinstance(filter_, FrequencySeries) and (
                    fftparams['method'] not in ['rayleigh']):
                specgram = apply_transfer_function_series(specgram, filter_)
            elif filter_ and fftparams['method'] not in ['rayleigh']:
                # manually setting x0 is a hack against precision error
                # somewhere inside the **(1/2.) operation (Quantity)
                x0 = specgram.x0.value
                specgram = (specgram ** (1/2.)).filter(*filter_,
                                                       inplace=True) ** 2
                specgram.x0 = x0
            if specgram.unit is None:
                specgram._unit = channel.unit
            elif len(globalv.SPECTROGRAMS[key]):
                specgram._unit = globalv.SPECTROGRAMS[key][-1].unit
            add_spectrogram(specgram, key=key)
            vprint('.')
        if len(timeserieslist):
            vprint('\n')

    if not return_:
        return

    # return correct data
    out = SpectrogramList()
    for specgram in globalv.SPECTROGRAMS[key]:
        for seg in segments:
            if abs(seg) < specgram.dt.value:
                continue
            if specgram.span.intersects(seg):
                common = specgram.span & type(seg)(seg[0],
                                                   seg[1] + specgram.dt.value)
                s = specgram.crop(*common)
                if format in ['amplitude', 'asd']:
                    s = s**(1/2.)
                elif format in ['rayleigh']:
                    # XXX FIXME: this corrects the bias offset in Rayleigh
                    med = numpy.median(s.value)
                    s /= med
                if s.shape[0]:
                    out.append(s)
    return out.coalesce()


def add_spectrogram(specgram, key=None, coalesce=True):
    """Add a `Spectrogram` to the global memory cache
    """
    if key is None:
        key = specgram.name or str(specgram.channel)
    globalv.SPECTROGRAMS.setdefault(key, SpectrogramList())
    globalv.SPECTROGRAMS[key].append(specgram)
    if coalesce:
        globalv.SPECTROGRAMS[key].coalesce()


@use_segmentlist
def get_spectrograms(channels, segments, config=None, cache=None, query=True,
                     nds=None, format='power', return_=True, frametype=None,
                     nproc=1, datafind_error='raise', **fftparams):
    """Get spectrograms for multiple channels
    """
    channels = map(get_channel, channels)

    # get timeseries data in bulk
    if query:
        # get underlying list of data channels to read
        qchannels = map(get_channel,
                        set([c for group in
                             map(split_channel_combination, channels)
                             for c in group]))

        # work out FFT params and storage keys for each data channel
        keys = []
        for channel in qchannels:
            fftparams_ = get_fftparams(channel, **fftparams)
            keys.append(make_globalv_key(channel, fftparams_))

        # restrict segments to those big enough to hold >= 1 stride
        strides = set([getattr(c, 'stride', 0) for c in qchannels])
        if len(strides) == 1:
            stride = strides.pop()
            segments = type(segments)(s for s in segments if abs(s) >= stride)

        # work out new segments for which to read data
        havesegs = reduce(operator.and_, (globalv.SPECTROGRAMS.get(
            key, SpectrogramList()).segments for key in keys))
        new = segments - havesegs

        # read data for new segments
        get_timeseries_dict(qchannels, new, config=config, cache=cache,
                            nproc=nproc, frametype=frametype,
                            datafind_error=datafind_error, nds=nds,
                            return_=False)
    # loop over channels and generate spectrograms
    out = OrderedDict()
    for channel in channels:
        out[channel] = get_spectrogram(
            channel, segments, config=config, cache=cache, query=query,
            nds=nds, format=format, nproc=nproc,
            return_=return_, datafind_error=datafind_error,
            **fftparams)
    return out


def size_for_spectrogram(size, stride, fftlength, overlap):
    if size < stride:
        return None
    x = size // stride * stride + overlap
    if x > size:
        x -= fftlength
    return x


def apply_transfer_function_series(specgram, tfunc):
    """Multiply a spectrogram by a transfer function `FrequencySeries`

    This method interpolates the transfer function onto the frequency vector
    of the spectrogram, so should work regardless of the inputs
    """
    # interpolate transfer function onto spectrogram frequency series
    interpolator = interpolate.interp1d(tfunc.frequencies.value, tfunc.value)
    itfunc = numpy.zeros((1, specgram.frequencies.size))
    known = specgram.frequencies.value >= tfunc.frequencies.value[0]
    known &= specgram.frequencies.value <= tfunc.frequencies.value[-1]
    itfunc[0, :][known] = interpolator(specgram.frequencies.value[known])
    # and multiply
    return (specgram ** (1/2.) * itfunc) ** 2


# -- spectrum -----------------------------------------------------------------

@use_segmentlist
def get_spectrum(channel, segments, config=None, cache=None,
                 query=True, nds=None, format='power', return_=True,
                 frametype=None, nproc=1, datafind_error='raise',
                 **fftparams):
    """Retrieve the time-series and generate a spectrogram of the given
    channel
    """
    channel = get_channel(channel)
    if isinstance(segments, DataQualityFlag):
        name = ','.join([channel.ndsname, segments.name])
        segments = segments.active
    else:
        name = channel.ndsname
    name += ',%s' % format

    # read data for all sub-channels
    specs = []
    channels = list(parse_math_definition(str(channel))[0])
    if len(channels) == 0:
        channels = [channel]
    for c in channels:
        specs.append(_get_spectrum(c, segments, config=config, cache=cache,
                                   query=query, nds=nds, format=format,
                                   return_=return_, frametype=frametype,
                                   nproc=nproc,
                                   datafind_error=datafind_error,
                                   **fftparams))
    if return_ and len(channels) == 1:
        return specs[0]
    elif return_:
        return [get_with_math(
                    channel, segments, _get_spectrum, _get_spectrum,
                    config=config, format=format, return_=True,
                    which=which)[0] for which in ['mean', 'min', 'max']]


def _get_spectrum(channel, segments, config=None, cache=None, query=True,
                  nds=None, format='power', return_=True, which='all',
                  **fftparams):
    """Retrieve the time-series and generate a spectrogram of the given
    channel
    """
    channel = get_channel(channel)
    if isinstance(segments, DataQualityFlag):
        name = ','.join([channel.ndsname, segments.name])
        segments = segments.active
    else:
        name = channel.ndsname
    name += ',%s' % format
    cmin = '%s.min' % name
    cmax = '%s.max' % name

    if name not in globalv.SPECTRUM:
        if os.path.isfile(channel.ndsname):
            globalv.SPECTRUM[name] = io.read_frequencyseries(channel.ndsname)
            globalv.SPECTRUM[cmin] = globalv.SPECTRUM[name]
            globalv.SPECTRUM[cmax] = globalv.SPECTRUM[name]
        else:
            fftparams.setdefault('fftlength', 1)
            fftparams.setdefault('overlap', 0.5)
            if 'stride' not in fftparams and 'fftlength' in fftparams:
                fftparams.setdefault('stride', fftparams['fftlength'])

            speclist = get_spectrogram(channel, segments, config=config,
                                       cache=cache, query=query, nds=nds,
                                       format=format, **fftparams)
            try:
                specgram = speclist.join(gap='ignore')
            except ValueError as e:
                if 'units do not match' in str(e):
                    warnings.warn(str(e))
                    for spec in speclist[1:]:
                        spec.unit = speclist[0].unit
                    specgram = speclist.join(gap='ignore')
                else:
                    raise
            try:
                globalv.SPECTRUM[name] = specgram.percentile(50)
            except (ValueError, IndexError):
                globalv.SPECTRUM[name] = FrequencySeries(
                    [], channel=channel, f0=0, df=1, unit=units.Unit(''))
                globalv.SPECTRUM[cmin] = globalv.SPECTRUM[name]
                globalv.SPECTRUM[cmax] = globalv.SPECTRUM[name]
            else:
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
