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

"""Utilities for data handling and display
"""

from __future__ import division

import operator
import warnings
from collections import OrderedDict

from six import string_types
from six.moves import (reduce, zip_longest)

import numpy

from astropy import units

from gwpy.segments import (DataQualityFlag, SegmentList, Segment)
from gwpy.frequencyseries import FrequencySeries
from gwpy.spectrogram import SpectrogramList

from .. import globalv
from ..utils import (vprint, safe_eval)
from ..channels import get_channel
from .utils import (use_segmentlist, get_fftparams, make_globalv_key)
from .timeseries import (get_timeseries, get_timeseries_dict)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


@use_segmentlist
def get_coherence_spectrogram(channel_pair, segments, config=None,
                              cache=None, query=True, nds=None,
                              return_=True, frametype=None, nproc=1,
                              datafind_error='raise', return_components=False,
                              **fftparams):
    """Retrieve the time-series and generate a coherence spectrogram of
    the two given channels
    """
    specs = _get_coherence_spectrogram(channel_pair, segments,
                                       config=config, cache=cache,
                                       query=query, nds=nds,
                                       return_=return_, frametype=frametype,
                                       nproc=nproc,
                                       datafind_error=datafind_error,
                                       return_components=return_components,
                                       **fftparams)

    return specs


@use_segmentlist
def _get_coherence_spectrogram(channel_pair, segments, config=None,
                               cache=None, query=True, nds=None,
                               return_=True, frametype=None, nproc=1,
                               datafind_error='raise', return_components=False,
                               **fftparams):

    channel1 = get_channel(channel_pair[0])
    channel2 = get_channel(channel_pair[1])

    # clean fftparams dict using channel 1 default values
    fftparams.setdefault('method', 'welch')
    fftparams = get_fftparams(channel1, **fftparams)

    # key used to store the coherence spectrogram in globalv
    key = make_globalv_key(channel_pair, fftparams)

    # keys used to store component spectrograms in globalv
    components = ('Cxy', 'Cxx', 'Cyy')
    ckeys = [
        make_globalv_key([channel1, channel2], fftparams),
        make_globalv_key(channel1, fftparams),
        make_globalv_key(channel2, fftparams),
    ]

    # convert fftparams to regular dict
    fftparams = fftparams.dict()

    # work out what new segments are needed
    # need to truncate to segments of integer numbers of strides
    stride = float(fftparams.pop('stride'))
    overlap = float(fftparams['overlap'])
    new = type(segments)()
    for seg in segments - globalv.SPECTROGRAMS.get(
            key, SpectrogramList()).segments:
        dur = float(abs(seg)) // stride * stride
        if dur < stride + overlap:
            continue
        new.append(type(seg)(seg[0], seg[0]+dur))

    # extract FFT params for TimeSeries.spectrogram
    spec_fftparams = fftparams.copy()
    for fftkey in ('method', 'scheme',):
        fftparams.pop(fftkey, None)

    # if there are no existing spectrogram, initialize as a list
    globalv.SPECTROGRAMS.setdefault(key, SpectrogramList())

    # XXX HACK: use dummy timeseries to find lower sampling rate
    if len(segments) > 0:
        s = segments[0].start
        dts1 = get_timeseries(channel1, SegmentList([Segment(s, s+1)]),
                              config=config, cache=cache, frametype=frametype,
                              nproc=nproc, query=query,
                              datafind_error=datafind_error, nds=nds)
        dts2 = get_timeseries(channel2, SegmentList([Segment(s, s+1)]),
                              config=config, cache=cache, frametype=frametype,
                              nproc=nproc, query=query,
                              datafind_error=datafind_error, nds=nds)
        sampling = min(dts1[0].sample_rate.value, dts2[0].sample_rate.value)
    else:
        sampling = None

    # initialize component lists if they don't exist yet
    for ck in ckeys:
        globalv.COHERENCE_COMPONENTS.setdefault(ck, SpectrogramList())

    # get data if query=True or if there are new segments
    query &= abs(new) != 0

    if query:

        # the intersecting segments will be calculated when needed
        intersection = None

        # loop over components needed to calculate coherence
        for comp in components:

            # key used to store this component in globalv (incl sample rate)
            ckey = ckeys[components.index(comp)]

            try:
                filter_ = channel1.frequency_response
            except AttributeError:
                filter_ = None
            else:
                if isinstance(filter_, string_types):
                    filter_ = safe_eval(filter_, strict=True)

            # check how much of this component still needs to be calculated
            req = new - globalv.COHERENCE_COMPONENTS.get(
                            ckey, SpectrogramList()).segments

            # get data if there are new segments
            if abs(req) != 0:

                # calculate intersection of segments lazily
                # this should occur on first pass (Cxy)
                if intersection is None:
                    total1 = get_timeseries(
                                 channel1, req, config=config,
                                 cache=cache, frametype=frametype,
                                 nproc=nproc, query=query,
                                 datafind_error=datafind_error,
                                 nds=nds)
                    total2 = get_timeseries(
                                 channel2, req, config=config,
                                 cache=cache, frametype=frametype,
                                 nproc=nproc, query=query,
                                 datafind_error=datafind_error,
                                 nds=nds)
                    intersection = total1.segments & total2.segments

                # get required timeseries data (using intersection)
                tslist1, tslist2 = [], []
                if comp in ('Cxy', 'Cxx'):
                    tslist1 = get_timeseries(
                                  channel1, intersection, config=config,
                                  cache=cache, frametype=frametype,
                                  nproc=nproc, query=query,
                                  datafind_error=datafind_error,
                                  nds=nds)
                if comp in ('Cxy', 'Cyy'):
                    tslist2 = get_timeseries(
                                  channel2, intersection, config=config,
                                  cache=cache, frametype=frametype,
                                  nproc=nproc, query=query,
                                  datafind_error=datafind_error,
                                  nds=nds)

                # calculate component
                if len(tslist1) + len(tslist2):
                    vprint("    Calculating component %s for coherence "
                           "spectrogram for %s and %s @ %d Hz" % (
                               comp, str(channel1), str(channel2), sampling))

                for ts1, ts2 in zip_longest(tslist1, tslist2):

                    # ensure there is enough data to do something with
                    if comp in ('Cxx', 'Cxy') and abs(ts1.span) < stride:
                        continue
                    elif comp in ('Cyy', 'Cxy') and abs(ts2.span) < stride:
                        continue

                    # downsample if necessary
                    if ts1 is not None and ts1.sample_rate.value != sampling:
                        ts1 = ts1.resample(sampling)
                    if ts2 is not None and ts2.sample_rate.value != sampling:
                        ts2 = ts2.resample(sampling)

                    # ignore units when calculating coherence
                    if ts1 is not None:
                        ts1._unit = units.Unit('count')
                    if ts2 is not None:
                        ts2._unit = units.Unit('count')

                    # calculate the component spectrogram
                    if comp == 'Cxy':
                        specgram = ts1.csd_spectrogram(
                            ts2, stride, nproc=nproc, **fftparams)
                    elif comp == 'Cxx':
                        specgram = ts1.spectrogram(stride, nproc=nproc,
                                                   **spec_fftparams)
                    elif comp == 'Cyy':
                        specgram = ts2.spectrogram(stride, nproc=nproc,
                                                   **spec_fftparams)

                    if filter_:
                        specgram = (specgram**(1/2.)).filter(*filter_,
                                                             inplace=True) ** 2
                    add_coherence_component_spectrogram(specgram, key=ckey)

                    vprint('.')

                if len(tslist1) + len(tslist2):
                    vprint('\n')

        # calculate coherence from the components and store in globalv
        for seg in new:
            cxy, cxx, cyy = [
                _get_from_list(globalv.COHERENCE_COMPONENTS[ck], seg) for
                ck in ckeys]
            csg = abs(cxy)**2 / cxx / cyy
            globalv.SPECTROGRAMS[key].append(csg)
            globalv.SPECTROGRAMS[key].coalesce()

    if not return_:
        return

    elif return_components:

        # return list of component spectrogram lists
        out = [SpectrogramList(), SpectrogramList(), SpectrogramList()]
        for comp in components:
            index = components.index(comp)
            ckey = ckeys[index]
            for specgram in globalv.COHERENCE_COMPONENTS[ckey]:
                for seg in segments:
                    if abs(seg) < specgram.dt.value:
                        continue
                    if specgram.span.intersects(seg):
                        common = specgram.span & type(seg)(
                                     seg[0], seg[1] + specgram.dt.value)
                        s = specgram.crop(*common)
                        if s.shape[0]:
                            out[index].append(s)
            out[index] = out[index].coalesce()
        return out

    else:

        # return list of coherence spectrograms
        out = SpectrogramList()
        for specgram in globalv.SPECTROGRAMS[key]:
            for seg in segments:
                if abs(seg) < specgram.dt.value:
                    continue
                if specgram.span.intersects(seg):
                    common = specgram.span & type(seg)(
                                 seg[0], seg[1] + specgram.dt.value)
                    s = specgram.crop(*common)
                    if s.shape[0]:
                        out.append(s)
        return out.coalesce()


def get_coherence_spectrum(channel_pair, segments, config=None,
                           cache=None, query=True, nds=None, return_=True,
                           **fftparams):
    """Retrieve the time-series and generate a coherence spectrogram of the given
    channel
    """

    channel1 = get_channel(channel_pair[0])
    channel2 = get_channel(channel_pair[1])

    if isinstance(segments, DataQualityFlag):
        name = ','.join([channel1.ndsname, channel2.ndsname, segments.name])
        segments = segments.active
    else:
        name = channel1.ndsname + ',' + channel2.ndsname
    name += ',%s' % format
    cmin = name + '.min'
    cmax = name + '.max'

    if name not in globalv.COHERENCE_SPECTRUM:
        vprint("    Calculating 5/50/95 percentile spectra for %s"
               % name.rsplit(',', 1)[0])

        # ask for a list of component spectrograms (a list of SpectrogramLists)
        speclist = get_coherence_spectrogram(
                       channel_pair, segments, config=config, cache=cache,
                       query=query, nds=nds, return_components=True,
                       **fftparams)

        cdict = {}
        components = ('Cxy', 'Cxx', 'Cyy')

        # join spectrograms in each list so we can average it
        for comp in components:
            index = components.index(comp)
            try:
                cdict[comp] = speclist[index].join(gap='ignore')
            except ValueError as e:
                if 'units do not match' in str(e):
                    warnings.warn(str(e))
                    for spec in speclist[index][1:]:
                        spec.unit = speclist[0].unit
                    cdict[comp] = speclist[index].join(gap='ignore')
                else:
                    raise

        # average spectrograms to get PSDs and CSD
        try:
            Cxy = complex_percentile(cdict['Cxy'], 50)
            Cxx = cdict['Cxx'].percentile(50)
            Cyy = cdict['Cyy'].percentile(50)
            globalv.COHERENCE_SPECTRUM[name] = FrequencySeries(
                numpy.abs(Cxy)**2 / Cxx / Cyy, f0=Cxx.f0, df=Cxx.df)
        except (ValueError, IndexError):
            globalv.COHERENCE_SPECTRUM[name] = FrequencySeries(
                [], channel=channel1, f0=0, df=1, unit=units.Unit(''))
            globalv.COHERENCE_SPECTRUM[cmin] = globalv.COHERENCE_SPECTRUM[name]
            globalv.COHERENCE_SPECTRUM[cmax] = globalv.COHERENCE_SPECTRUM[name]
        else:
            # FIXME: how to calculate percentiles correctly?
            globalv.COHERENCE_SPECTRUM[cmin] = FrequencySeries(
                abs(complex_percentile(cdict['Cxy'], 5))**2 /
                cdict['Cxx'].percentile(95) / cdict['Cyy'].percentile(95),
                f0=Cxx.f0, df=Cxx.df)
            globalv.COHERENCE_SPECTRUM[cmax] = FrequencySeries(
                abs(complex_percentile(cdict['Cxy'], 95))**2 /
                cdict['Cxx'].percentile(5) / cdict['Cyy'].percentile(5),
                f0=Cxx.f0, df=Cxx.df)

        # set the spectrum's name manually; this will be used for the legend
        globalv.COHERENCE_SPECTRUM[name].name = (
            channel1.ndsname + '\n' + channel2.ndsname)

        vprint(".\n")

    if not return_:
        return

    cmin = '%s.min' % name
    cmax = '%s.max' % name
    out = (globalv.COHERENCE_SPECTRUM[name], globalv.COHERENCE_SPECTRUM[name],
           globalv.COHERENCE_SPECTRUM[name])
    return out


def add_coherence_component_spectrogram(specgram, key=None, coalesce=True):
    """Add a `Coherence spectrogram` to the global memory cache
    """
    if key is None:
        key = specgram.name or str(specgram.channel)
    globalv.COHERENCE_COMPONENTS.setdefault(key, SpectrogramList())
    globalv.COHERENCE_COMPONENTS[key].append(specgram)
    if coalesce:
        globalv.COHERENCE_COMPONENTS[key].coalesce()


@use_segmentlist
def get_coherence_spectrograms(channel_pairs, segments, config=None,
                               cache=None, query=True, nds=None,
                               return_=True, frametype=None, nproc=1,
                               datafind_error='raise', **fftparams):
    """Get coherence spectrograms for multiple channels
    """
    if fftparams.get('method', 'welch') != 'welch':
        raise ValueError("Cannot process coherence data with method=%r"
                         % fftparams.get('method'))
    fftparams['method'] = 'welch'
    channels = list(map(get_channel, channel_pairs))
    pairs = list(zip(channels[0::2], channels[1::2]))

    # get timeseries data in bulk
    if query:
        qchannels = []
        havesegs = []
        for c1, c2 in pairs:
            c1 = get_channel(c1)
            c2 = get_channel(c2)
            fftparams_ = get_fftparams(c1, **fftparams)
            key = make_globalv_key((c1, c2), fftparams_)
            qchannels.extend((c1, c2))
            havesegs.append(globalv.SPECTROGRAMS.get(
                key, SpectrogramList()).segments)
        havesegs = reduce(operator.and_, havesegs)
        new = segments - havesegs
        strides = set([getattr(c, 'stride', 0) for c in qchannels])
        if len(strides) == 1:
            stride = strides.pop()
            new = type(new)([s for s in new if abs(s) >= stride])
        get_timeseries_dict(qchannels, new, config=config, cache=cache,
                            nproc=nproc, frametype=frametype,
                            datafind_error=datafind_error, nds=nds,
                            return_=False)
    # loop over channels and generate spectrograms
    out = OrderedDict()

    for channel_pair in pairs:
        out[channel_pair] = get_coherence_spectrogram(
            channel_pair, segments, config=config, cache=cache, query=query,
            nds=nds, nproc=nproc, return_=return_,
            datafind_error=datafind_error, **fftparams)
    return out


def _get_from_list(serieslist, segment):
    """Internal function to crop a series from a serieslist

    Should only be used in situations where the existence of the target
    data within the list is guaranteed
    """
    for series in serieslist:
        if segment in series.span:
            return series.crop(*segment)
    raise ValueError("Cannot crop series for segment %s from list"
                     % str(segment))


def complex_percentile(array, percentile):
    re = array.real.percentile(percentile)
    im = array.imag.percentile(percentile) * 1j
    return re + im
