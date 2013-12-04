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

import os
from math import (floor, ceil)
try:
    from configparser import (ConfigParser, NoSectionError, NoOptionError)
except ImportError:
    from ConfigParser import (ConfigParser, NoSectionError, NoOptionError)

import numpy
import nds2

from glue import datafind
from glue.lal import Cache

from gwpy.detector import Channel
from gwpy.segments import (DataQualityFlag, SegmentList)
from gwpy.timeseries import (TimeSeries, TimeSeriesList)
from gwpy.spectrum import Spectrum
from gwpy.spectrogram import SpectrogramList
from gwpy.io import nds as ndsio

from . import (globalv, version)
from .utils import *

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


def get_channel(channel):
    """Define a new :class:`~gwpy.detector.channel.Channel`

    Parameters
    ----------
    channel : `str`
        name of new channel

    Returns
    -------
    Channel : :class:`~gwpy.detector.channel.Channel`
        new channel.
    """
    if ',' in str(channel):
        name, type_ = str(channel).rsplit(',', 1)
        if type_.isdigit():
            type_ = int(type_)
        else:
            type_ = ndsio.NDS2_CHANNEL_TYPE[type_]
        found = globalv.CHANNELS.sieve(name=name, type=type_)
    else:
        name = str(channel)
        found = globalv.CHANNELS.sieve(name=str(channel))
    if len(found) == 1:
        return found[0]
    elif len(found) > 1:
        raise ValueError("Ambiguous channel request '%s', multiple existing "
                         "channels recovered" % str(channel))
    else:
        try:
            new = Channel.query(name)
        except ValueError:
            new = Channel(channel)
        else:
            new.name = str(channel)
        globalv.CHANNELS.append(new)
        return get_channel(str(channel))


def find_frames(ifo, frametype, gpsstart, gpsend, config=ConfigParser(),
                urltype='file', gaps='warn'):
    """Query the datafind server for GWF files for the given type
    """
    # find datafind host:port
    try:
        host = config.get('datafind', 'server')
    except (NoOptionError, NoSectionError):
        try:
            host = os.environ['LIGO_DATAFIND_SERVER']
        except KeyError:
            host = None
            port = None
        else:
            host, port = host.split(':')
            port = int(port)
    else:
        port = config.getint('datafind', 'port')
    # check certificates
    if not port == 80:
        cert, key = datafind.find_credential()
    else:
        cert, key = None, None
    # connect with security
    if cert and key:
        dfconn = datafind.GWDataFindHTTPSConnection(host=host, port=port,
                                                    cert_file=cert,
                                                    key_file=key)
    else:
        dfconn = datafind.GWDataFindHTTPConnection(host=host, port=port)
    # query frames
    ifo = ifo[0].upper()
    gpsstart = int(floor(gpsstart))
    gpsend = int(ceil(gpsend))
    return dfconn.find_frame_urls(ifo[0].upper(), frametype, gpsstart, gpsend,
                                  urltype=urltype, on_gaps=gaps)


def find_cache_segments(*caches):
    """Construct a :class:`~gwpy.segments.segments.SegmentList` describing
    the validity of a given :class:`~glue.lal.Cache`, or list of them.

    Parameters
    ----------
    cache : :class:`~glue.lal.Cache`
        Cache of frame files to check

    Returns
    -------
    segments : :class:`~gwpy.segments.segments.SegmentList`
        list of segments containing in cache
    """
    out = SegmentList()
    nframes = sum(len(c) for c in caches)
    if nframes == 0:
        return out
    for cache in caches:
        # build segment for this cache
        if not len(cache):
            continue
        seg = cache[0].segment
        for e in cache:
            # if new segment doesn't overlap, append and start again
            if e.segment.disjoint(seg):
                out.append(seg)
                seg = e.segment
            # otherwise, append to current segment
            else:
                seg |= e.segment
    # append final segment and return
    out.append(seg)
    return out


def get_timeseries(channel, segments, config=ConfigParser(), cache=Cache(),
                   query=True, nds='guess'):
    """Retrieve the data (time-series) for a given channel
    """
    if isinstance(segments, DataQualityFlag):
        segments = segments.active
    channel = get_channel(channel)
    name = str(channel)

    # read segments from global memory
    havesegs = globalv.DATA.get(str(channel), TimeSeriesList()).segments
    new = segments - havesegs

    # read channel information
    try:
        filter_ = channel.filter
    except AttributeError:
        filter_ = None
    try:
        resample = float(channel.resample)
    except AttributeError:
        resample = None
    resample = None

    # work out whether to use NDS or not
    if nds == 'guess':
        nds = 'LIGO_DATAFIND_SERVER' not in os.environ

    # read new data
    globalv.DATA.setdefault(name, TimeSeriesList())
    query &= (abs(new) > 0)
    if query:
        # open NDS connection
        if nds and config.has_option('nds', 'host'):
            host = config.get('nds', 'host')
            port = config.getint('nds', 'port')
            try:
                ndsconnection = nds2.connection(host, port)
            except RuntimeError as e:
                if 'SASL authentication' in str(e):
                    from gwpy.io.nds import kinit
                    kinit()
                    ndsconnection = nds2.connection(host, port)
            source = 'nds'
        elif nds:
            ndsconnection = None
            source = 'nds'
        # or find frame type and check cache
        else:
            try:
                ftype = channel.frametype
            except AttributeError:
                try:
                    ndstype = channel.type
                except AttributeError:
                    ndstype = nds2.channel.CHANNEL_TYPE_RAW
                if ndstype == nds2.channel.CHANNEL_TYPE_MTREND:
                    new = type(new)([s for s in new if abs(s) >= 60.])
                    ftype = 'M'
                elif ndstype == nds2.channel.CHANNEL_TYPE_STREND:
                    new = type(new)([s for s in new if abs(s) >= 1.])
                    ftype = 'T'
                elif ndstype == nds2.channel.CHANNEL_TYPE_RDS:
                    ftype = 'LDAS_C02_L2'
                elif ndstype == nds2.channel.CHANNEL_TYPE_ONLINE:
                    ftype = 'lldetchar'
                else:
                    ftype = 'R'
                ftype = '%s_%s' % (channel.ifo, ftype)
                # XXX: remove me when L1 moves frame type
                if channel.ifo[0] == 'L' and len(ftype) == 4:
                    ftype = ftype[-1]
            if cache is not None:
                fcache = cache.sieve(description=ftype, exact_match=True)
            if cache is None or len(fcache) == 0 and len(new):
                span = new.extent()
                fcache = find_frames(channel.ifo, ftype, span[0], span[1],
                                     config=config)
            # parse discontiguous cache blocks and rebuild segment list
            cachesegments = find_cache_segments(fcache)
            new &= cachesegments
            source = 'frames'

        # loop through segments, recording data for each
        if len(new):
            vprint("    Fetching data (from %s) for %s"
                   % (source, name))
        for segment in new:
            if nds:
                data = TimeSeries.fetch(channel, segment[0], segment[1],
                                        connection=ndsconnection,
                                        ndschanneltype=channel.type)
            else:
                segcache = fcache.sieve(segment=segment)
                data = TimeSeries.read(segcache, channel, float(segment[0]),
                                       float(segment[1]),
                                       verbose=globalv.VERBOSE)
            data.channel = channel
            if not channel.sample_rate:
                channel.sample_rate = data.sample_rate
            if filter_:
                data = data.filter(*filter_)
            if resample:
                factor = data.sample_rate.value / resample
                if numpy.isclose(factor, int(factor)):
                    data = data.decimate(factor)
                else:
                    data = data.resample(resample)
            globalv.DATA[name].append(data)
            globalv.DATA[name].coalesce()
            vprint(".")
        if len(new):
            vprint("\n")

    # return correct data
    out = TimeSeriesList()
    for ts in globalv.DATA[name]:
        for seg in segments:
            if abs(seg) == 0:
                continue
            if ts.span.intersects(seg):
                cropped = ts.crop(*seg)
                if cropped.size:
                    out.append(cropped)
    return out.coalesce()


def get_spectrogram(channel, segments, config=ConfigParser(), cache=None,
                    query=True, nds='guess', format='power', **fftparams):
    """Retrieve the time-series and generate a spectrogram of the given
    channel
    """
    if isinstance(segments, DataQualityFlag):
        segments = segments.active
    # read segments from global memory
    havesegs = globalv.SPECTROGRAMS.get(str(channel),
                                        SpectrogramList()).segments
    new = segments - havesegs

    globalv.SPECTROGRAMS.setdefault(str(channel), SpectrogramList())

    query &= abs(new) != 0
    if query:
        # read channel information
        try:
            filter_ = channel.frequency_response
        except AttributeError:
            filter_ = None

        # read FFT params
        fftparams = fftparams.copy()
        fftparams.setdefault('method', 'medianmean')
        for param in ['fftlength', 'fftstride']:
            if hasattr(channel, param):
                fftparams[param] = float(getattr(channel, param))
            elif param in fftparams:
                fftparams[param] = float(fftparams[param])
        if hasattr(channel, 'stride'):
            stride = float(channel.stride)
        elif 'stride' in fftparams:
            stride = float(fftparams.pop('stride', 0))
        else:
            stride = None
        # get time-series data
        timeserieslist = get_timeseries(channel, new, config=config,
                                        cache=cache, query=query, nds=nds)
        # calculate spectrograms
        vprint("    Calculating spectrograms for %s" % str(channel))
        for ts in timeserieslist:
            fftparams.setdefault('fftlength', int(4096 * ts.dx.value))
            fftparams.setdefault('fftstride', fftparams['fftlength'] / 2)
            if not stride and fftparams['fftstride'] != fftparams['fftlength']:
                stride = fftparams['fftlength'] * 1.5
            elif not stride:
                stride = fftparams['fftlength']
            if abs(ts.span) < stride:
                continue
            specgram = ts.spectrogram(stride, **fftparams)
            if filter_:
                specgram = (specgram ** (1/2.)).filter(*filter_, inplace=True) ** 2
            globalv.SPECTROGRAMS[str(channel)].append(specgram)
            globalv.SPECTROGRAMS[str(channel)].coalesce()
            vprint('.')
        vprint('\n')

    # return correct data
    out = SpectrogramList()
    for specgram in globalv.SPECTROGRAMS[str(channel)]:
        for seg in segments:
            if abs(seg) < specgram.dt.value:
                continue
            if specgram.span.intersects(seg):
                if format in ['amplitude', 'asd']:
                    s = specgram.crop(*seg) ** (1/2.)
                else:
                    s = specgram.crop(*seg)
                if s.shape[0]:
                    out.append(s)
    return out.coalesce()


def get_spectrum(channel, segments, config=ConfigParser(), cache=None,
                 query=True, nds='guess', format='power', **fftparams):
    """Retrieve the time-series and generate a spectrogram of the given
    channel
    """
    channel = get_channel(channel)
    if isinstance(segments, DataQualityFlag):
        name = ','.join([str(channel), segments.name])
        segments = segments.active
    else:
        name = str(channel)
    cmin = '%s.min' % name
    cmax = '%s.max' % name

    if name not in globalv.SPECTRUM:
        speclist = get_spectrogram(channel, segments, config=config,
                                   cache=cache, query=query, nds=nds,
                                   **fftparams)
        specgram = speclist.join(gap='ignore')
        try:
            globalv.SPECTRUM[name] = specgram.percentile(50)
        except ValueError:
            globalv.SPECTRUM[name] = Spectrum([], channel=channel, f0=0, df=1)
            globalv.SPECTRUM[cmin] = globalv.SPECTRUM[name]
            globalv.SPECTRUM[cmax] = globalv.SPECTRUM[name]
        else:
            globalv.SPECTRUM[cmin] = specgram.percentile(5)
            globalv.SPECTRUM[cmax] = specgram.percentile(95)

    cmin = '%s.min' % name
    cmax = '%s.max' % name
    out = (globalv.SPECTRUM[name], globalv.SPECTRUM[cmin],
           globalv.SPECTRUM[cmax])
    if format in ['amplitude', 'asd']:
        out = [s ** (1/2.) for s in out]
    return out
