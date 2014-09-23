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

import operator
import os
import urllib2
from math import (floor, ceil, pi, sqrt)
from Queue import Queue
import threading
try:
    from configparser import (ConfigParser, NoSectionError, NoOptionError)
except ImportError:
    from ConfigParser import (ConfigParser, NoSectionError, NoOptionError)

import numpy
import nds2

from astropy import units

from glue import datafind
from glue.lal import Cache

from gwpy.detector import Channel
from gwpy.segments import (DataQualityFlag, SegmentList)
from gwpy.timeseries import (TimeSeries, TimeSeriesList, TimeSeriesDict,
                             StateVector, StateVectorDict)
from gwpy.spectrum import Spectrum
from gwpy.spectrogram import SpectrogramList
from gwpy.io import nds as ndsio

from . import (globalv, version)
from .mode import *
from .utils import *

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


class ThreadChannelQuery(threading.Thread):
    """Threaded CIS `Channel` query.
    """
    def __init__(self, inqueue, outqueue):
        threading.Thread.__init__(self)
        self.in_ = inqueue
        self.out = outqueue

    def run(self):
        i, channel = self.in_.get()
        self.in_.task_done()
        try:
            self.out.put((i, get_channel(channel, False)))
        except Exception as e:
            self.out.put(e)
        self.out.task_done()


def get_channel(channel, find_trend_source=True):
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
        found = globalv.CHANNELS.sieve(name=name, type=type_, exact_match=True)
    else:
        type_ = isinstance(channel, Channel) and channel.type or None
        sr = isinstance(channel, Channel) and channel.sample_rate or None
        name = str(channel)
        found = globalv.CHANNELS.sieve(name=str(channel), type=type_,
                                       sample_rate=sr, exact_match=True)
    if len(found) == 1:
        return found[0]
    elif len(found) > 1:
        raise ValueError("Ambiguous channel request '%s', multiple existing "
                         "channels recovered:\n    %s"
                         % (str(channel),
                            '\n    '.join([c.ndsname for c in found])))
    else:
        try:
            if not re_channel.match(name):
                raise ValueError()
            # trends are not stored in CIS, but try and get their raw source
            if re.search('.[a-z]+\Z', name):
                raise TypeError()
            new = Channel.query(name)
        except TypeError:
            # set default trend type based on mode
            if type_ is None and globalv.MODE == SUMMARY_MODE_GPS:
                type_ = 's-trend'
            elif type_ is None:
                type_ = 'm-trend'
            name += ',%s' % type_
            new = Channel(name)
            if find_trend_source:
                try:
                    source = get_channel(name.rsplit('.', 1)[0])
                except ValueError:
                    pass
                else:
                    new.url = source.url
            # determine sample rate for trends
            if type_ == 'm-trend':
                new.sample_rate = 1/60.
            elif type_ == 's-trend':
                new.sample_rate = 1
        except (ValueError, urllib2.URLError):
            new = Channel(str(channel))
        else:
            new.name = str(channel)
        globalv.CHANNELS.append(new)
        try:
            return get_channel(new)
        except RuntimeError as e:
            if 'maximum recursion depth' in str(e):
                raise RuntimeError("Recursion error while access channel "
                                   "information for %s" % str(channel))
            else:
                raise


def get_channels(channels):
    """Multi-threaded channel query
    """
    if len(channels) == 0:
        return []

    # set up Queues
    inqueue = Queue()
    outqueue = Queue()

    # open threads
    for i in range(len(channels)):
        t = ThreadChannelQuery(inqueue, outqueue)
        t.setDaemon(True)
        t.start()

    # populate input queue
    for i, c in enumerate(channels):
        inqueue.put((i, c))

    # block
    inqueue.join()
    outqueue.join()
    result = []
    for i in range(len(channels)):
        c = outqueue.get()
        if isinstance(c, Exception):
            raise c
        else:
            result.append(c)
    return zip(*sorted(result, key=lambda (idx, chan): idx))[1]


def override_sample_rate(channel, rate):
    try:
        cid = globalv.CHANNELS.find(channel.name)
    except ValueError:
        print('Failed to reset sample_rate for %s' % channel.name)
        pass
    else:
        globalv.CHANNELS[cid].sample_rate = rate


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

    # XXX HACK: LLO changed frame types on Dec 6 2013:
    LLOCHANGE = 1070291904

    if re.match('L1_{CRMT}', frametype) and gpsstart < LLOCHANGE:
        frametype = frametype[-1]

    # query frames
    ifo = ifo[0].upper()
    gpsstart = int(floor(gpsstart))
    gpsend = int(ceil(min(globalv.NOW, gpsend)))
    if gpsend <= gpsstart:
        return Cache()

    try:
        cache = dfconn.find_frame_urls(ifo[0].upper(), frametype, gpsstart,
                                       gpsend, urltype=urltype, on_gaps=gaps)
    except RuntimeError as e:
        if 'Invalid GPS times' in str(e):
            e2 = str(e2) + ': %d ... %d' % (gpsstart, gpsend)
            raise RuntimeError(e2)
        else:
            raise

    # XXX: if querying for day of LLO frame type change, do both
    if (ifo[0].upper() == 'L' and frametype in ['C', 'R', 'M', 'T'] and
            gpsstart < LLOCHANGE < gpsend):
        start = len(cache) and cache[-1].segment[1] or gpsstart
        if start < gpsend:
            cache.extend(dfconn.find_frame_urls(ifo[0].upper(),
                                                'L1_%s' % frametype, start,
                                                gpsend, urltype=urltype,
                                                on_gaps=gaps)[1:])
    cache, _ = cache.checkfilesexist()
    return cache


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


def find_frame_type(channel):
    try:
        return channel.frametype
    except AttributeError:
        try:
            ndstype = ndsio.NDS2_CHANNEL_TYPE[channel.type]
        except (AttributeError, KeyError):
            ndstype = channel.type = nds2.channel.CHANNEL_TYPE_RAW
        if ndstype == nds2.channel.CHANNEL_TYPE_MTREND:
            ftype = 'M'
        elif ndstype == nds2.channel.CHANNEL_TYPE_STREND:
            ftype = 'T'
        elif ndstype == nds2.channel.CHANNEL_TYPE_RDS:
            ftype = 'LDAS_C02_L2'
        elif ndstype == nds2.channel.CHANNEL_TYPE_ONLINE:
            ftype = 'lldetchar'
        else:
            ftype = 'R'
        channel.frametype = '%s1_%s' % (channel.ifo[0], ftype)
        return channel.frametype


def get_timeseries_dict(channels, segments, config=ConfigParser(),
                        cache=None, query=True, nds='guess',
                        multiprocess=True, statevector=False, return_=True):
    """Retrieve the data for a set of channels
    """
    # separate channels by type
    frametypes = dict()
    for channel in channels:
        channel = get_channel(channel)
        ifo = channel.ifo
        ftype = find_frame_type(channel)
        id_ = (ifo, ftype)
        if id_ in frametypes:
            frametypes[id_].append(channel)
        else:
            frametypes[id_] = [channel]
    out = dict()
    for channellist in frametypes.itervalues():
        data = _get_timeseries_dict(channellist, segments, config=config,
                                    cache=cache, query=query, nds=nds,
                                    multiprocess=multiprocess,
                                    statevector=statevector)
        if not return_:
            continue
        for (name, tslist) in data.iteritems():
            if name in out:
                out[name].extend(tslist)
            else:
                out[name] = tslist
    return out


def _get_timeseries_dict(channels, segments, config=ConfigParser(),
                         cache=None, query=True, nds='guess',
                         multiprocess=True, return_=True, statevector=False):
    """Internal method to retrieve the data for a set of like-typed
    channels using the :meth:`TimeSeriesDict.read` accessor.
    """
    if isinstance(segments, DataQualityFlag):
        segments = segments.active
    channels = map(get_channel, channels)

    # set classes
    if statevector:
        SeriesClass = StateVector
        ListClass = TimeSeriesList
        DictClass = StateVectorDict
    else:
        SeriesClass = TimeSeries
        ListClass = TimeSeriesList
        DictClass = TimeSeriesDict

    # read segments from global memory
    havesegs = reduce(operator.and_,
                      (globalv.DATA.get(channel.ndsname,
                                        ListClass()).segments
                       for channel in channels))
    new = segments - havesegs

    # get processes
    if multiprocess is True:
        nproc = count_free_cores()
    elif multiprocess is False:
        nproc = 1
    else:
        nproc = count_free_cores(multiprocess)

    if globalv.VERBOSE and not multiprocess:
        verbose = '    '
    else:
        verbose = False

    # read channel information
    filter_ = dict()
    resample = dict()
    for channel in channels:
        try:
            filter_[channel.ndsname] = channel.filter
        except AttributeError:
            pass
        try:
            resample[channel] = float(channel.resample)
        except AttributeError:
            pass

    # work out whether to use NDS or not
    if nds == 'guess':
        nds = 'LIGO_DATAFIND_SERVER' not in os.environ

    # read new data
    for channel in channels:
        globalv.DATA.setdefault(channel.ndsname, ListClass())
    query &= (abs(new) > 0)
    if cache is not None:
        query &= len(cache) > 0
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
            ndstype = channels[0].type
        elif nds:
            ndsconnection = None
            source = 'nds'
            ndstype = channels[0].type
        # or find frame type and check cache
        else:
            ifo = channels[0].ifo
            ftype = channels[0].frametype
            if ftype == '%s_M' % ifo:
                new = type(new)([s for s in new if abs(s) >= 60.])
            elif ftype == '%s_T' % ifo:
                new = type(new)([s for s in new if abs(s) >= 1.])
            elif ((globalv.NOW - new[0][0]) < 86400 * 20 and
                  ftype == '%s_R' % ifo):
                ftype = '%s_C' % ifo
            if cache is not None:
                fcache = cache.sieve(ifos=ifo[0], description=ftype,
                                     exact_match=True)
            else:
                fcache = Cache()
            if (cache is None or len(fcache) == 0) and len(new):
                span = new.extent().protract(8)
                fcache = find_frames(ifo, ftype, span[0], span[1],
                                     config=config, gaps='ignore')
            # parse discontiguous cache blocks and rebuild segment list
            cachesegments = find_cache_segments(fcache)
            new &= cachesegments
            source = 'frames'

        # check whether each channel exists for all new times already
        qchannels = []
        qresample = {}
        for channel in channels:
            oldsegs = globalv.DATA.get(channel.ndsname,
                                       ListClass()).segments
            if abs(new - oldsegs) != 0:
                qchannels.append(channel)
                if channel in resample:
                    qresample[channel] = resample[channel]

        # find channel type
        if not nds:
            ctype = set()
            for channel in qchannels:
                try:
                    ctype.add(channel.ctype)
                except AttributeError:
                    channel.ctype = 'adc'
                    ctype.add(channel.ctype)
                    continue
            if len(ctype) == 1:
                ctype = list(ctype)[0]
            else:
                ctype = None
        # loop through segments, recording data for each
        if len(new) and nproc > 1:
            vprint("    Fetching data (from %s) for %d channels"
                   % (source, len(qchannels)))
        for segment in new:
            if nds:
                tsd = DictClass.fetch(qchannels, segment[0], segment[1],
                                      connection=ndsconnection, type=ndstype)
            else:
                # pad resampling
                if segment[1] == cachesegments[-1][1] and qresample:
                    resamplepad = 8
                    if abs(segment) <= resamplepad:
                        continue
                    segment = type(segment)(segment[0],
                                            segment[1] - resamplepad)
                    segcache = fcache.sieve(
                                   segment=segment.protract(resamplepad))
                else:
                    segcache = fcache.sieve(segment=segment)
                tsd = DictClass.read(segcache, qchannels, format='lcf',
                                     start=float(segment[0]),
                                     end=float(segment[1]), type=ctype,
                                     nproc=nproc,
                                     resample=qresample, verbose=verbose)
            for (channel, data) in tsd.iteritems():
                if channel.sample_rate is None:
                    override_sample_rate(channel, data.sample_rate)
                if (channel.ndsname in globalv.DATA and
                    data.span in globalv.DATA[channel.ndsname].segments):
                    continue
                for seg in globalv.DATA[channel.ndsname].segments:
                    if seg.intersects(data.span):
                        data = data.crop(*(data.span - seg))
                        break
                data.channel = channel
                if channel.ndsname in filter_:
                    data = data.filter(*filter_[channel.ndsname])
                if isinstance(data, StateVector):
                    data.unit = units.dimensionless_unscaled
                    if hasattr(channel, 'bits'):
                        data.bits = channel.bits
                # XXX: HACK for failing unit check
                globalv.DATA[channel.ndsname].append(data)
                try:
                    globalv.DATA[channel.ndsname].coalesce()
                except ValueError as e:
                    if not 'units do not match' in str(e):
                        raise
                    warnings.warn(str(e))
                    globalv.DATA[channel.ndsname][-1].unit = (
                        globalv.DATA[channel.ndsname][0].unit)
                    globalv.DATA[channel.ndsname].coalesce()
            vprint('.')
        if len(new):
            vprint("\n")

    if not return_:
        return

    # return correct data
    out = dict()
    for channel in channels:
        data = ListClass()
        for ts in globalv.DATA[channel.ndsname]:
            for seg in segments:
                if abs(seg) == 0 or abs(seg) < ts.dt.value:
                    continue
                if ts.span.intersects(seg):
                    cropped = ts.crop(float(seg[0]), float(seg[1]), copy=False)
                    if cropped.size:
                        data.append(cropped)
        out[channel.ndsname] = data.coalesce()
    return out


def get_timeseries(channel, segments, config=ConfigParser(), cache=None,
                   query=True, nds='guess', multiprocess=True,
                   statevector=False, return_=True):
    """Retrieve the data (time-series) for a given channel
    """
    channel = get_channel(channel)
    out = get_timeseries_dict([channel.ndsname], segments, config=config,
                              cache=cache, query=query, nds=nds,
                              multiprocess=multiprocess,
                              statevector=statevector, return_=return_)
    if return_:
        return out[channel.ndsname]
    return


def get_spectrogram(channel, segments, config=ConfigParser(), cache=None,
                    query=True, nds='guess', format='power', return_=True,
                    multiprocess=True, **fftparams):
    """Retrieve the time-series and generate a spectrogram of the given
    channel
    """
    if isinstance(segments, DataQualityFlag):
        segments = segments.active
    channel = get_channel(channel)
    # read segments from global memory
    havesegs = globalv.SPECTROGRAMS.get(channel.ndsname,
                                        SpectrogramList()).segments
    new = segments - havesegs

    # get processes
    if multiprocess is True:
        nproc = count_free_cores()
    elif multiprocess is False:
        nproc = 1
    else:
        nproc = count_free_cores(multiprocess)

    globalv.SPECTROGRAMS.setdefault(channel.ndsname, SpectrogramList())

    query &= abs(new) != 0
    if query:
        # read channel information
        try:
            filter_ = channel.frequency_response
        except AttributeError:
            filter_ = None
        else:
            if isinstance(filter_, str):
                filter_ = eval(filter_)

        # read FFT params
        fftparams = fftparams.copy()
        fftparams.setdefault('method', 'median-mean')
        for param in ['fftlength', 'overlap']:
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
                                        cache=cache, query=False, nds=nds)
        # calculate spectrograms
        if len(timeserieslist):
            vprint("    Calculating spectrograms for %s" % str(channel))
        for ts in timeserieslist:
            fftparams.setdefault('fftlength', int(4096 * ts.dx.value))
            fftparams.setdefault('overlap', fftparams['fftlength'] / 2)
            if not stride and fftparams['overlap'] != fftparams['fftlength']:
                stride = fftparams['fftlength'] * 1.5
            elif not stride:
                stride = fftparams['fftlength']
            if abs(ts.span) < stride:
                continue
            try:
                specgram = ts.spectrogram(stride, nproc=nproc,
                                          **fftparams)
            except ZeroDivisionError:
                if stride == 0:
                    raise ZeroDivisionError("Spectrogram stride is 0")
                elif fftparams['fftlength'] == 0:
                    raise ZeroDivisionError("FFT length is 0")
                else:
                    raise
            if filter_:
                specgram = (specgram ** (1/2.)).filter(*filter_, inplace=True) ** 2
            globalv.SPECTROGRAMS[channel.ndsname].append(specgram)
            globalv.SPECTROGRAMS[channel.ndsname].coalesce()
            vprint('.')
        if len(timeserieslist):
            vprint('\n')

    if not return_:
        return

    # return correct data
    out = SpectrogramList()
    for specgram in globalv.SPECTROGRAMS[channel.ndsname]:
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
                 query=True, nds='guess', format='power', return_=True,
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
    cmin = '%s.min' % name
    cmax = '%s.max' % name

    if name not in globalv.SPECTRUM:
        vprint("    Calculating 5/50/95 percentile spectra for %s"
               % name.rsplit(',', 1)[0])
        speclist = get_spectrogram(channel, segments, config=config,
                                   cache=cache, query=False, nds=nds,
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
        vprint(".\n")

    if not return_:
        return

    cmin = '%s.min' % name
    cmax = '%s.max' % name
    out = (globalv.SPECTRUM[name], globalv.SPECTRUM[cmin],
           globalv.SPECTRUM[cmax])
    if format in ['amplitude', 'asd']:
        out = [s ** (1/2.) for s in out]
    return out


def add_timeseries(timeseries, key=None, coalesce=True):
    """Add a `TimeSeries` to the global memory cache
    """
    if key is None:
        key = timeseries.name or timeseries.channel.ndsname
    globalv.DATA.setdefault(key, TimeSeriesList())
    globalv.DATA[key].append(timeseries)
    if coalesce:
        globalv.DATA[key].coalesce()


def add_spectrogram(specgram, key=None, coalesce=True):
    """Add a `Spectrogram` to the global memory cache
    """
    if key is None:
        key = specgram.name or str(specgram.channel)
    globalv.SPECTROGRAMS.setdefault(key, SpectrogramList())
    globalv.SPECTROGRAMS[key].append(specgram)
    if coalesce:
        globalv.SPECTROGRAMS[key].coalesce()

