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
try:
    from configparser import (ConfigParser, NoSectionError, NoOptionError)
except ImportError:
    from ConfigParser import (ConfigParser, NoSectionError, NoOptionError)
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

import decorator
import numpy
import nds2
import warnings
import operator

from astropy import units

from glue import datafind
from glue.lal import Cache
from glue.segments import segmentlist

from gwpy.detector import Channel
from gwpy.segments import (DataQualityFlag, SegmentList, Segment)
from gwpy.timeseries import (TimeSeries, TimeSeriesList, TimeSeriesDict,
                             StateVector, StateVectorDict)
from gwpy.spectrum import Spectrum
from gwpy.spectrogram import SpectrogramList
from gwpy.io import nds as ndsio

from . import (globalv, version)
from .mode import *
from .utils import *
from .channels import get_channel

OPERATOR = {
    '*': operator.mul,
    '-': operator.sub,
    '+': operator.add,
    '/': operator.div,
}

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version


# -----------------------------------------------------------------------------
# data access

@decorator.decorator
def use_segmentlist(f, arg1, segments, *args, **kwargs):
    """Decorator a method to convert incoming segments into a `SegmentList`
    """
    if isinstance(segments, DataQualityFlag):
        segments = segments.active
    elif not isinstance(segments, segmentlist):
        segments = SegmentList([Segment(*x) for x in segments])
    return f(arg1, segments, *args, **kwargs)


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
    vprint('    Finding %s-%s frames for [%d, %d)...'
           % (ifo[0], frametype, int(gpsstart), int(gpsend)))
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
    vprint(' %d found.\n' % len(cache))
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
    if channel.frametype is None:
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


def find_types(site=None, match=None):
    """Query the DataFind server for frame types matching the given options
    """
    conn = datafind.GWDataFindHTTPConnection()
    return conn.find_types(site=site, match=match)


@use_segmentlist
def get_timeseries_dict(channels, segments, config=ConfigParser(),
                        cache=None, query=True, nds='guess', multiprocess=True,
                        frametype=None, statevector=False, return_=True,
                        **ioargs):
    """Retrieve the data for a set of channels
    """
    # separate channels by type
    if query:
        if frametype is not None:
            frametypes = {frametype: channels}
        else:
            frametypes = dict()
            allchannels = set([
                c for group in
                    map(lambda x: re_channel.findall(Channel(x).ndsname),
                        channels)
                for c in group])
            for channel in allchannels:
                channel = get_channel(channel)
                ifo = channel.ifo
                ftype = find_frame_type(channel)
                id_ = (ifo, ftype)
                if id_ in frametypes:
                    frametypes[id_].append(channel)
                else:
                    frametypes[id_] = [channel]
        for ftype, channellist in frametypes.iteritems():
            _get_timeseries_dict(channellist, segments, config=config,
                                 cache=cache, query=query, nds=nds,
                                 multiprocess=multiprocess, frametype=ftype[1],
                                 statevector=statevector, return_=False,
                                 **ioargs)
    if not return_:
        return
    else:
        out = OrderedDict()
        for channel in channels:
            channel = Channel(channel)
            chanstrs = re_channel.findall(channel.ndsname)
            chans = map(get_channel, chanstrs)
            tsdict = _get_timeseries_dict(chans, segments, config=config,
                                          query=False, statevector=statevector,
                                          **ioargs)
            tslist = [tsdict[Channel(c).ndsname] for c in chans]
            # if only one channel, simply append
            if len(chanstrs) == 1 and len(channel.ndsname) == len(chanstrs[0]):
                out[channel.ndsname] = tslist[0]
            # if one channel and some operator, do calculation
            elif len(chanstrs) == 1:
                stub = channel.ndsname[len(chanstrs[0]):].strip(' ')
                try:
                   op = OPERATOR[stub[0]]
                except KeyError as e:
                   print(channel.ndsname, chanstrs, stub)
                   e.args = ('Cannot parse math operator %r' % stub[0],)
                   raise
                value = float(stub[1:])
                out[channel.ndsname] = type(tslist[0])(*[op(ts, value)
                                                         for ts in tslist[0]])
            # if multiple channels
            else:
                # get union of segments for all sub-channels
                datasegs = reduce(operator.and_,
                                  [tsl.segments for tsl in tslist])
                # build meta-timeseries for all interseceted segments
                meta = TimeSeriesList()
                operators = [channel.name[m.span()[1]] for m in
                             list(re_channel.finditer(channel.ndsname))[:-1]]
                for seg in datasegs:
                    ts = get_timeseries(
                        chanstrs[0], SegmentList([seg]), config=config,
                        query=False, return_=True)[0]
                    ts.name = str(channel)
                    for op, ch in zip(operators, chanstrs[1:]):
                        try:
                            op = OPERATOR[op]
                        except KeyError as e:
                            e.args = ('Cannot parse math operator %r' % op,)
                            raise
                        data = get_timeseries(ch, SegmentList([seg]),
                                              config=config, query=False,
                                              return_=True)
                        ts = op(ts, data[0])
                    meta.append(ts)
                out[channel.ndsname] = meta
        return out


@use_segmentlist
def _get_timeseries_dict(channels, segments, config=ConfigParser(),
                         cache=None, query=True, nds='guess', frametype=None,
                         multiprocess=True, return_=True, statevector=False,
                         archive=True, **ioargs):
    """Internal method to retrieve the data for a set of like-typed
    channels using the :meth:`TimeSeriesDict.read` accessor.
    """
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
    dtype_ = dict()
    for channel in channels:
        try:
            filter_[channel.ndsname] = channel.filter
        except AttributeError:
            pass
        try:
            resample[channel] = float(channel.resample)
        except AttributeError:
            pass
        if channel.dtype is None:
            dtype_[channel] = ioargs.get('dtype')
        else:
            dtype_[channel] = channel.dtype

    # work out whether to use NDS or not
    if nds == 'guess':
        nds = 'LIGO_DATAFIND_SERVER' not in os.environ

    # read new data
    query &= (abs(new) > 0)
    if cache is not None:
        query &= len(cache) > 0
    if query:
        for channel in channels:
            globalv.DATA.setdefault(channel.ndsname, ListClass())
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
                else:
                    raise
            frametype = source = 'nds'
            ndstype = channels[0].type
        elif nds:
            ndsconnection = None
            frametype = source = 'nds'
            ndstype = channels[0].type
        # or find frame type and check cache
        else:
            ifo = channels[0].ifo
            frametype = frametype or channels[0].frametype
            if frametype is not None and frametype.endswith('%s_M' % ifo):
                new = type(new)([s for s in new if abs(s) >= 60.])
            elif frametype is not None and frametype.endswith('%s_T' % ifo):
                new = type(new)([s for s in new if abs(s) >= 1.])
            #elif ((globalv.NOW - new[0][0]) < 86400 * 10 and
            #      frametype == '%s_R' % ifo and
            #      find_types(site=ifo[0], match='_C\Z')):
            #    frametype = '%s_C' % ifo
            if cache is not None:
                fcache = cache.sieve(ifos=ifo[0], description=frametype,
                                     exact_match=True)
            else:
                fcache = Cache()
            if (cache is None or len(fcache) == 0) and len(new):
                span = new.extent().protract(8)
                fcache = find_frames(ifo, frametype, span[0], span[1],
                                     config=config, gaps='ignore')
            # parse discontiguous cache blocks and rebuild segment list
            cachesegments = find_cache_segments(fcache)
            new &= cachesegments
            source = 'frames'
        for channel in channels:
            channel.frametype = frametype

        # check whether each channel exists for all new times already
        qchannels = []
        qresample = {}
        qdtype = {}
        for channel in channels:
            oldsegs = globalv.DATA.get(channel.ndsname,
                                       ListClass()).segments
            if abs(new - oldsegs) != 0:
                qchannels.append(channel)
                if channel in resample:
                    qresample[channel] = resample[channel]
                qdtype[channel] = dtype_.get(channel, ioargs.get('dtype'))
        ioargs['dtype'] = qdtype

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
            vprint("    Fetching data (from %s) for %d channels [%s]"
                   % (source, len(qchannels), nds and ndstype or frametype))
        for segment in new:
            # force reading integer-precision segments
            segment = type(segment)(int(segment[0]), int(segment[1]))
            if abs(segment) < 1:
                continue
            if nds:
                tsd = DictClass.fetch(qchannels, segment[0], segment[1],
                                      connection=ndsconnection, type=ndstype,
                                      **ioargs)
            else:
                if nproc > 1:
                    for c in qchannels:
                        if c.ndsname in filter_:
                            del c.filter
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
                # set minute trend times modulo 60 from GPS 0
                if re.match('(?:(.*)_)?[A-Z]\d_M', frametype):
                    segstart = int(segment[0]) // 60 * 60
                    segend = int(segment[1]) // 60 * 60
                    if segend >= segment[1]:
                        segend -= 60
                    # and ignore segments shorter than 1 full average
                    if (segend - segstart) < 60:
                        continue
                else:
                    segstart, segend = map(float, segment)
                tsd = DictClass.read(segcache, qchannels, format='lcf',
                                     start=segstart, end=segend, type=ctype,
                                     nproc=nproc, resample=qresample,
                                     verbose=verbose, **ioargs)
                for c in qchannels:
                    if c.ndsname in filter_:
                        c.filter = filter_[c.ndsname]
            for (channel, data) in tsd.iteritems():
                if (channel.ndsname in globalv.DATA and
                    data.span in globalv.DATA[channel.ndsname].segments):
                    continue
                for seg in globalv.DATA[channel.ndsname].segments:
                    if seg.intersects(data.span):
                        data = data.crop(*(data.span - seg))
                        break
                data.channel = channel
                if channel.ndsname in filter_:
                    if callable(filter_[channel.ndsname]):
                        data = filter_[channel.ndsname](data)
                    else:
                        data = data.filter(*filter_[channel.ndsname])
                if isinstance(data, StateVector) or ':GRD-' in str(channel):
                    try:
                        data.unit = units.dimensionless_unscaled
                    except AttributeError:
                        data._unit = units.dimensionless_unscaled
                    if hasattr(channel, 'bits'):
                        data.bits = channel.bits
                elif data.unit is None:
                    data._unit = channel.unit
                # XXX: HACK for failing unit check
                if len(globalv.DATA[channel.ndsname]):
                    data._unit = globalv.DATA[channel.ndsname][-1].unit
                # append and coalesce
                globalv.DATA[channel.ndsname].append(data)
                globalv.DATA[channel.ndsname].coalesce()
            vprint('.')
        if len(new):
            vprint("\n")

    if not return_:
        return

    # return correct data
    out = OrderedDict()
    for channel in channels:
        data = ListClass()
        if channel.ndsname not in globalv.DATA:
            out[channel.ndsname] = ListClass()
        else:
            for ts in globalv.DATA[channel.ndsname]:
                for seg in segments:
                    if abs(seg) == 0 or abs(seg) < ts.dt.value:
                        continue
                    if ts.span.intersects(seg):
                        cropped = ts.crop(float(seg[0]), float(seg[1]),
                                          copy=False)
                        if cropped.size:
                            data.append(cropped)
        out[channel.ndsname] = data.coalesce()
    return out


@use_segmentlist
def get_timeseries(channel, segments, config=ConfigParser(), cache=None,
                   query=True, nds='guess', multiprocess=True,
                   frametype=None, statevector=False, return_=True,
                   **ioargs):
    """Retrieve the data (time-series) for a given channel
    """
    channel = get_channel(channel)
    out = get_timeseries_dict([channel.ndsname], segments, config=config,
                              cache=cache, query=query, nds=nds,
                              multiprocess=multiprocess, frametype=frametype,
                              statevector=statevector, return_=return_,
                              **ioargs)
    if return_:
        return out[channel.ndsname]
    return


@use_segmentlist
def get_spectrogram(channel, segments, config=ConfigParser(), cache=None,
                    query=True, nds='guess', format='power', return_=True,
                    frametype=None, multiprocess=True, **fftparams):
    """Retrieve the time-series and generate a spectrogram of the given
    channel
    """
    channel = get_channel(channel)

    # read data for all sub-channels
    specs = []
    channels = re_channel.findall(channel.ndsname)
    for c in channels:
        specs.append(_get_spectrogram(c, segments, config=config, cache=cache,
                                      query=query, nds=nds, format=format,
                                      return_=return_,
                                      multiprocess=multiprocess,
                                      **fftparams))
    if return_ and len(channels) == 1:
        return specs[0]
    elif return_:
        # get union of segments for all sub-channels
        datasegs = reduce(operator.and_, [sgl.segments for sgl in specs])
        # build meta-spectrogram for all interseceted segments
        out = SpectrogramList()
        operators = [channel.name[m.span()[1]] for m in
                     list(re_channel.finditer(channel.ndsname))[:-1]]
        for seg in datasegs:
            sg = _get_spectrogram(channels[0], SegmentList([seg]),
                                  config=config, query=False, format=format,
                                  return_=True)[0]
            sg.name = str(channel)
            for op, ch in zip(operators, channels[1:]):
                try:
                    op = OPERATOR[op]
                except KeyError as e:
                    e.args = ('Cannot parse math operator %r' % op,)
                    raise
                data = _get_spectrogram(ch, SegmentList([seg]),
                                        config=config, query=False,
                                        format=format, return_=True)
                sg = op(sg, data[0])
            out.append(sg)
        return out


@use_segmentlist
def _get_spectrogram(channel, segments, config=ConfigParser(), cache=None,
                     query=True, nds='guess', format='power', return_=True,
                     frametype=None, multiprocess=True, method='median-mean',
                     **fftparams):
    channel = get_channel(channel)
    if format in ['rayleigh']:
        method = format
    key = '%s,%s' % (channel.ndsname, method)
    # read segments from global memory
    havesegs = globalv.SPECTROGRAMS.get(key,
                                        SpectrogramList()).segments
    new = segments - havesegs

    # get processes
    if multiprocess is True:
        nproc = count_free_cores()
    elif multiprocess is False:
        nproc = 1
    else:
        nproc = count_free_cores(multiprocess)

    globalv.SPECTROGRAMS.setdefault(key, SpectrogramList())

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
        if stride is not None:
            tmp = type(new)()
            for s in new:
                if abs(s) < stride:
                    continue
                else:
                    d = float(abs(s))
                    tmp.append(type(s)(s[0], s[0] + d//stride * stride))
            new = tmp
        timeserieslist = get_timeseries(channel, new, config=config,
                                        cache=cache, frametype=frametype,
                                        query=query, nds=nds)
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
                specgram = ts.spectrogram(stride, nproc=nproc, method=method,
                                          **fftparams)
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
            if filter_ and method not in ['rayleigh']:
                specgram = (specgram ** (1/2.)).filter(*filter_, inplace=True) ** 2
            if specgram.unit is None:
                specgram._unit = channel.unit
            elif len(globalv.SPECTROGRAMS[key]):
                specgram._unit = globalv.SPECTROGRAMS[key][-1].unit
            globalv.SPECTROGRAMS[key].append(specgram)
            globalv.SPECTROGRAMS[key].coalesce()
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
                if format in ['amplitude', 'asd']:
                    s = specgram.crop(*seg) ** (1/2.)
                else:
                    s = specgram.crop(*seg)
                    # XXX FIXME: this corrects the bias offset in Rayleigh
                    if format in ['rayleigh']:
                        med = numpy.median(s.value)
                        s /= med
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
    name += ',%s' % format
    cmin = '%s.min' % name
    cmax = '%s.max' % name

    if name not in globalv.SPECTRUM:
        vprint("    Calculating 5/50/95 percentile spectra for %s"
               % name.rsplit(',', 1)[0])
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
            globalv.SPECTRUM[name] = Spectrum([], channel=channel, f0=0, df=1,
                                              unit=units.Unit(''))
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


@use_segmentlist
def get_spectrograms(channels, segments, config=ConfigParser(), cache=None,
                     query=True, nds='guess', format='power', return_=True,
                     method='median-mean', frametype=None, multiprocess=True,
                     **fftparams):
    """Get spectrograms for multiple channels
    """
    channels = map(get_channel, channels)
    # get timeseries data in bulk
    if query:
        if format in ['rayleigh']:
            method_ = format
        else:
            method_ = method
        keys = ['%s,%s' % (channel.ndsname, method_) for channel in channels]
        havesegs = reduce(operator.and_, (globalv.SPECTROGRAMS.get(
            key, SpectrogramList()).segments for key in keys))
        new = segments - havesegs
        strides = set([getattr(c, 'stride', 0) for c in channels])
        if len(strides) == 1:
            stride = strides.pop()
            new = type(new)([s for s in new if abs(s) >= stride])
        get_timeseries_dict(channels, new, config=config, cache=cache,
                            frametype=frametype, nds=nds, return_=False)
    # loop over channels and generate spectrograms
    out = OrderedDict()
    for channel in channels:
         out[channel] = get_spectrogram(
             channel, segments, config=config, cache=cache, query=query,
             nds=nds, format=format, multiprocess=multiprocess,
             return_=return_, method=method, **fftparams)
    return out
