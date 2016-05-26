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
import re
import os
import urllib2
from math import (floor, ceil, pi, sqrt)
from time import sleep
try:
    from configparser import (ConfigParser, NoSectionError, NoOptionError)
except ImportError:
    from ConfigParser import (ConfigParser, NoSectionError, NoOptionError)
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict
from itertools import (izip, izip_longest)

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
from gwpy import astro
from gwpy.segments import (DataQualityFlag, SegmentList, Segment)
try:
    from gwpy.timeseries import (TimeSeries, TimeSeriesList, TimeSeriesDict,
                                 StateVector, StateVectorList, StateVectorDict)
except ImportError:
    from gwpy.timeseries import (TimeSeries, TimeSeriesList, TimeSeriesDict,
                                 StateVector, StateVectorDict)
    StateVectorList = TimeSeriesList
try:
    from gwpy.frequencyseries import FrequencySeries
except ImportError:
    from gwpy.spectrum import Spectrum as FrequencySeries
from gwpy.spectrogram import SpectrogramList
from gwpy.io import nds as ndsio

from . import globalv
from .mode import *
from .utils import *
from .channels import (get_channel, update_missing_channel_params)

OPERATOR = {
    '*': operator.mul,
    '-': operator.sub,
    '+': operator.add,
    '/': operator.div,
}

FRAMETYPE_REGEX = {
    'commissioning': re.compile('[A-Z][0-9]_C\Z'),
    'science': re.compile('[A-Z][0-9]_R\Z'),
    'second-trend': re.compile('[A-Z][0-9]_T\Z'),
    'minute-trend': re.compile('[A-Z][0-9]_M\Z'),
    'calibrated h(t) version 0': re.compile('[A-Z][0-9]_HOFT_C00\Z'),
    'calibrated h(t) version 1': re.compile('[A-Z][0-9]_HOFT_C01\Z'),
    'calibrated h(t) version 2': re.compile('[A-Z][0-9]_HOFT_C02\Z'),
    'DMT SenseMon on GDS h(t)': re.compile('SenseMonitor_hoft_[A-Z][0-9]_M\Z'),
    'DMT SenseMon on front-end h(t)': re.compile(
        'SenseMonitor_CAL_[A-Z][0-9]_M\Z'),
}

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


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
                urltype='file', gaps='warn', onerror='raise'):
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
            try:
                host, port = host.rsplit(':', 1)
            except ValueError:
                port = None
            else:
                port = int(port)
    else:
        port = config.getint('datafind', 'port')
    # get credentials
    if port == 80:
        cert = None
        key = None
    else:
        cert, key = datafind.find_credential()

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

    def _query():
        if cert is not None:
            dfconn = datafind.GWDataFindHTTPSConnection(
                host=host, port=port, cert_file=cert, key_file=key)
        else:
            dfconn = datafind.GWDataFindHTTPConnection(host=host, port=port)
        return dfconn.find_frame_urls(ifo[0].upper(), frametype, gpsstart,
                                      gpsend, urltype=urltype, on_gaps=gaps)
    try:
        cache = _query()
    except RuntimeError as e:
        sleep(1)
        try:
            cache = _query()
        except RuntimeError:
            if 'Invalid GPS times' in str(e):
                e.args = ('%s: %d ... %s' % (str(e), gpsstart, gpsend),)
            if onerror in ['ignore', None]:
                pass
            elif onerror in ['warn']:
                warnings.warn('Caught %s: %s'
                              % (type(e).__name__, str(e)))
            else:
                raise
            cache = Cache()

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
        if channel.ifo == 'C1':
            channel.frametype = ftype
        else:
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
                        datafind_error='raise', **ioargs):
    """Retrieve the data for a set of channels
    """
    # separate channels by type
    if query:
        if frametype is not None:
            frametypes = {(None, frametype): channels}
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
                                 datafind_error=datafind_error, **ioargs)
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
                meta = type(globalv.DATA[chans[0].ndsname])()
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
                         archive=True, datafind_error='raise', **ioargs):
    """Internal method to retrieve the data for a set of like-typed
    channels using the :meth:`TimeSeriesDict.read` accessor.
    """
    channels = map(get_channel, channels)

    # set classes
    if statevector:
        SeriesClass = StateVector
        ListClass = StateVectorList
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
        nproc = multiprocess

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
                                     config=config, gaps='ignore',
                                     onerror=datafind_error)
                if len(fcache) == 0 and frametype == '%s_R' % ifo:
                    frametype = '%s_C' % ifo
                    vprint("    Moving to backup frametype %s\n" % frametype)
                    fcache = find_frames(ifo, frametype, span[0], span[1],
                                         config=config, gaps='ignore',
                                         onerror=datafind_error)

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
                if (re.match('(?:(.*)_)?[A-Z]\d_M', frametype) or
                        (ifo == 'C1' and frametype == 'M')):
                    segstart = int(segment[0]) // 60 * 60
                    segend = int(segment[1]) // 60 * 60
                    if segend >= segment[1]:
                        segend -= 60
                    # and ignore segments shorter than 1 full average
                    if (segend - segstart) < 60:
                        continue
                else:
                    segstart, segend = map(float, segment)
                # pull filters out because they can break multiprocessing
                if nproc > 1:
                    for c in qchannels:
                        if c.ndsname in filter_:
                            del c.filter
                # read data
                tsd = DictClass.read(segcache, qchannels,
                                     start=segstart, end=segend, type=ctype,
                                     nproc=nproc, resample=qresample,
                                     verbose=verbose, **ioargs)
                # put filters back
                for c in qchannels:
                    if c.ndsname in filter_:
                        c.filter = filter_[c.ndsname]
            for (channel, data) in tsd.iteritems():
                channel = update_missing_channel_params(data.channel)
                if (channel.ndsname in globalv.DATA and
                    data.span in globalv.DATA[channel.ndsname].segments):
                    continue
                if data.unit is None:
                    data.unit = 'undef'
                for seg in globalv.DATA[channel.ndsname].segments:
                    if seg.intersects(data.span):
                        data = data.crop(*(data.span - seg))
                        break
                try:
                    filt = filter_[channel.ndsname]
                except KeyError:
                    pass
                else:
                    # filter with function
                    if callable(filt):
                        data = filt(data)
                    # filter with gain
                    elif (isinstance(filt, tuple) and len(filt) == 3 and
                              len(filt[0] + filt[1]) == 0):
                        try:
                            data *= filt[2]
                        except TypeError:
                            data = data * filt[2]
                    # filter zpk
                    elif isinstance(filt, tuple):
                        data = data.filter(*filt)
                    # filter fail
                    else:
                        raise ValueError("Cannot parse filter for %s: %r"
                                         % (channel.ndsname,
                                            filt))
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
                        common = map(float, ts.span & seg)
                        cropped = ts.crop(*common, copy=False)
                        if cropped.size:
                            data.append(cropped)
        out[channel.ndsname] = data.coalesce()
    return out


@use_segmentlist
def get_timeseries(channel, segments, config=ConfigParser(), cache=None,
                   query=True, nds='guess', multiprocess=True,
                   frametype=None, statevector=False, return_=True,
                   datafind_error='raise', **ioargs):
    """Retrieve the data (time-series) for a given channel
    """
    channel = get_channel(channel)
    out = get_timeseries_dict([channel.ndsname], segments, config=config,
                              cache=cache, query=query, nds=nds,
                              multiprocess=multiprocess, frametype=frametype,
                              statevector=statevector, return_=return_,
                              datafind_error=datafind_error, **ioargs)
    if return_:
        return out[channel.ndsname]
    return


@use_segmentlist
def get_spectrogram(channel, segments, config=ConfigParser(), cache=None,
                    query=True, nds='guess', format='power', return_=True,
                    frametype=None, multiprocess=True, datafind_error='raise',
                    **fftparams):
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
                                      return_=return_, frametype=frametype,
                                      multiprocess=multiprocess,
                                      datafind_error=datafind_error,
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
                try:
                    sg = op(sg, data[0])
                except ValueError as e:
                    if 'could not be broadcast together' in str(e):
                        s = min(sg.shape[0], data[0].shape[0])
                        sg = op(sg[:s], data[0][:s])
                    else:
                        raise
            out.append(sg)
        return out


@use_segmentlist
def _get_spectrogram(channel, segments, config=ConfigParser(), cache=None,
                     query=True, nds='guess', format='power', return_=True,
                     frametype=None, multiprocess=True, method=None,
                     datafind_error='raise', **fftparams):

    channel = get_channel(channel)

    if method is None:
        methods = set([key.split(',', 1)[1] for key in globalv.SPECTROGRAMS
                       if key.startswith('%s,' % channel.ndsname)])
        if len(methods) == 1:
            method = list(methods)[0]
        else:
            method = 'median-mean'
    if format in ['rayleigh']:
        method = format

    # clean fftparams dict using channel default values
    fftparams = _clean_fftparams(fftparams, channel)

    # key used to store the coherence spectrogram in globalv
    key = _make_key(channel, fftparams, method=method)

    # extract spectrogram stride from dict
    stride = float(fftparams.pop('stride'))

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
        nproc = multiprocess

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
                                        multiprocess=nproc, query=query,
                                        datafind_error=datafind_error, nds=nds)
        # calculate spectrograms
        if len(timeserieslist):
            vprint("    Calculating spectrograms for %s" % str(channel))
        for ts in timeserieslist:
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
                specgram = (specgram ** (1/2.)).filter(*filter_,
                                                       inplace=True) ** 2
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


@use_segmentlist
def get_coherence_spectrogram(channel_pair, segments, config=ConfigParser(),
                              cache=None, query=True, nds='guess',
                              return_=True, frametype=None, multiprocess=True,
                              datafind_error='raise', return_components=False,
                              **fftparams):
    """Retrieve the time-series and generate a coherence spectrogram of
    the two given channels
    """

    specs = _get_coherence_spectrogram(channel_pair, segments,
                                       config=config, cache=cache,
                                       query=query, nds=nds,
                                       return_=return_, frametype=frametype,
                                       multiprocess=multiprocess,
                                       datafind_error=datafind_error,
                                       return_components=return_components,
                                       **fftparams)

    return specs


@use_segmentlist
def _get_coherence_spectrogram(channel_pair, segments, config=ConfigParser(),
                               cache=None, query=True, nds='guess',
                               return_=True, frametype=None, multiprocess=True,
                               datafind_error='raise', return_components=False,
                               **fftparams):

    channel1 = get_channel(channel_pair[0])
    channel2 = get_channel(channel_pair[1])

    # clean fftparams dict using channel 1 default values
    fftparams = _clean_fftparams(fftparams, channel1)

    # key used to store the coherence spectrogram in globalv
    key = _make_key(channel_pair, fftparams)

    # check how much of the coherence spectrogram still needs to be calculated
    new = segments - globalv.COHERENCE_SPECTROGRAMS.get(
              key, SpectrogramList()).segments

    # get processes
    if multiprocess is True:
        nproc = count_free_cores()
    elif multiprocess is False:
        nproc = 1
    else:
        nproc = multiprocess

    # if there are no existing spectrogram, initialize as a list
    globalv.COHERENCE_SPECTROGRAMS.setdefault(key, SpectrogramList())

    # XXX HACK: use dummy timeseries to find lower sampling rate
    if len(segments) > 0:
        s = segments[0].start
        dts1 = get_timeseries(channel1, SegmentList([Segment(s, s+1)]),
                              config=config, cache=cache, frametype=frametype,
                              multiprocess=nproc, query=query,
                              datafind_error=datafind_error, nds=nds)
        dts2 = get_timeseries(channel2, SegmentList([Segment(s, s+1)]),
                              config=config, cache=cache, frametype=frametype,
                              multiprocess=nproc, query=query,
                              datafind_error=datafind_error, nds=nds)
        sampling = min(dts1[0].sample_rate.value, dts2[0].sample_rate.value)
    else:
        sampling = None;

    # keys used to store component spectrograms in globalv
    components = ('Cxy', 'Cxx', 'Cyy')
    ckeys = [_make_key([channel1, channel2], fftparams, sampling=sampling),
             _make_key(channel1, fftparams, sampling=sampling),
             _make_key(channel2, fftparams, sampling=sampling)]

    # initialize component lists if they don't exist yet
    for ck in ckeys:
        globalv.COHERENCE_COMPONENTS.setdefault(ck, SpectrogramList())

    # get data if query=True or if there are new segments
    query &= abs(new) != 0

    if query:

        # the intersecting segments will be calculated when needed
        intersection = None

        # loop over components needed to calculate coherence
        stride = float(fftparams.pop('stride'))
        for comp in components:

            # key used to store this component in globalv (incl sample rate)
            ckey = ckeys[components.index(comp)]

            try:
                filter_ = channel1.frequency_response
            except AttributeError:
                filter_ = None
            else:
                if isinstance(filter_, str):
                    filter_ = eval(filter_)

            # check how much of this component still needs to be calculated
            req = new - globalv.COHERENCE_COMPONENTS.get(
                            ckey, SpectrogramList()).segments

            # get data if there are new segments
            if abs(req) != 0:

                # get time-series data
                if stride is not None:
                    tmp = type(req)()
                    for s in req:
                        if abs(s) < stride:
                            continue
                        else:
                            d = float(abs(s))
                            tmp.append(type(s)(s[0], s[0]+d//stride * stride))
                    req = tmp

                # calculate intersection of segments lazily
                # this should occur on first pass (Cxy)
                if intersection is None:
                    total1 = get_timeseries(
                                 channel1, req, config=config,
                                 cache=cache, frametype=frametype,
                                 multiprocess=nproc, query=query,
                                 datafind_error=datafind_error,
                                 nds=nds)
                    total2 = get_timeseries(
                                 channel2, req, config=config,
                                 cache=cache, frametype=frametype,
                                 multiprocess=nproc, query=query,
                                 datafind_error=datafind_error,
                                 nds=nds)
                    intersection = total1.segments & total2.segments

                # get required timeseries data (using intersection)
                tslist1, tslist2 = [], []
                if comp in ('Cxy', 'Cxx'):
                    tslist1 = get_timeseries(
                                  channel1, intersection, config=config,
                                  cache=cache, frametype=frametype,
                                  multiprocess=nproc, query=query,
                                  datafind_error=datafind_error,
                                  nds=nds)
                if comp in ('Cxy', 'Cyy'):
                    tslist2 = get_timeseries(
                                  channel2, intersection, config=config,
                                  cache=cache, frametype=frametype,
                                  multiprocess=nproc, query=query,
                                  datafind_error=datafind_error,
                                  nds=nds)

                # calculate component
                if len(tslist1) + len(tslist2):
                    vprint("    Calculating component %s for coherence "
                           "spectrogram for %s and %s @ %d Hz" % (
                               comp, str(channel1), str(channel2), sampling))

                for ts1, ts2 in izip_longest(tslist1, tslist2):

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
                        specgram = ts1.spectrogram(stride, nproc=nproc,
                                                   cross=ts2, **fftparams)
                    elif comp == 'Cxx':
                        specgram = ts1.spectrogram(stride, nproc=nproc,
                                                   **fftparams)
                    elif comp == 'Cyy':
                        specgram = ts2.spectrogram(stride, nproc=nproc,
                                                   **fftparams)

                    if filter_:
                        specgram = (specgram**(1/2.)).filter(*filter_,
                                                             inplace=True) ** 2

                    globalv.COHERENCE_COMPONENTS[ckey].append(specgram)
                    globalv.COHERENCE_COMPONENTS[ckey].coalesce()

                    vprint('.')

                if len(tslist1) + len(tslist2):
                    vprint('\n')

    # calculate coherence from the components and store in globalv
    for (cxy, cxx, cyy) in izip(
            *[globalv.COHERENCE_COMPONENTS[ck] for ck in ckeys]):
        csg = abs(cxy)**2 / cxx / cyy
        globalv.COHERENCE_SPECTROGRAMS[key].append(csg)
        globalv.COHERENCE_SPECTROGRAMS[key].coalesce()

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
        for specgram in globalv.COHERENCE_SPECTROGRAMS[key]:
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
            globalv.SPECTRUM[name] = FrequencySeries([], channel=channel, f0=0,
                                                     df=1, unit=units.Unit(''))
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


def get_coherence_spectrum(channel_pair, segments, config=ConfigParser(),
                           cache=None, query=True, nds='guess', return_=True,
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
                        spec.unit = cxxlist[0].unit
                    cdict[comp] = speclist[index].join(gap='ignore')
                else:
                    raise

        # average spectrograms to get PSDs and CSD
        try:
            Cxy = cdict['Cxy'].percentile(50)
            Cxx = cdict['Cxx'].percentile(50)
            Cyy = cdict['Cyy'].percentile(50)
            globalv.COHERENCE_SPECTRUM[name] = FrequencySeries(
                abs(Cxy)**2 / Cxx / Cyy, f0=Cxx.f0, df=Cxx.df)
        except (ValueError, IndexError):
            globalv.COHERENCE_SPECTRUM[name] = FrequencySeries(
                [], channel=channel1, f0=0, df=1, unit=units.Unit(''))
            globalv.COHERENCE_SPECTRUM[cmin] = globalv.COHERENCE_SPECTRUM[name]
            globalv.COHERENCE_SPECTRUM[cmax] = globalv.COHERENCE_SPECTRUM[name]
        else:
            # FIXME: how to calculate percentiles correctly?
            globalv.COHERENCE_SPECTRUM[cmin] = FrequencySeries(
                abs(cdict['Cxy'].percentile(5))**2 /
                    cdict['Cxx'].percentile(95) /
                    cdict['Cyy'].percentile(95), f0=Cxx.f0, df=Cxx.df)
            globalv.COHERENCE_SPECTRUM[cmax] = FrequencySeries(
                abs(cdict['Cxy'].percentile(95))**2 /
                    cdict['Cxx'].percentile(5) /
                    cdict['Cyy'].percentile(5), f0=Cxx.f0, df=Cxx.df)

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


def add_timeseries(timeseries, key=None, coalesce=True):
    """Add a `TimeSeries` to the global memory cache
    """
    update_missing_channel_params(timeseries.channel)
    if key is None:
        key = timeseries.name or timeseries.channel.ndsname
    if isinstance(timeseries, StateVector):
        globalv.DATA.setdefault(key, StateVectorList())
    else:
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


def add_coherence_spectrogram(specgram, key=None, coalesce=True):
    """Add a `Coherence spectrogram` to the global memory cache
    """
    if key is None:
        key = specgram.name or str(specgram.channel)
    globalv.COHERENCE_SPECTROGRAMS.setdefault(key, SpectrogramList())
    globalv.COHERENCE_SPECTROGRAMS[key].append(specgram)
    if coalesce:
        globalv.COHERENCE_SPECTROGRAMS[key].coalesce()


@use_segmentlist
def get_spectrograms(channels, segments, config=ConfigParser(), cache=None,
                     query=True, nds='guess', format='power', return_=True,
                     method='median-mean', frametype=None, multiprocess=True,
                     datafind_error='raise', **fftparams):
    """Get spectrograms for multiple channels
    """
    channels = map(get_channel, channels)
    # get timeseries data in bulk
    if query:
        qchannels = map(
            get_channel, set(c2 for c in
                map(lambda x: re_channel.findall(x.ndsname), channels)
                for c2 in c))
        if format in ['rayleigh']:
            method_ = format
        else:
            method_ = method
        keys = []
        for channel in qchannels:
            fftparams_ = _clean_fftparams(fftparams, channel)
            keys.append(_make_key(channel, fftparams, method=method))
        havesegs = reduce(operator.and_, (globalv.SPECTROGRAMS.get(
            key, SpectrogramList()).segments for key in keys))
        new = segments - havesegs
        strides = set([getattr(c, 'stride', 0) for c in qchannels])
        if len(strides) == 1:
            stride = strides.pop()
            new = type(new)([s for s in new if abs(s) >= stride])
        get_timeseries_dict(qchannels, new, config=config, cache=cache,
                            multiprocess=multiprocess, frametype=frametype,
                            datafind_error=datafind_error, nds=nds,
                            return_=False)
    # loop over channels and generate spectrograms
    out = OrderedDict()
    for channel in channels:
        out[channel] = get_spectrogram(
            channel, segments, config=config, cache=cache, query=query,
            nds=nds, format=format, multiprocess=multiprocess,
            return_=return_, method=method, datafind_error=datafind_error,
            **fftparams)
    return out


@use_segmentlist
def get_coherence_spectrograms(channel_pairs, segments, config=ConfigParser(),
                               cache=None, query=True, nds='guess',
                               return_=True, frametype=None, multiprocess=True,
                               datafind_error='raise', **fftparams):
    """Get coherence spectrograms for multiple channels
    """

    channels = map(get_channel, channel_pairs)

    # get timeseries data in bulk
    if query:
        qchannels = map(
            get_channel, set(
                c2 for c in
                map(lambda x: re_channel.findall(x.ndsname), channels)
                for c2 in c))
        keys = [channel.ndsname for channel in qchannels]
        havesegs = reduce(operator.and_, (globalv.COHERENCE_SPECTROGRAMS.get(
            key, SpectrogramList()).segments for key in keys))
        new = segments - havesegs
        strides = set([getattr(c, 'stride', 0) for c in qchannels])
        if len(strides) == 1:
            stride = strides.pop()
            new = type(new)([s for s in new if abs(s) >= stride])
        get_timeseries_dict(qchannels, new, config=config, cache=cache,
                            multiprocess=multiprocess, frametype=frametype,
                            datafind_error=datafind_error, nds=nds,
                            return_=False)
    # loop over channels and generate spectrograms
    out = OrderedDict()

    for channel_pair in zip(channels[0::2], channels[1::2]):
        out[channel] = get_coherence_spectrogram(
                            channel_pair, segments, config=config, cache=cache,
                            query=query, nds=nds, multiprocess=multiprocess,
                            return_=return_, datafind_error=datafind_error,
                            **fftparams)
    return out


def get_range_channel(channel, **rangekwargs):
    """Return the meta-channel name used to store range data
    """
    if not rangekwargs:
        rangekwargs = {'mass1': 1.4, 'mass2': 1.4}
    rangetype = 'energy' in rangekwargs and 'burst' or 'inspiral'
    re_float = re.compile('[.-]')
    rkey = '_'.join(['%s_%s' % (re_cchar.sub('_', key),
                                re_float.sub('_', str(val))) for key, val in
                    rangekwargs.iteritems()])
    channel = get_channel(channel)
    return '%s_%s' % (channel.ndsname, rkey)


@use_segmentlist
def get_range(channel, segments, config=ConfigParser(), cache=None,
              query=True, nds='guess', return_=True, multiprocess=True,
              datafind_error='raise', frametype=None,
              stride=None, fftlength=None, overlap=None,
              method=None, **rangekwargs):
    """Calculate the sensitive distance for a given strain channel
    """
    if not rangekwargs:
        rangekwargs = {'mass1': 1.4, 'mass2': 1.4}
    rangetype = 'energy' in rangekwargs and 'burst' or 'inspiral'
    if rangetype == 'burst':
        range_func = astro.burst_range
    else:
        range_func = astro.inspiral_range
    channel = get_channel(channel)
    key = get_range_channel(channel, **rangekwargs)
    # get old segments
    havesegs = globalv.DATA.get(key, TimeSeriesList()).segments
    new = segments - havesegs
    query &= abs(new) != 0
    # calculate new range
    out = TimeSeriesList()
    if query:
        # get spectrograms
        spectrograms = get_spectrogram(channel, new, config=config,
                                       cache=cache, multiprocess=multiprocess,
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


def _make_key(channels, fftparams, method=None, sampling=None):

    # generate the key string used to store spectrograms in `globalv`

    # foundation is channel name(s), input is a single channel or a list
    if isinstance(channels, basestring):
        key = channels
    elif isinstance(channels, Channel):
        key = channels.ndsname
    elif all(isinstance(chan, basestring) for chan in channels):
        key = ','.join(channels)
    elif all(isinstance(chan, Channel) for chan in channels):
        key = ','.join(map(lambda chan: chan.ndsname, channels))
    else:
        raise TypeError('Input is not a channel or list of channels.')

    # parameters to append to channel names
    p = []
    p.append(method)
    p.append(fftparams.get('window', None))
    for param in ['stride', 'fftlength', 'overlap']:
        try:
            p.append(float(fftparams[param]))
        except KeyError:
            p.append(None)
    p.append(sampling)

    # build string
    key = key + ';' + ';'.join(map(str, p))

    return key


def _clean_fftparams(fftparams, channel, **defaults):

    channel = get_channel(channel)
    # replace missing fftparams with defaults or values from given channel

    # set defaults
    out = {
        'window': None,
        'fftlength': 1,
        'overlap': .5,
        'method': 'median-mean'
    }
    out.update(defaults)
    out.update(fftparams)
    if out['overlap'] != 0:
        out.setdefault('stride', ceil(out['fftlength'] * 1.5))
    else:
        out.setdefault('stride', out['fftlength'])


    # channel values take precedence
    for param in ['window', 'method']:
        if hasattr(channel, param):
            out[param] = getattr(channel, param)
    for param in ['fftlength', 'overlap', 'stride']:
        if hasattr(channel, param):
            out[param] = float(getattr(channel, param))
        # if channel doesn't have parameter, set it at the earliest oppotunity
        else:
            setattr(channel, param, out[param])

    # checks
    if out['stride'] == 0:
        raise ZeroDivisionError("Spectrogram stride is 0")
    elif out['fftlength'] == 0:
        raise ZeroDivisionError("FFT length is 0")

    return out
