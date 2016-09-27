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
import warnings
from math import (floor, ceil)
from time import sleep
try:
    from collections import OrderedDict
except ImportError:
    from ordereddict import OrderedDict

import nds2

from astropy import units

from glue import datafind
from glue.lal import Cache

from gwpy.segments import SegmentList
from gwpy.timeseries import (TimeSeriesList, TimeSeriesDict,
                             StateVector, StateVectorDict)
try:
    from gwpy.timeseries import StateVectorList
except ImportError:
    StateVectorList = TimeSeriesList
from gwpy.io import nds as ndsio

from .. import globalv
from ..utils import (vprint, count_free_cores)
from ..config import (GWSummConfigParser, NoSectionError, NoOptionError)
from ..channels import (get_channel, update_missing_channel_params,
                        split_combination as split_channel_combination)
from .utils import (use_configparser, use_segmentlist, make_globalv_key)
from .mathutils import get_with_math


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
    'low-latency h(t)': re.compile('[A-Z][0-9]_DMT_C00\Z'),
    'calibrated h(t) version 0': re.compile('[A-Z][0-9]_HOFT_C00\Z'),
    'calibrated h(t) version 1': re.compile('[A-Z][0-9]_HOFT_C01\Z'),
    'calibrated h(t) version 2': re.compile('[A-Z][0-9]_HOFT_C02\Z'),
    'DMT SenseMon on GDS h(t)': re.compile('SenseMonitor_hoft_[A-Z][0-9]_M\Z'),
    'DMT SenseMon on front-end h(t)': re.compile(
        'SenseMonitor_CAL_[A-Z][0-9]_M\Z'),
}

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


# -- utilities ----------------------------------------------------------------

@use_configparser
def find_frames(ifo, frametype, gpsstart, gpsend, config=GWSummConfigParser(),
                urltype='file', gaps='warn', onerror='raise'):
    """Query the datafind server for GWF files for the given type

    Parameters
    ----------
    ifo : `str`
        prefix for the IFO of interest (either one or two characters)

    frametype : `str`
        name of the frametype to find

    gpsstart : `int`
        GPS start time of the query

    gpsend : `int`
        GPS end time of the query

    config : `~ConfigParser.ConfigParser`, optional
        configuration with `[datafind]` section containing `server`
        specification, otherwise taken from the environment

    urltype : `str`, optional
        what type of file paths to return, default: `file`

    gaps : `str`, optional
        what to do when gaps are detected, one of

        - `ignore` : do nothing
        - `warn` : display the existence of gaps but carry on
        - `raise` : raise an exception

    onerror : `str`, optional
        what to do when the `~glue.datafind` query itself fails, same
        options as for ``gaps``

    Returns
    -------
    cache : `~glue.lal.Cache`
        a list of structured frame file descriptions matching the ifo and
        frametype requested
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

    # parse match
    try:
        frametype, match = frametype.split('|', 1)
    except ValueError:
        match = None

    def _query():
        if cert is not None:
            dfconn = datafind.GWDataFindHTTPSConnection(
                host=host, port=port, cert_file=cert, key_file=key)
        else:
            dfconn = datafind.GWDataFindHTTPConnection(host=host, port=port)
        return dfconn.find_frame_urls(ifo[0].upper(), frametype, gpsstart,
                                      gpsend, urltype=urltype, on_gaps=gaps,
                                      match=match)
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
    """Return the segments covered by one or more data caches

    Parameters
    ----------
    *cache : `~glue.lal.Cache`
        one or more file caches

    Returns
    -------
    segments : `~gwpy.segments.SegmentList`
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
    """Find the frametype associated with the given channel

    If the input channel has a `frametype` attribute, that will be used,
    otherwise the frametype will be guessed based on the channel name and
    any trend options given
    """
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


def get_channel_type(name):
    """Returns the probable type of this channel, based on the name

    Parameters
    ----------
    name : `str`
        the name of the channel

    Returns
    -------
    type : `str`
        one of ``'adc'``, ``'proc'``, or ``'sim'``
    """
    channel = get_channel(name)
    # find GEO h(t) channels
    if channel.ifo == 'G1' and channel.system in ['DER']:
        return 'proc'
    # find LIGO h(t) or LOSC channels
    if channel.ifo in ['L1', 'H1'] and channel.system in ['GDS', 'LOSC']:
        return 'proc'
    # default to ADC
    return 'adc'


# -- data accessors -----------------------------------------------------------

@use_configparser
@use_segmentlist
def get_timeseries_dict(channels, segments, config=GWSummConfigParser(),
                        cache=None, query=True, nds=None, multiprocess=True,
                        frametype=None, statevector=False, return_=True,
                        datafind_error='raise', **ioargs):
    """Retrieve the data for a set of channels

    Parameters
    ----------
    channels : `list` of `str` or `~gwpy.detector.Channel`
        the channels you want to get

    segments : `~gwpy.segments.SegmentList`
        the data segments of interest

    config : `~gwsumm.config.GWSummConfigParser`
        the configuration for this analysis

    query : `bool`, optional
        whether you want to retrieve new data from the source if it hasn't
        been loaded already

    nds : `bool`, optional
        whether to try and use NDS2 for data access, default is to guess
        based on other arguments and the environment

    multiprocess : `bool`, `int`, optional
        whether to use multiprocessing for file reading

    frametype : `str`, optional`
        the frametype of the target channels, if not given, this will be
        guessed based on the channel name(s)

    statevector : `bool`, optional
        whether you want to load `~gwpy.timeseries.StateVector` rather than
        `~gwpy.timeseries.TimeSeries` data

    datafind_error : `str`, optional
        what to do in the event of a datafind error, one of

        - 'raise' : stop immediately upon error
        - 'warn' : print warning and continue as if no frames had been found
        - 'ignore' : print nothing and continue with no frames

    return_ : `bool`, optional
        whether you actually want anything returned to you, or you are just
        calling this function to load data for use later

    **ioargs
        all other keyword arguments are passed to the relevant data
        reading method (either `~gwpy.timeseries.TimeSeriesDict.read` or
        `~gwpy.timeseries.TimeSeriesDict.fetch` or state-vector equivalents)

    Returns
    -------
    datalist : `dict` of `~gwpy.timeseries.TimeSeriesList`
        a set of `(channel, TimeSeriesList`) pairs
    """
    # separate channels by type
    if query:
        if frametype is not None or cache is not None:
            frametypes = {(None, frametype): channels}
        else:
            frametypes = dict()
            allchannels = set([c for group in
                map(split_channel_combination, channels) for c in group])
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
        for name in channels:
            channel = get_channel(name)
            out[channel.ndsname] = get_with_math(
                name, segments, _get_timeseries_dict, get_timeseries,
                config=config, query=False, statevector=statevector, **ioargs)
        return out


@use_segmentlist
def _get_timeseries_dict(channels, segments, config=None,
                         cache=None, query=True, nds=None, frametype=None,
                         multiprocess=True, return_=True, statevector=False,
                         archive=True, datafind_error='raise', **ioargs):
    """Internal method to retrieve the data for a set of like-typed
    channels using the :meth:`TimeSeriesDict.read` accessor.
    """
    channels = map(get_channel, channels)

    # set classes
    if statevector:
        ListClass = StateVectorList
        DictClass = StateVectorDict
    else:
        ListClass = TimeSeriesList
        DictClass = TimeSeriesDict

    # check we have a configparser
    if config is None:
        config = GWSummConfigParser()

    # read segments from global memory
    keys = dict((c.ndsname, make_globalv_key(c)) for c in channels)
    havesegs = reduce(operator.and_,
                      (globalv.DATA.get(keys[channel.ndsname],
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
    if nds is None and cache is not None:
        nds = False
    elif nds is None:
        nds = 'LIGO_DATAFIND_SERVER' not in os.environ

    # read new data
    query &= (abs(new) > 0)
    if cache is not None:
        query &= len(cache) > 0
    if query:
        for channel in channels:
            globalv.DATA.setdefault(keys[channel.ndsname], ListClass())
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
                cachesegments = find_cache_segments(fcache)
                gaps = SegmentList([span]) - cachesegments
                if len(fcache) == 0 and frametype == '%s_R' % ifo:
                    frametype = '%s_C' % ifo
                    vprint("    Moving to backup frametype %s\n" % frametype)
                    fcache = find_frames(ifo, frametype, span[0], span[1],
                                         config=config, gaps='ignore',
                                         onerror=datafind_error)
                elif abs(gaps) and frametype == '%s_HOFT_C00' % ifo:
                    frametype = '%s_DMT_C00' % ifo
                    vprint("    Gaps discovered in aggregated h(t) type "
                           "%s_HOFT_C00, checking %s\n" % (ifo, frametype))
                    c2 = find_frames(ifo, frametype, span[0], span[1],
                                     config=config, gaps='ignore',
                                     onerror=datafind_error)
                    g2 = SegmentList([span]) - find_cache_segments(c2)
                    if abs(g2) < abs(gaps):
                        vprint("    Greater coverage with frametype %s\n"
                               % frametype)
                        fcache = c2
                    else:
                        vprint("    No extra coverage with frametype %s\n"
                               % frametype)

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
                    ctype.add(get_channel_type(channel))
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
                if (re.match('(?:(.*)_)?[A-Z]\d_M', str(frametype)) or
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
                key = keys[channel.ndsname]
                if (key in globalv.DATA and
                        data.span in globalv.DATA[key].segments):
                    continue
                if data.unit is None:
                    data.unit = 'undef'
                for seg in globalv.DATA[key].segments:
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
                        try:
                            data = filt(data)
                        except TypeError as e:
                            if 'Can only apply' in str(e):
                                data.value[:] = filt(data.value)
                            else:
                                raise
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
                if len(globalv.DATA[key]):
                    data._unit = globalv.DATA[key][-1].unit
                # update channel type for trends
                if (data.channel.type is None and
                       data.channel.trend is not None):
                    if data.dt.to('s').value == 1:
                        data.channel.type = 's-trend'
                    elif data.dt.to('s').value == 60:
                        data.channel.type = 'm-trend'
                # append and coalesce
                add_timeseries(data, key=key, coalesce=True)
            if multiprocess:
                vprint('.')
        if len(new):
            vprint("\n")

    if not return_:
        return

    # return correct data
    out = OrderedDict()
    for channel in channels:
        data = ListClass()
        if keys[channel.ndsname] not in globalv.DATA:
            out[channel.ndsname] = ListClass()
        else:
            for ts in globalv.DATA[keys[channel.ndsname]]:
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
def get_timeseries(channel, segments, config=None, cache=None,
                   query=True, nds=None, multiprocess=True,
                   frametype=None, statevector=False, return_=True,
                   datafind_error='raise', **ioargs):
    """Retrieve data for channel

    Parameters
    ----------
    channel : `str` or `~gwpy.detector.Channel`
        the name of the channel you want

    segments : `~gwpy.segments.SegmentList`
        the data segments of interest

    config : `~gwsumm.config.GWSummConfigParser`
        the configuration for this analysis

    cache : `~glue.lal.Cache`
        a cache of data files from which to read

    query : `bool`, optional
        whether you want to retrieve new data from the source if it hasn't
        been loaded already

    nds : `bool`, optional
        whether to try and use NDS2 for data access, default is to guess
        based on other arguments and the environment

    multiprocess : `bool`, `int`, optional
        whether to use multiprocessing for file reading

    frametype : `str`, optional`
        the frametype of the target channels, if not given, this will be
        guessed based on the channel name(s)

    statevector : `bool`, optional
        whether you want to load `~gwpy.timeseries.StateVector` rather than
        `~gwpy.timeseries.TimeSeries` data

    datafind_error : `str`, optional
        what to do in the event of a datafind error, one of

        - 'raise' : stop immediately upon error
        - 'warn' : print warning and continue as if no frames had been found
        - 'ignore' : print nothing and continue with no frames

    return_ : `bool`, optional
        whether you actually want anything returned to you, or you are just
        calling this function to load data for use later

    **ioargs
        all other keyword arguments are passed to the relevant data
        reading method (either `~gwpy.timeseries.TimeSeries.read` or
        `~gwpy.timeseries.TimeSeries.fetch` or state-vector equivalents)

    Returns
    -------
    data : `~gwpy.timeseries.TimeSeriesList`
        a list of `TimeSeries`
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


def add_timeseries(timeseries, key=None, coalesce=True):
    """Add a `TimeSeries` to the global memory cache

    Parameters
    ----------
    timeseries : `~gwpy.timeseries.TimeSeries` or `~gwpy.timeseries.StateVector`
        the data series to add

    key : `str`, optional
        the key with which to store these data, defaults to the
        `~gwpy.timeseries.TimeSeries.name` of the series

    coalesce : `bool`, optional
        coalesce contiguous series after adding, defaults to `True`
    """
    if timeseries.channel is not None:
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
