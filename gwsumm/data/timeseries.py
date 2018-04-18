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

from six.moves import reduce

from astropy import units

from lal.utils import CacheEntry

from glue import (datafind, lal as glue_lal)

from gwpy.io import nds2 as io_nds2
from gwpy.io.cache import cache_segments
from gwpy.segments import (Segment, SegmentList)
from gwpy.timeseries import (TimeSeriesList, TimeSeriesDict,
                             StateVector, StateVectorList, StateVectorDict)
from gwpy.timeseries.io.gwf import get_default_gwf_api

from .. import globalv
from ..utils import vprint
from ..config import (GWSummConfigParser, NoSectionError, NoOptionError)
from ..channels import (get_channel, update_missing_channel_params,
                        split_combination as split_channel_combination)
from .utils import (use_configparser, use_segmentlist, make_globalv_key)
from .mathutils import get_with_math

warnings.filterwarnings("ignore", "LAL has no unit corresponding")

# avoid DeprecationWarnings from glue re: CacheEntry
Cache = glue_lal.Cache
glue_lal.CacheEntry = Cache.entry_class = CacheEntry

OPERATOR = {
    '*': operator.mul,
    '-': operator.sub,
    '+': operator.add,
    '/': operator.truediv,
}

FRAMETYPE_REGEX = {
    'commissioning': re.compile('[A-Z][0-9]_C\Z'),
    'raw data': re.compile('([A-Z][0-9]_R\Z|raw)'),
    'second-trend': re.compile('[A-Z][0-9]_T\Z'),
    'minute-trend': re.compile('[A-Z][0-9]_M\Z'),
    'low-latency h(t)': re.compile('([A-Z][0-9]_DMT_C00\Z|[A-Z][0-9]_llhoft)'),
    'calibrated h(t) version 0': re.compile('[A-Z][0-9]_HOFT_C00\Z'),
    'calibrated h(t) version 1': re.compile(
        '([A-Z][0-9]_HOFT_C01|G1_RDS_C01_L3)\Z'),
    'calibrated h(t) version 2': re.compile('[A-Z][0-9]_HOFT_C02\Z'),
    'DMT SenseMon on GDS h(t)': re.compile('SenseMonitor_hoft_[A-Z][0-9]_M\Z'),
    'DMT SenseMon on front-end h(t)': re.compile(
        'SenseMonitor_CAL_[A-Z][0-9]_M\Z'),
}

# list of GWF frametypes that contain only ADC channels
#     allows big memory/time savings when reading with frameCPP
try:
    GWF_API = get_default_gwf_api()
except ImportError:
    GWF_API = None

# frameCPP I/O optimisations
ADC_TYPES = [
    'R', 'C',  # old LIGO raw and commissioning types
    'T', 'M',  # old LIGO trend types
    'H1_R', 'H1_C', 'L1_R', 'L1_C',  # new LIGO raw and commissioning types
    'H1_T', 'H1_M', 'L1_T', 'L1_M',  # new LIGO trend types
    'raw',  # Virgo raw type
]

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


# -- utilities ----------------------------------------------------------------

re_gwf_gps_epoch = re.compile('[-\/](?P<gpsepoch>\d+)$')


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

    # extend cache beyond datafind's knowledge to reduce latency
    try:
        latest = cache[-1]
        ngps = len(re_gwf_gps_epoch.search(
            os.path.dirname(latest.path)).groupdict()['gpsepoch'])
    except (IndexError, AttributeError):
        pass
    else:
        while True:
            s, e = latest.segment
            if s >= gpsend:
                break
            # replace GPS time of file basename
            new = latest.path.replace('-%d-' % s, '-%d-' % e)
            # replace GPS epoch in dirname
            new = new.replace('%s/' % str(s)[:ngps], '%s/' % str(e)[:ngps])
            if os.path.isfile(new):
                latest = CacheEntry.from_T050017(new)
                cache.append(latest)
            else:
                break

    # validate files existing and return
    cache, _ = cache.checkfilesexist()
    vprint(' %d found.\n' % len(cache))
    return cache


def find_best_frames(ifo, frametype, start, end, **kwargs):
    """Find frames for the given type, replacing with a better type if needed
    """
    # find cache for this frametype
    cache = find_frames(ifo, frametype, start, end, **kwargs)

    # check for gaps in current cache
    span = SegmentList([Segment(start, end)])
    gaps = span - cache_segments(cache)

    # if gaps and using aggregated h(t), check short files
    if abs(gaps) and frametype == '%s_HOFT_C00' % ifo:
        f2 = '%s_DMT_C00' % ifo
        vprint("    Gaps discovered in aggregated h(t) type "
               "%s, checking %s\n" % (frametype, f2))
        kwargs['gaps'] = 'ignore'
        c2 = find_frames(ifo, f2, start, end, **kwargs)
        g2 = span - cache_segments(c2)
        if abs(g2) < abs(gaps):
            vprint("    Greater coverage with frametype %s\n" % f2)
            return c2, f2
        vprint("    No extra coverage with frametype %s\n" % f2)

    return cache, frametype


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
            ndstype = io_nds2.Nds2ChannelType.find(channel.type)
        except (AttributeError, KeyError, ValueError):
            ndstype = channel.type = 'raw'
        if ndstype == io_nds2.Nds2ChannelType.MTREND:
            ftype = 'M'
        elif ndstype == io_nds2.Nds2ChannelType.STREND:
            ftype = 'T'
        elif ndstype == io_nds2.Nds2ChannelType.RDS:
            ftype = 'LDAS_C02_L2'
        elif ndstype == io_nds2.Nds2ChannelType.ONLINE:
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
                        cache=None, query=True, nds=None, nproc=1,
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

    nproc : `int`, optional
        number of parallel cores to use for file reading, default: ``1``

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
            allchannels = set([
                c for group in map(split_channel_combination, channels) for
                c in group])
            for channel in allchannels:
                channel = get_channel(channel)
                ifo = channel.ifo
                ftype = find_frame_type(channel)
                id_ = (ifo, ftype)
                if id_ in frametypes:
                    frametypes[id_].append(channel)
                else:
                    frametypes[id_] = [channel]
        for ftype, channellist in frametypes.items():
            _get_timeseries_dict(channellist, segments, config=config,
                                 cache=cache, query=query, nds=nds,
                                 nproc=nproc, frametype=ftype[1],
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
                         nproc=1, return_=True, statevector=False,
                         archive=True, datafind_error='raise', **ioargs):
    """Internal method to retrieve the data for a set of like-typed
    channels using the :meth:`TimeSeriesDict.read` accessor.
    """
    channels = list(map(get_channel, channels))

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

    # read channel information
    filter_ = dict()
    resample = dict()
    dtype_ = dict()
    for channel in channels:
        name = str(channel)
        try:
            filter_[name] = channel.filter
        except AttributeError:
            pass
        try:
            resample[name] = float(channel.resample)
        except AttributeError:
            pass
        if channel.dtype is None:
            dtype_[name] = ioargs.get('dtype')
        else:
            dtype_[name] = channel.dtype

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
            ndsconnection = io_nds2.connect(host, port)
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
            if cache is not None:
                fcache = cache.sieve(ifos=ifo[0], description=frametype,
                                     exact_match=True)
            else:
                fcache = Cache()
            if (cache is None or len(fcache) == 0) and len(new):
                span = new.extent().protract(8)
                fcache, frametype = find_best_frames(
                    ifo, frametype, span[0], span[1],
                    config=config, gaps='ignore', onerror=datafind_error)

            # parse discontiguous cache blocks and rebuild segment list
            new &= cache_segments(fcache)
            source = 'frames'

            if cache is None and GWF_API == 'framecpp':
                # set ctype if reading with framecpp (using datafind)
                if frametype in ADC_TYPES:
                    ioargs['type'] = 'adc'

        for channel in channels:
            channel.frametype = frametype

        # check whether each channel exists for all new times already
        qchannels = []
        qresample = {}
        qdtype = {}
        for channel in channels:
            name = str(channel)
            oldsegs = globalv.DATA.get(name, ListClass()).segments
            if abs(new - oldsegs) != 0:
                qchannels.append(name)
                if name in resample:
                    qresample[name] = resample[name]
                qdtype[name] = dtype_.get(name, ioargs.get('dtype'))
        ioargs['dtype'] = qdtype

        # loop through segments, recording data for each
        if len(new):
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
                if segment[1] == cache_segments(fcache)[-1][1] and qresample:
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
                    segcache = segcache.sieve(
                        segment=type(segment)(segstart, segend))
                else:
                    segstart, segend = map(float, segment)
                # read data
                tsd = DictClass.read(segcache, qchannels,
                                     start=segstart, end=segend,
                                     nproc=nproc, resample=qresample,
                                     **ioargs)

            for c, data in tsd.items():
                channel = get_channel(c)
                key = keys[channel.ndsname]
                if (key in globalv.DATA and
                        data.span in globalv.DATA[key].segments):
                    continue
                if data.unit is None:
                    data.unit = 'undef'
                for i, seg in enumerate(globalv.DATA[key].segments):
                    if seg in data.span:
                        # new data completely covers existing segment
                        # (and more), so just remove the old stuff
                        globalv.DATA[key].pop(i)
                        break
                    elif seg.intersects(data.span):
                        # new data extends existing segment, so only keep
                        # the really new stuff
                        data = data.crop(*(data.span - seg))
                        break
                try:
                    filt = filter_[str(channel)]
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
                if data.channel.type is None and (
                        data.channel.trend is not None):
                    if data.dt.to('s').value == 1:
                        data.channel.type = 's-trend'
                    elif data.dt.to('s').value == 60:
                        data.channel.type = 'm-trend'
                # append and coalesce
                add_timeseries(data, key=key, coalesce=True)
            if nproc > 1:
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
                   query=True, nds=None, nproc=1,
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

    nproc : `int`, optional
        number of parallel cores to use for file reading, default: ``1``

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
                              nproc=nproc, frametype=frametype,
                              statevector=statevector, return_=return_,
                              datafind_error=datafind_error, **ioargs)
    if return_:
        return out[channel.ndsname]
    return


def add_timeseries(timeseries, key=None, coalesce=True):
    """Add a `TimeSeries` to the global memory cache

    Parameters
    ----------
    timeseries : `TimeSeries` or `StateVector`
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
