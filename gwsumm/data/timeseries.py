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
from collections import OrderedDict
from configparser import (NoSectionError, NoOptionError)

from six.moves import reduce
from six.moves.urllib.parse import urlparse

from astropy import units

import gwdatafind

from gwpy.io import nds2 as io_nds2
from gwpy.io.cache import (file_segment, cache_segments)
from gwpy.io.gwf import data_segments
from gwpy.segments import (Segment, SegmentList)
from gwpy.timeseries import (TimeSeriesList, TimeSeriesDict,
                             StateVector, StateVectorList, StateVectorDict)
from gwpy.timeseries.io.gwf import get_default_gwf_api
from gwpy.utils.mp import multiprocess_with_queues

from .. import globalv
from ..utils import vprint
from ..config import GWSummConfigParser
from ..channels import (get_channel, update_missing_channel_params,
                        split_combination as split_channel_combination,
                        update_channel_params)
from .utils import (use_configparser, use_segmentlist, make_globalv_key)
from .mathutils import get_with_math

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

warnings.filterwarnings("ignore", "LAL has no unit corresponding")

OPERATOR = {
    '*': operator.mul,
    '-': operator.sub,
    '+': operator.add,
    '/': operator.truediv,
}

FRAMETYPE_REGEX = {
    'commissioning': re.compile(
        r'[A-Z][0-9]_C\Z'),
    'raw data': re.compile(
        r'([A-Z][0-9]_R\Z|raw)'),
    'second-trend': re.compile(
        r'[A-Z][0-9]_T\Z'),
    'minute-trend': re.compile(
        r'[A-Z][0-9]_M\Z'),
    'low-latency h(t)': re.compile(
        r'([A-Z][0-9]_DMT_C00\Z|[A-Z][0-9]_llhoft)'),
    'calibrated h(t) version 0': re.compile(
        r'[A-Z][0-9]_HOFT_C00\Z'),
    'calibrated h(t) version 1': re.compile(
        r'([A-Z][0-9]_HOFT_C01|G1_RDS_C01_L3)\Z'),
    'calibrated h(t) version 2': re.compile(
        r'[A-Z][0-9]_HOFT_C02\Z'),
    'DMT SenseMon on GDS h(t)': re.compile(
        r'SenseMonitor_hoft_[A-Z][0-9]_M\Z'),
    'DMT SenseMon on front-end h(t)': re.compile(
        r'SenseMonitor_CAL_[A-Z][0-9]_M\Z'),
}

# list of GWF frametypes that contain only ADC channels
#     allows big memory/time savings when reading with frameCPP
try:
    GWF_API = get_default_gwf_api()
except ImportError:
    GWF_API = None

# frameCPP I/O optimisations
ADC_TYPES = {
    'R', 'C',  # old LIGO raw and commissioning types
    'T', 'M',  # old LIGO trend types
    'H1_R', 'H1_C', 'L1_R', 'L1_C',  # new LIGO raw and commissioning types
    'H1_T', 'H1_M', 'L1_T', 'L1_M',  # new LIGO trend types
    'raw',  # Virgo raw type
}

SHORT_HOFT_TYPES = {  # map aggregated h(t) type to short h(t) type
    'H1_HOFT_C00': 'H1_DMT_C00',
    'L1_HOFT_C00': 'L1_DMT_C00',
    'V1Online': 'V1_llhoft',
}

VIRGO_HOFT_CHANNELS = {
    "V1:Hrec_hoft_16384Hz",
    "V1:DQ_ANALYSIS_STATE_VECTOR",
}


# -- utilities ----------------------------------------------------------------

re_gwf_gps_epoch = re.compile(r'[-\/](?P<gpsepoch>\d+)$')


def _urlpath(url):
    return urlparse(url).path


def sieve_cache(cache, ifo=None, tag=None, segment=None):
    def _sieve(url):
        try:
            uifo, utag, useg = gwdatafind.utils.filename_metadata(url)
        except (AttributeError, TypeError):  # CacheEntry
            uifo = url.observatory
            utag = url.description
            useg = url.segment
        if segment is not None:  # cast type to prevent TypeErrors
            useg = type(segment)(*useg)
        return (
            (ifo is None or uifo == ifo) and
            (tag is None or utag == tag) and
            (segment is None or useg.intersects(segment))
        )
    return type(cache)(filter(_sieve, cache))


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
        what to do when the `gwdatafind` query itself fails, same
        options as for ``gaps``

    Returns
    -------
    cache : `list` of `str`
        a list of file paths pointing at GWF files matching the request
    """
    vprint('    Finding %s-%s frames for [%d, %d)...'
           % (ifo[0], frametype, int(gpsstart), int(gpsend)))
    # find datafind host:port
    try:
        host = config.get('datafind', 'server')
    except (NoOptionError, NoSectionError):
        host = None
        port = None
    else:
        port = config.getint('datafind', 'port')

    # XXX HACK: LLO changed frame types on Dec 6 2013:
    LLOCHANGE = 1070291904
    if re.match(r'L1_{CRMT}', frametype) and gpsstart < LLOCHANGE:
        frametype = frametype[-1]

    # query frames
    ifo = ifo[0].upper()
    gpsstart = int(floor(gpsstart))
    gpsend = int(ceil(min(globalv.NOW, gpsend)))
    if gpsend <= gpsstart:
        return []

    # parse match
    try:
        frametype, match = frametype.split('|', 1)
    except ValueError:
        match = None

    def _query():
        return gwdatafind.find_urls(ifo[0].upper(), frametype, gpsstart,
                                    gpsend, urltype=urltype, on_gaps=gaps,
                                    match=match, host=host, port=port)
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
            cache = []

    # XXX: if querying for day of LLO frame type change, do both
    if (ifo[0].upper() == 'L' and frametype in ['C', 'R', 'M', 'T'] and
            gpsstart < LLOCHANGE < gpsend):
        start = len(cache) and cache[-1].segment[1] or gpsstart
        if start < gpsend:
            cache.extend(gwdatafind.find_urls(
                ifo[0].upper(), 'L1_%s' % frametype, start, gpsend,
                urltype=urltype, on_gaps=gaps, host=host, port=port)[1:])

    # extend cache beyond datafind's knowledge to reduce latency
    try:
        latest = cache[-1]
        ngps = len(re_gwf_gps_epoch.search(
            os.path.dirname(latest)).groupdict()['gpsepoch'])
    except (IndexError, AttributeError):
        pass
    else:
        while True:
            s, e = file_segment(latest)
            if s >= gpsend:
                break
            # replace GPS time of file basename
            new = latest.replace('-%d-' % s, '-%d-' % e)
            # replace GPS epoch in dirname
            new = new.replace('%s/' % str(s)[:ngps], '%s/' % str(e)[:ngps])
            if os.path.isfile(new):
                cache.append(new)
            else:
                break

    # validate files existing and return
    cache = list(filter(os.path.exists, map(_urlpath, cache)))
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
    if abs(gaps) and frametype in SHORT_HOFT_TYPES:
        f2 = SHORT_HOFT_TYPES[frametype]
        vprint("    Gaps discovered in aggregated h(t) type "
               "%s, checking %s\n" % (frametype, f2))
        kwargs['gaps'] = 'ignore'
        cache.extend(filter(lambda e: file_segment(e) in gaps,
                            find_frames(ifo, f2, start, end, **kwargs)))
        new = int(abs(gaps - cache_segments(cache)))
        if new:
            vprint("    %ss extra coverage with frametype %s\n" % (new, f2))
        else:
            vprint("    No extra coverage with frametype %s\n" % f2)

    return cache, frametype


def find_frame_type(channel):
    """Find the frametype associated with the given channel

    If the input channel has a `frametype` attribute, that will be used,
    otherwise the frametype will be guessed based on the channel name and
    any trend options given
    """
    if (channel.frametype is None and
            channel.type is None and
            channel.trend is not None):
        raise ValueError("Cannot automatically determine frametype for {0}, "
                         "please manually select frametype or give NDS-style "
                         "channel suffix".format(channel.name))
    if channel.frametype is None:
        try:
            ndstype = io_nds2.Nds2ChannelType.find(channel.type)
        except (AttributeError, KeyError, ValueError):
            ndstype = None
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


def frame_trend_type(ifo, frametype):
    """Returns the trend type of based on the given frametype
    """
    if ifo == 'C1' and frametype == 'M':
        return 'minute'
    if re.match(r'(?:(.*)_)?[A-Z]\d_M', str(frametype)):
        return 'minute'
    if re.match(r'(?:(.*)_)?[A-Z]\d_T', str(frametype)):
        return 'second'
    return None


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


def exclude_short_trend_segments(segments, ifo, frametype):
    """Remove segments from a list shorter than 1 trend sample
    """
    frametype = frametype or ''
    trend = frame_trend_type(ifo, frametype)
    if trend == 'minute':
        mindur = 60.
    elif trend == 'second':
        mindur = 1.
    else:
        mindur = 0.

    return type(segments)([s for s in segments if abs(s) >= mindur])


def all_adc(cache):
    """Returns `True` if all cache entries point to GWF file known to
    contain only ADC channels

    This is useful to set `type='adc'` when reading with frameCPP, which
    can greatly speed things up.
    """
    for path in cache:
        try:
            tag = os.path.basename(path).split('-')[1]
        except (AttributeError, TypeError):  # CacheEntry
            tag = path.description
            path = path.path
        if not path.endswith('.gwf') or tag not in ADC_TYPES:
            return False
    return True


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
                         archive=True, datafind_error='raise', dtype=None,
                         **ioargs):
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
        if channel.dtype is not None:
            dtype_[name] = channel.dtype
        elif dtype is not None:
            dtype_[name] = dtype

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

        ifo = channels[0].ifo

        # open NDS connection
        if nds:
            if config.has_option('nds', 'host'):
                host = config.get('nds', 'host')
                port = config.getint('nds', 'port')
                ndsconnection = io_nds2.connect(host, port)
            else:
                ndsconnection = None
            frametype = source = 'nds'
            ndstype = channels[0].type

            # get NDS channel segments
            if ndsconnection is not None and ndsconnection.get_protocol() > 1:
                span = map(int, new.extent())
                avail = io_nds2.get_availability(
                    channels, *span, connection=ndsconnection
                )
                new &= avail.intersection(avail.keys())

        # or find frame type and check cache
        else:
            frametype = frametype or channels[0].frametype
            new = exclude_short_trend_segments(new, ifo, frametype)

            if cache is not None:
                fcache = sieve_cache(cache, ifo=ifo[0], tag=frametype)
            else:
                fcache = []

            if (cache is None or len(fcache) == 0) and len(new):
                span = new.extent()
                fcache, frametype = find_best_frames(
                    ifo, frametype, span[0], span[1],
                    config=config, gaps='ignore', onerror=datafind_error)

            # parse discontiguous cache blocks and rebuild segment list
            new &= cache_segments(fcache)
            source = 'files'

            # if reading Virgo h(t) GWF data, filter out files that don't
            # contain the channel (Virgo state-vector only)
            _names = set(map(str, channels))
            _virgohoft = _names.intersection(VIRGO_HOFT_CHANNELS)
            if _virgohoft:
                vprint("    Determining available segments for "
                       "Virgo h(t) data...")
                new &= data_segments(fcache, _virgohoft.pop())

            # set channel type if reading with frameCPP
            if fcache and all_adc(fcache):
                ioargs['type'] = 'adc'

        # store frametype for display in Channel Information tables
        for channel in channels:
            channel.frametype = frametype

        # check whether each channel exists for all new times already
        qchannels = []
        for channel in channels:
            oldsegs = globalv.DATA.get(keys[channel.ndsname],
                                       ListClass()).segments
            if abs(new - oldsegs) != 0 and nds:
                qchannels.append(channel.ndsname)
            elif abs(new - oldsegs) != 0:
                qchannels.append(str(channel))

        # loop through segments, recording data for each
        if len(new):
            vprint("    Fetching data (from %s) for %d channels [%s]:\n"
                   % (source, len(qchannels),
                      nds and ndstype or frametype or ''))
        vstr = "        [{0[0]}, {0[1]})"
        for segment in new:
            # force reading integer-precision segments
            segment = type(segment)(int(segment[0]), int(segment[1]))
            if abs(segment) < 1:
                continue

            # reset to minute trend sample times
            if frame_trend_type(ifo, frametype) == 'minute':
                segment = Segment(*io_nds2.minute_trend_times(*segment))
                if abs(segment) < 60:
                    continue

            if nds:  # fetch
                tsd = DictClass.fetch(qchannels, segment[0], segment[1],
                                      connection=ndsconnection, type=ndstype,
                                      verbose=vstr.format(segment), **ioargs)
            else:  # read
                # NOTE: this sieve explicitly casts our segment to
                #       ligo.segments.segment to prevent `TypeError` from
                #       a mismatch with ligo.segments.segment
                segcache = sieve_cache(fcache, segment=segment)
                segstart, segend = map(float, segment)
                tsd = DictClass.read(segcache, qchannels, start=segstart,
                                     end=segend, nproc=nproc,
                                     verbose=vstr.format(segment), **ioargs)

            vprint("        post-processing...\n")

            # apply type casting (copy=False means same type just returns)
            for chan, ts in tsd.items():
                tsd[chan] = ts.astype(dtype_.get(chan, ts.dtype),
                                      casting='unsafe', copy=False)

            # apply resampling
            tsd = resample_timeseries_dict(tsd, nproc=1, **resample)

            # post-process
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

                # filter
                try:
                    filt = filter_[str(channel)]
                except KeyError:
                    pass
                else:
                    data = filter_timeseries(data, filt)

                if isinstance(data, StateVector) or ':GRD-' in str(channel):
                    data.override_unit(units.dimensionless_unscaled)
                    if hasattr(channel, 'bits'):
                        data.bits = channel.bits
                elif data.unit is None:
                    data.override_unit(channel.unit)

                # update channel type for trends
                if data.channel.type is None and (
                        data.channel.trend is not None):
                    if data.dt.to('s').value == 1:
                        data.channel.type = 's-trend'
                    elif data.dt.to('s').value == 60:
                        data.channel.type = 'm-trend'

                # append and coalesce
                add_timeseries(data, key=key, coalesce=True)

    # rebuilt global channel list with new parameters
    update_channel_params()

    if not return_:
        return

    return locate_data(channels, segments, list_class=ListClass)


def locate_data(channels, segments, list_class=TimeSeriesList):
    """Find and return available (already loaded) data
    """
    keys = dict((c.ndsname, make_globalv_key(c)) for c in channels)

    # return correct data
    out = OrderedDict()
    for channel in channels:
        data = list_class()
        if keys[channel.ndsname] not in globalv.DATA:
            out[channel.ndsname] = list_class()
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

    cache : `~glue.lal.Cache` or `list` of `str`
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
        # transfer parameters from timeseries.channel to the globalv channel
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


def resample_timeseries_dict(tsd, nproc=1, **sampling_dict):
    """Resample a `TimeSeriesDict`

    Parameters
    ----------
    tsd : `~gwpy.timeseries.TimeSeriesDict`
        the input dict to resample

    nproc : `int`, optional
        the number of parallel processes to use

    **sampling_dict
        ``<name>=<sampling frequency>`` pairs defining new
        sampling frequencies for keys of ``tsd``

    Returns
    -------
    resampled : `~gwpy.timeseries.TimeSeriesDict`
        a new dict with the keys from ``tsd`` and resampled values, if
        that key was included in ``sampling_dict``, or the original value
    """
    # define resample function (must take single argument)
    def _resample(args):
        ts, fs = args
        if fs and units.Quantity(fs, "Hz") == ts.sample_rate:
            warnings.warn(
                "requested resample rate for {0} matches native rate ({1}), "
                "please update configuration".format(ts.name, ts.sample_rate),
                UserWarning,
            )
        elif fs:
            return ts.resample(fs, ftype='fir', window='hamming')
        return ts

    # group timeseries with new sampling frequencies
    inputs = ((ts, sampling_dict.get(name)) for name, ts in tsd.items())

    # apply resampling
    resampled = multiprocess_with_queues(nproc, _resample, inputs)

    # map back to original dict keys
    return dict(zip(list(tsd), resampled))


def filter_timeseries(ts, filt):
    """Filter a `TimeSeris` using a function or a ZPK definition.
    """
    # filter with function
    if callable(filt):
        try:
            return filt(ts)
        except TypeError as e:
            if 'Can only apply' in str(e):  # units error
                ts.value[:] = filt(ts.value)
                return ts
            else:
                raise

    # filter with gain
    else:
        return ts.filter(*filt)
