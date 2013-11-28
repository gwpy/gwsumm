
"""Utilities for data handling and display
"""

from math import (floor, ceil)
try:
    from configparser import (ConfigParser, NoSectionError, NoOptionError)
except ImportError:
    from ConfigParser import (ConfigParser, NoSectionError, NoOptionError)

import numpy
import nds2

from glue import datafind
from glue.lal import Cache

from gwpy.segments import (DataQualityFlag, SegmentList)
from gwpy.timeseries import (TimeSeries, TimeSeriesList)
from gwpy.spectrogram import SpectrogramList

from . import (globalv, version)
from .utils import *

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__version__ = version.version

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
                   query=True, nds=False):
    """Retrieve the data (time-series) for a given channel
    """
    if isinstance(segments, DataQualityFlag):
        segments = segments.active
    # read segments from global memory
    havesegs = globalv.DATA.get(str(channel), TimeSeriesList()).segments
    new = segments - havesegs

    # read channel information
    filter_ = None
    if config.has_section(channel):
        channel = channel
    if config.has_section(channel):
        if config.has_option(channel, 'unit'):
            channel.unit = config.get(channel, 'unit')
        if config.has_option(channel, 'filter'):
            filter_ = eval(config.get(channel, 'filter'))

    # read new data
    globalv.DATA.setdefault(str(channel), TimeSeriesList())
    query &= (abs(new) != 0)
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
            fcache = cache.sieve(description=ftype, exact_match=True)
            if len(fcache) == 0 and len(new):
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
                   % (source, str(channel)))
        for segment in new:
            if nds:
                data = TimeSeries.fetch(channel, segment[0], segment[1],
                                        connection=ndsconnection,
                                        ndschanneltype=channel.type)
            else:
                segcache = fcache.sieve(segment=segment)
                data = TimeSeries.read(segcache, channel, segment[0],
                                       segment[1], verbose=globalv.VERBOSE)
            if not channel.sample_rate:
                channel.sample_rate = data.sample_rate
            if filter_:
                data = data.filter(*filter_)
            globalv.DATA[str(channel)].append(data)
            globalv.DATA[str(channel)].coalesce()
            vprint(".")
        vprint("\n")

    # return correct data
    out = TimeSeriesList()
    for seg in segments:
        for ts in globalv.DATA[str(channel)]:
            if ts.span.intersects(seg):
                out.append(ts.crop(*seg))
    return out.coalesce()


def get_spectrogram(channel, segments, config=ConfigParser(), cache=None,
                    query=True, nds=False, format='power', **fftparams):
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
        if hasattr(channel, 'stride'):
            stride = channel.stride
        elif 'stride' in fftparams:
            stride = fftparams.pop('stride', 0)
        else:
            stride = fftparams['fftlength']
        stride = float(stride)
        fftparams.setdefault('method', 'medianmean')
        for param in ['fftlength', 'fftstride']:
            if hasattr(channel, param):
                fftparams[param] = float(getattr(channel, param))
            elif param in fftparams:
                fftparams[param] = float(fftparams[param])
        if hasattr(channel, 'stride'):
            stride = channel.stride
        elif 'stride' in fftparams:
            stride = fftparams.pop('stride', 0)
        if stride:
            stride = float(stride)
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
            specgram = ts.spectrogram(stride, **fftparams)
            if filter_:
                specgram = (specgram ** (1/2.)).filter(*filter_, inplace=True) ** 2
            globalv.SPECTROGRAMS[str(channel)].append(specgram)
            globalv.SPECTROGRAMS[str(channel)].coalesce()
            vprint('.')
        vprint('\n')

    # return correct data
    out = SpectrogramList()
    for seg in segments:
        for specgram in globalv.SPECTROGRAMS[str(channel)]:
            if specgram.span.intersects(seg):
                if format in ['amplitude', 'asd']:
                    out.append(specgram.crop(*seg) ** (1/2.))
                else:
                    out.append(specgram.crop(*seg))
    return out.coalesce()


def get_spectrum(channel, segments, config=ConfigParser(), cache=None,
                 query=True, nds=False, format='power',
                 **fftparams):
    """Retrieve the time-series and generate a spectrogram of the given
    channel
    """
    if isinstance(segments, DataQualityFlag):
        segments = segments.active

    if not (str(channel) in globalv.SPECTRUM):
        speclist = get_spectrogram(channel, segments, config=config,
                                   cache=cache, query=query, nds=nds,
                                   **fftparams)
        specgram = speclist.join(gap='ignore')
        globalv.SPECTRUM[str(channel)] = specgram.percentile(50)
        cmin = '%s.min' % str(channel)
        globalv.SPECTRUM[cmin] = specgram.percentile(5)
        cmax = '%s.max' % str(channel)
        globalv.SPECTRUM[cmax] = specgram.percentile(95)

    cmin = '%s.min' % str(channel)
    cmax = '%s.max' % str(channel)
    out = (globalv.SPECTRUM[str(channel)], globalv.SPECTRUM[cmin],
           globalv.SPECTRUM[cmax])
    if format in ['amplitude', 'asd']:
        out = [s ** (1/2.) for s in out]
    return out
