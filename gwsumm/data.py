
"""Utilities for data handling and display
"""

from .version import version as __version__
__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

try:
    from configparser import (ConfigParser, NoSectionError, NoOptionError)
except ImportError:
    from ConfigParser import (ConfigParser, NoSectionError, NoOptionError)

import numpy
import nds2

from gwpy.segments import DataQualityFlag
from gwpy.timeseries import (TimeSeries, TimeSeriesList)
from gwpy.spectrogram import SpectrogramList

import globalv
from .utils import *


def get_timeseries(channel, segments, config=ConfigParser(), cache=None,
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
    if query and nds:
        if config.has_option('nds', 'host'):
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
            ndsconnection = None
        vprint("    Fetching data for %s" % str(channel))
        type_ = channel.type
        for segment in segments:
            data = TimeSeries.fetch(channel, segment[0], segment[1],
                                    connection=ndsconnection)
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