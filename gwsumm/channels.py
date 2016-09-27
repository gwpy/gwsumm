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

"""Utilities for channel access
"""

import threading
import re
from Queue import Queue

try:
    from kerberos import GSSError
except ImportError:
    GSSError = None

try:
    from gwpy.io.nds import NDS2_CHANNEL_TYPE
except ImportError:
    NDS2_TYPES = [
        'm-trend', 'online', 'raw', 'rds', 'reduced',
        's-trend', 'static', 'test-pt',
    ]
else:
    NDS2_TYPES = NDS2_CHANNEL_TYPE.keys()

from gwpy.detector import Channel

from . import globalv
from .mode import Mode
from .utils import re_quote

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

CIS_URL = 'https://cis.ligo.org'

re_channel = re.compile(r'[A-Z]\d:[a-zA-Z0-9]+'  # core channel section L1:TEST
                         '(?:[-_][a-zA-Z0-9_]+)?'  # underscore-delimiter parts
                         '(?:\.[a-z]+)?'  # trend type
                         '(?:,[a-z-]+)?')  # NDS channel type


class ThreadChannelQuery(threading.Thread):
    """Threaded CIS `Channel` query.
    """
    def __init__(self, inqueue, outqueue, find_trend_source=False, timeout=5):
        threading.Thread.__init__(self)
        self.in_ = inqueue
        self.out = outqueue
        self.find_trends = find_trend_source
        self.timeout = timeout

    def run(self):
        i, channel = self.in_.get()
        self.in_.task_done()
        try:
            self.out.put((i, get_channel(
                channel, find_trend_source=self.find_trends,
                timeout=self.timeout)))
        except Exception as e:
            self.out.put(e)
        self.out.task_done()


def get_channel(channel, find_trend_source=True, timeout=5):
    """Define a new :class:`~gwpy.detector.channel.Channel`

    Parameters
    ----------
    channel : `str`
        name of new channel
    find_trend_source : `bool`, optional, default: `True`
        query for raw version of trend channel (trends not in CIS)
    timeout : `float`, optional, default: `5`
        number of seconds to wait before connection times out

    Returns
    -------
    Channel : :class:`~gwpy.detector.channel.Channel`
        new channel.
    """
    chans = re_channel.findall(str(channel))
    nchans = len(chans)
    if nchans > 1 or (nchans == 1 and chans[0] != str(channel)):
        name = str(channel)
        try:
            type_ = Channel.MATCH.match(name).groupdict()['type']
        except AttributeError:
            type_ = None
        # handle special characters in channel name
        rename = name
        for cchar in ['+', '*', '^', '|']:
            rename = rename.replace(cchar, '\%s' % cchar)
        found = globalv.CHANNELS.sieve(name=rename, exact_match=True)
    elif ',' in str(channel):
        name, type_ = str(channel).rsplit(',', 1)
        found = globalv.CHANNELS.sieve(name=name, type=type_, exact_match=True)
    else:
        type_ = isinstance(channel, Channel) and channel.type or None
        name = str(channel)
        found = globalv.CHANNELS.sieve(name=name, type=type_, exact_match=True)
    if len(found) == 1:
        return found[0]
    elif len(found) > 1:
        cstrings = set(['%s [%s, %s]' % (c.ndsname, c.sample_rate, c.unit)
                        for c in found])
        if len(cstrings) == 1:
            return found[0]
        else:
            raise ValueError("Ambiguous channel request '%s', multiple "
                             "existing channels recovered:\n    %s"
                         % (str(channel), '\n    '.join(cstrings)))
    else:
        matches = list(Channel.MATCH.finditer(name))
        # match single raw channel
        if len(matches) == 1 and not re.search('\.[a-z]+\Z', name):
            new = Channel(str(channel))
            new.url = '%s/channel/byname/%s' % (CIS_URL, str(new))
        # match single trend
        elif len(matches) == 1:
            # set default trend type based on mode
            if type_ is None and ':DMT-' in name:  # DMT is always m-trend
                type_ = 'm-trend'
            elif type_ is None and globalv.MODE == Mode.gps:
                type_ = 's-trend'
            elif type_ is None:
                type_ = 'm-trend'
            name += ',%s' % type_
            new = Channel(name)
            if find_trend_source:
                try:
                    source = get_channel(new.name.rsplit('.')[0])
                except ValueError:
                    pass
                else:
                    new.url = source.url
                    new.unit = source.unit
                    try:
                        new.bits = source.bits
                    except AttributeError:
                        pass
                    try:
                        new.filter = source.filter
                    except AttributeError:
                        pass
                    for param in filter(lambda x: x.endswith('_range') and
                                                  not hasattr(new, x),
                                        vars(source)):
                        setattr(new, param, getattr(source, param))
            # determine sample rate for trends
            if type_ == 'm-trend':
                new.sample_rate = 1/60.
            elif type_ == 's-trend':
                new.sample_rate = 1
        # match composite channel
        else:
            parts = get_channels([m.group() for m in matches])
            new = Channel(name)
            new.subchannels = parts
            new._ifo = "".join(set(p.ifo for p in parts if p.ifo))
        globalv.CHANNELS.append(new)
        try:
            return get_channel(new)
        except RuntimeError as e:
            if 'maximum recursion depth' in str(e):
                raise RuntimeError("Recursion error while accessing channel "
                                   "information for %s" % str(channel))
            else:
                raise


def get_channels(channels, **kwargs):
    """Multi-threaded channel query
    """
    if len(channels) == 0:
        return []

    # set up Queues
    inqueue = Queue()
    outqueue = Queue()

    # open threads
    for i in range(len(channels)):
        t = ThreadChannelQuery(inqueue, outqueue, **kwargs)
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


def update_missing_channel_params(channel, **kwargs):
    """Update empty channel parameters using the given input

    This method will only set parameters in the channel if the target
    parameter is `None`.

    Parameters
    ----------
    channel : `~gwpy.detector.Channel`
        channel to update
    **kwargs
        `(key, value)` pairs to set
    """
    target = get_channel(channel)
    for param in ['unit', 'sample_rate', 'frametype']:
        if getattr(target, param) is None:
            setattr(target, param, getattr(channel, param))
    for param in kwargs:
        if getattr(target, param) is None:
            setattr(target, param, kwargs[param])
    return target


def split(channelstring):
    """Split a comma-separated list of channels that may, or may not
    contain NDS2 channel types as well
    """
    out = []
    channelstring = re_quote.sub('', channelstring)
    while True:
        channelstring = channelstring.strip('\' \n')
        if ',' not in channelstring:
            break
        # check for complete line without NDS type
        line = channelstring.split('\n')[0].rstrip('\', \n')
        if not ',' in line:
            try:
                channelstring = channelstring.split('\n', 1)[1]
            except IndexError:
                pass
            else:
                out.append(line)
                continue
        # check for channel name with optional nds type
        for nds2type in NDS2_TYPES + ['']:
            if nds2type and ',%s' % nds2type in channelstring:
                try:
                    channel, ctype, channelstring = channelstring.split(',', 2)
                except ValueError:
                    channel, ctype = channelstring.split(',')
                    channelstring = ''
                out.append('%s,%s' % (channel, ctype))
                break
            elif nds2type == '' and ',' in channelstring:
                channel, channelstring = channelstring.split(',', 1)
                out.append(channel)
                break
    if channelstring:
        out.append(channelstring)
    return out


def split_combination(channelstring):
    """Split a math-combination of channels
    """
    channel = Channel(channelstring)
    if channel.ifo == 'G1':
        return channel.ndsname.split(' ')
    else:
        return re_channel.findall(channel.ndsname)
