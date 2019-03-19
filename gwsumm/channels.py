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

import re
from functools import wraps

from astropy.units import Unit

from gwpy.detector import (Channel, ChannelList)
from gwpy.io.nds2 import Nds2ChannelType

from . import globalv
from .utils import re_quote

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

NDS2_TYPES = Nds2ChannelType.names()
CIS_URL = 'https://cis.ligo.org'
re_channel = re.compile(r'[A-Z]\d:[a-zA-Z0-9]+'  # core channel section L1:TEST
                        r'(?:[-_][a-zA-Z0-9_]+)?'  # underscore-delimited parts
                        r'(?:\.[a-z]+)?'  # trend type
                        r'(?:,[a-z-]+)?')  # NDS channel type


# -- channel creation ---------------------------------------------------------

def _match(channel):
    channel = Channel(channel)
    name = str(channel)
    type_ = channel.type
    found = globalv.CHANNELS.sieve(
        name=name, type=type_, exact_match=True)

    # if match, return now
    if found:
        return found

    # if no matches, try again without matching type
    found = globalv.CHANNELS.sieve(name=name, exact_match=True)
    if len(found) == 1:
        # found single match that is less specific, so we make it more
        # specific. If someone else wants another type for the sme channel
        # this is fine, and it will create another channel.
        found[0].type = type_
    return found


def _find_parent(channel):
    """Find the parent for a given channel

    This is either the raw channel from which a trend was generated, or the
    real channel for a mathematically-manipulated channel.

    Raises
    ------
    ValueError
        if the parent cannot be parsed
    """
    if channel.trend:
        parent = str(channel).rsplit('.')[0]
    else:
        parent, = re_channel.findall(str(channel))
    if parent == str(channel):
        raise ValueError("Cannot find parent for '{0!s}'".format(channel))
    return get_channel(parent)


def _update_dependent(channel):
    """Update a trend channel from its parent
    """
    try:
        source = _find_parent(channel)
    except (ValueError, IndexError):
        return channel
    if source is channel:
        return channel

    channel.url = source.url
    channel.unit = source.unit

    # copy custom named params
    #     this works because the upstream Channel stores all
    #     attributes in private variable names
    # NOTE: we need to exclude 'resample' to not attempt to upsample trends
    for param in filter(
            lambda x: not x.startswith('_') and
            not hasattr(channel, x) and x not in ('resample',),
            vars(source)):
        try:
            setattr(channel, param, getattr(source, param))
        except AttributeError:
            pass
    return channel


def _with_update_dependent(func):
    """Decorate ``func`` to call `_update_dependent()` upon exit
    """
    @wraps(func)
    def wrapped_func(*args, **kwargs):
        _update = kwargs.pop('find_parent', True)
        out = func(*args, **kwargs)
        if _update and out.trend:
            out = _update_dependent(out)
        return out
    return wrapped_func


def _new(channel, find_parent=True):
    """Create a new `~gwpy.detector.Channel` in the globalv cache.
    """
    # convert to Channel
    if isinstance(channel, Channel):
        new = channel
    else:
        new = Channel(channel)
    name = str(channel)
    type_ = new.type

    # work out what kind of channel it is
    parts = re_channel.findall(name)

    # match single raw channel for LIGO
    if (
            len(parts) == 1 and
            new.ifo in ('H1', 'L1') and
            not re.search(r'\.[a-z]+\Z', name)
    ):
        new.url = '%s/channel/byname/%s' % (CIS_URL, str(new))

    # match single trend
    elif len(parts) == 1:
        # set default trend type based on mode
        if type_ is None and ':DMT-' in name:  # DMT is always m-trend
            new.type = 'm-trend'
        # match parameters from 'raw' version of this channel

    # match composite channel
    else:
        new.subchannels = parts
        new._ifo = "".join(set(p.ifo for p in map(Channel, parts) if p.ifo))

    if find_parent:
        _update_dependent(new)

    # store new channel and return
    globalv.CHANNELS.append(new)
    try:
        return get_channel(new)
    except RuntimeError as e:
        if 'maximum recursion depth' in str(e):
            raise RuntimeError("Recursion error while accessing channel "
                               "information for %s" % str(channel))
        raise


@_with_update_dependent
def get_channel(channel, find_parent=True, timeout=5):
    """Find (or create) a :class:`~gwpy.detector.Channel`.

    If ``channel`` has already been created, the cached copy will be
    returned, otherwise a new `~gwpy.detector.Channel` object will be created.

    Parameters
    ----------
    channel : `str`
        name of new channel

    find_parent : `bool`, optional, default: `True`
        query for raw version of trend channel (trends not in CIS)

    timeout : `float`, optional, default: `5`
        number of seconds to wait before connection times out

    Returns
    -------
    Channel : :class:`~gwpy.detector.Channel`
        new channel.
    """
    chans = re_channel.findall(str(channel))
    nchans = len(chans)

    # match compound channel name
    if nchans > 1 or (nchans == 1 and chans[0] != str(channel)):
        name = str(channel)
        # handle special characters in channel name
        rename = name
        for cchar in ['+', '*', '^', '|']:
            rename = rename.replace(cchar, r'\%s' % cchar)
        found = globalv.CHANNELS.sieve(name=rename, exact_match=True)
    # match normal channel
    else:
        found = _match(channel)

    # if single match, return it
    if len(found) == 1:
        return found[0]

    # if multiple matches raise error (unless weird circumstances)
    elif len(found) > 1:
        cstrings = set(['%s [%s, %s]' % (c.ndsname, c.sample_rate, c.unit)
                        for c in found])
        if len(cstrings) == 1:  # somehow all channels are the same
            return found[0]
        raise ValueError("Ambiguous channel request '%s', multiple "
                         "existing channels recovered:\n    %s"
                         % (str(channel), '\n    '.join(cstrings)))

    # otherwise there were no matches, so we create a new channel
    return _new(channel, find_parent=find_parent)


def get_channels(channels, **kwargs):
    """Find (or create) multiple channels.

    See Also
    --------
    get_channel
    """
    return ChannelList(get_channel(c, **kwargs) for c in channels)


# -- channel manipulation -----------------------------------------------------

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
    target = get_channel(str(channel))
    if isinstance(channel, Channel):
        for param in ['unit', 'sample_rate', 'frametype']:
            if getattr(target, param) is None or (
                    param == 'unit' and
                    getattr(target, param) is Unit('undef')):
                setattr(target, param, getattr(channel, param))
    for param in kwargs:
        if getattr(target, param) is None:
            setattr(target, param, kwargs[param])
    return target


def update_channel_params():
    """Update the `globalv.CHANNELS` list based on internal parameter changes

    This is required to update `Channel.type` based on `Channel.frametype`,
    and similar.
    """
    for c in globalv.CHANNELS:
        # update type based on frametype
        if c.type is None and c.frametype == '{0.ifo}_M'.format(c):
            c.type = 'm-trend'
        elif c.type is None and c.frametype == '{0.ifo}_T'.format(c):
            c.type = 's-trend'

        # update sample_rate based on trend type
        if c.type == 'm-trend' and c.sample_rate is None:
            c.sample_rate = 1/60.
        elif c.type == 's-trend' and c.sample_rate is None:
            c.sample_rate = 1.
    return


# -- string parsing -----------------------------------------------------------

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
        if ',' not in line:
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

    Returns
    -------
    channels : `~gwpy.detector.ChannelList`
    """
    return get_channels(re_channel.findall(str(channelstring)),
                        find_parent=False)
