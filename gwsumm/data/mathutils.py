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

"""Handle arbitrary mathematical operations applied to data series
"""

import operator
import re

from gwpy.segments import SegmentList

from ..channels import (get_channel, re_channel)

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


# -- parse channel names that include mathematical operations -----------------

OPERATOR = {
    '*': operator.mul,
    '-': operator.sub,
    '+': operator.add,
    '/': operator.div,
    '^': operator.pow,
    '**': operator.pow,
}

re_math = re.compile('(?P<operator>.+?)'
                     '(?P<value>[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?)$')


def parse_math_definition(definition):
    """Parse the definition for a channel combination

    This method can only handle commutative operations, no fancy stuff
    with parentheses. Something like ``A*B`` is fine, but not ``(A+B)^2``

    Returns
    -------
    channels : `list` of `tuple`
        a list of 2-tuples containing the name of each channel, and any
        mathematical operations to be applied to that channel only
    operators : `list` of `callable`
        the list of functions that combine one channel and the previous,
        if `channels` is a list of length ``N``, then the `operators` list
        will have length ``N-1``

    Examples
    --------
    >>> parse_math_definition('H1:TEST * L1:TEST^2')
    ([('H1:TEST', None), ('L1:TEST', (<built-in function pow>, 2.0))], [<built-in function mul>])
    """
    breaks = re_channel.finditer(definition)
    channels = []
    operators = []
    try:
        match = next(breaks)
    except StopIteration:  # no channel names parsed at all, just return
        return [(definition, None)], []
    while True:
        # find channel
        a, b = match.span()
        cname = definition[a:b]
        channels.append((cname, None))

        # find next channel and parse whatever's inbetween
        try:
            match = next(breaks)
        except StopIteration:
            c = None
        else:
            c = match.span()[0]

        # parse math string
        mathstr = definition[b:c].strip()
        if not mathstr and c is None:  # if at end, break
            break
        elif not mathstr:  # otherwise no operator, panic
            raise ValueError("Cannot parse math operator between channels "
                             "in definition %r" % definition)
        if c is not None:  # if between channels, find the combiner
            operators.append(get_operator(mathstr[-1]))
            mathstr = mathstr[:-1].strip()
        try:  # then parse some math to apply to the current channel
            cop = re_math.match(mathstr).groupdict()
        except AttributeError:
            if mathstr:
                raise ValueError("Cannot parse math operation %r" % mathstr)
        else:
            op = get_operator(cop['operator'].strip())
            value = float(cop['value'])
            channels[-1] = (cname, (op, value))  # record math to be done
        if c is None:  # if at the end, break
            break
    return channels, operators


def get_operator(opstr):
    try:
        return OPERATOR[opstr]
    except KeyError:
        raise ValueError("Cannot parse math operator %r" % opstr)


def get_with_math(channel, segments, load_func, get_func, **ioargs):
    """Get data with optional arbitrary math definitions

    Parameters
    ----------
    channel : `str`
        name of the meta-channel to create
    segments : `~gwpy.segments.SegmentList`
        segments over which to create the new channel
    load_func : `callable`
        method to call to load data from disk
    get_func : `callable`
        method to call to return channel data
    **ioargs
        all other kwargs are passed to the `load_func` and `get_func`

    Returns
    -------
    datalist : `TimeSeriesList`, or similar
        a structured list of data objects, probably either for `TimeSeries`
        or `Spectrogram`
    """
    # parse definition
    singleops, joinoperators = parse_math_definition(str(channel))
    channel = get_channel(channel)
    names = zip(*singleops)[0]
    chans = map(get_channel, names)
    # get raw data
    if load_func is get_func:  # if load_func returns a single channel
        tsdict = dict((c.ndsname, load_func(c, segments, **ioargs))
                      for c in chans)
    else:
        tsdict = load_func(chans, segments, **ioargs)
    # shortcut single channel with no math
    if len(names) == 1 and singleops[0][1] is None:
        return tsdict.values()[0]
    # get union of segments for all sub-channels
    tslist = [tsdict[c.ndsname] for c in chans]
    datasegs = reduce(operator.and_, [tsl.segments for tsl in tslist])
    # build meta-timeseries for all intersected segments
    meta = type(tsdict.values()[0])()
    for seg in datasegs:
        # get data for first channel
        ts, = get_func(names[0], SegmentList([seg]), **ioargs)
        ts.name = str(channel)
        # apply math to this channel
        cmath = singleops[0][1]
        if cmath is not None:
            ts = cmath[0](ts, cmath[1])
        # for each other channel do the same
        for joinop, ch in zip(joinoperators, singleops[1:]):
            name, cmath = ch
            data, = get_func(name, SegmentList([seg]), **ioargs)
            if cmath is not None:  # apply simple math
                data = cmath[0](data, cmath[1])
            ts = joinop(ts, data)  # apply combination math
        meta.append(ts)
    return meta
