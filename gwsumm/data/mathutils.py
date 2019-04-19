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

import numbers
import operator
import re
from collections import OrderedDict

from six.moves import reduce

import numpy

from gwpy.segments import SegmentList
from gwpy.frequencyseries import FrequencySeries
from gwpy.types import Series

from ..channels import get_channel

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'


# -- parse channel names that include mathematical operations -----------------

OPERATOR = {
    '*': operator.mul,
    '-': operator.sub,
    '+': operator.add,
    '/': operator.truediv,
    '^': operator.pow,
    '**': operator.pow,
}

re_operator = re.compile(r'\s+[+/^\*-]+\s+')
re_value = re.compile(r'[-+]?(\d+(\.\d*)?|\.\d+)([eE][-+]?\d+)?')


def parse_math_definition(definition):
    """Parse the definition for a channel combination

    This method can only handle commutative operations, no fancy stuff
    with parentheses. Something like ``A * B`` is fine, but not ``(A + B) ^ 2``

    All operands, operators, and values should be space-separated.

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
    ([('H1:TEST', None), ('L1:TEST', (<built-in function pow>, 2.0))],
     [<built-in function mul>])
    """
    channels = OrderedDict()
    operators = []
    ops = re_operator.finditer(definition)
    try:
        match = next(ops)
    except StopIteration:
        return OrderedDict([(definition, None)]), []
    x = 0
    while True:
        # parse operator position and method
        a, b = match.span()
        op = get_operator(match.group().strip())

        # parse before operator
        before = definition[x:a]

        # record first channel
        if not x:
            key = before
            channels[key] = []

        # find next operator
        try:
            match = next(ops)
        except StopIteration:  # no further operators
            c = None
        else:  # operator found
            c, d = match.span()

        # parse after operator
        after = definition[b:c]
        x = b  # truncate definition for next time

        # match value after operator
        try:
            a2, b2 = re_value.match(after).span()
        except AttributeError:  # not a value (must be channel name)
            # record joining operator
            operators.append(op)
            # record new channel
            key = definition[b:c]
            channels[key] = []
        else:
            # parse value as a float and pair with operator
            value = float(after[a2:b2])
            math = (op, value)
            channels[key].append(math)

        # no more operators, so we're done
        if c is None:
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
    names = list(singleops)
    chans = list(map(get_channel, names))
    # get raw data
    if load_func is get_func:  # if load_func returns a single channel
        tsdict = dict((c.ndsname, load_func(c, segments, **ioargs))
                      for c in chans)
    else:
        tsdict = load_func(chans, segments, **ioargs)
    # shortcut single channel with no math
    if len(names) == 1 and not singleops[names[0]]:
        vals = list(tsdict.values())
        if isinstance(vals[0], list):
            return vals[0]
        else:
            return vals
    # get union of segments for all sub-channels
    tslist = [tsdict[c.ndsname] for c in chans]
    if isinstance(tslist[0], FrequencySeries):
        datasegs = segments
        meta = []
    else:
        datasegs = reduce(operator.and_, [tsl.segments for tsl in tslist])
        meta = type(list(tsdict.values())[0])()
    for seg in datasegs:
        # get data for first channel
        if isinstance(tslist[0], FrequencySeries):
            ts = get_func(names[0], SegmentList([seg]), **ioargs)
        else:
            ts, = get_func(names[0], SegmentList([seg]), **ioargs)
        ts.name = str(channel)
        # apply math to this channel
        cmath = singleops[names[0]]
        for op_, val_ in cmath:
            if op_ in {operator.add, operator.sub} and \
                    isinstance(val_, numbers.Number):
                val_ *= ts.unit
            ts = op_(ts, val_)
        # for each other channel do the same
        for joinop, name in zip(joinoperators, names[1:]):
            cmath = singleops[name]
            if isinstance(tslist[0], FrequencySeries):
                data = get_func(name, SegmentList([seg]), **ioargs)
            else:
                data, = get_func(name, SegmentList([seg]), **ioargs)
            for op_, val_ in cmath:
                data = op_(data, val_)
            ts = _join(ts, data, joinop)
        meta.append(ts)
    if not meta and isinstance(tslist[0], Series):
        meta.append(type(tslist[0])([], channel=channel))
    if not meta:
        return type(tslist[0])()
    return meta


def _join(a, b, op):
    """Method to combine two data structures, handling shape mismatches

    This method is for internal use only, and should not be called from
    outside
    """
    # crop time-axis to select overlapping data
    if a.xunit.physical_type == 'time' and (
            abs(a.xspan) != abs(b.xspan) or a.shape[0] != b.shape[0]
    ):
        overlap = a.xspan & b.xspan
        a = a.crop(*overlap)
        b = b.crop(*overlap)
    # try and join now
    try:
        return op(a, b)
    # handle mismatched frequency scale
    except ValueError as e:
        # if error is _not_ a shape mismatch, raise
        if 'operands could not be broadcast' not in str(e):
            raise
        # if the FFTlength is not the same, raise (no interpolation)
        if a.df != b.df or a.f0 != b.f0:
            raise
        # otherwise, lengthen the shorted array in frequency
        if a.ndim == 1:  # frequencyseries
            nf = a.size - b.size
            if nf > 0:
                c = numpy.require(b, requirements=['O'])
                b = numpy.resize(b, (c.size - nf,))
            elif nf < 0:
                c = numpy.require(a, requirements=['O'])
                a = numpy.resize(a, (c.size - nf,))
        else:  # spectrogram
            nf = a.shape[1] - b.shape[1]
            new = numpy.zeros((a.shape[0], abs(nf)))
            if nf > 0:  # reshape 'b'
                del b.frequencies
                b2 = numpy.concatenate((b, new), axis=1)
                b2.__array_finalize__(b)
                b = b2
            elif nf < 0:  # reshape 'a'
                del a.frequencies
                a2 = numpy.concatenate((a, new), axis=1)
                a2.__array_finalize__(a)
                a = a2
        return op(a, b)
