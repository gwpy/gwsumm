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

"""Utilities for GWSumm
"""

import os
import sys
import time
import re
from multiprocessing import (cpu_count, active_children)

from gwpy.io import nds as ndsio

from . import globalv

re_cchar = re.compile("[\W_]+")
re_channel = re.compile(r'(?:[A-Z]\d:[A-Z1-9]+-\w+)(?:\.[a-zA-Z0-9]+)?')
re_quote = re.compile(r'^[\s\"\']+|[\s\"\']+$')
re_flagdiv = re.compile("(&|!=|!|\|)")

# define some colours
WARNC = '\033[93m'
ERRC = '\033[91m'
ENDC = '\033[0m'


def elapsed_time():
    """Return the time (seconds) since this job started
    """
    time.time() - globalv.START


def vprint(message, verbose=True, stream=sys.stdout, profile=True):
    """Prints the given message to the stream.

    Parameters
    ----------
    message : `str`
        string to print
    verbose : `bool`, optional, default: `True`
        flag to print or not, default: print
    stream : `file`, optional, default: `stdout`
        file object stream in which to print, default: stdout
    profile : `bool`, optional, default: `True`
        flag to print timestamp for debugging and profiling purposes
    """
    if stream != sys.stderr:
        profile &= globalv.PROFILE
        verbose &= globalv.VERBOSE
    if profile and message.endswith("\n"):
        message = "%s (%.2f)\n" % (message.rstrip("\n"), elapsed_time())
    if verbose:
        stream.write(message)
        stream.flush()


def mkdir(*paths):
    """Conditional mkdir operation, for convenience
    """
    for path in paths:
        path = os.path.normpath(path)
        if not os.path.isdir(path):
            os.makedirs(path)


def nat_sorted(l, key=None):
    """Sorted a list in the way that humans expect.

    Parameters
    ----------
    l : `iterable`
        iterable to sort
    key : `callable`
        sorting key

    Returns
    -------
    lsorted : `list`
        sorted() version of input ``l``
    """
    k = key and map(key, l) or l
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [convert(c) for c in
                                re.split('([0-9]+)', k[l.index(key)])]
    return sorted(l, key=alphanum_key)


def which(program):
    """Find full path of executable program
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file


def split_channels(channelstring):
    """Split a comma-separated list of channels that may, or may not
    contain NDS2 channel types as well
    """
    out = []
    channelstring = re_quote.sub('', channelstring)
    while True:
        channelstring = channelstring.strip('\' \n')
        if ',' not in channelstring:
            break
        for nds2type in ndsio.NDS2_CHANNEL_TYPE.keys() + ['']:
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


def count_free_cores(max=cpu_count()):
    """Count the number of CPU cores not currently used by this job.
    """
    active = 1 #len(active_children()
    return max - (active + 1)
