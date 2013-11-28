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

from . import globalv

re_cchar = re.compile("[\W_]+")


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


def mkdir(path):
    """Conditional mkdir operation, for convenience
    """
    path = os.path.normpath(path)
    if not os.path.isdir(path):
        os.makedirs(path)
