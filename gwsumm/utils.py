
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