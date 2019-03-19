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
import re
from socket import getfqdn

from six import string_types

# import filter evals
from math import pi  # noqa: F401
import numpy  # noqa: F401

from . import globalv

re_cchar = re.compile(r"[\W_]+")
re_quote = re.compile(r'^[\s\"\']+|[\s\"\']+$')
re_flagdiv = re.compile(r"(&|!=|!|\|)")

# define some colours
WARNC = r'\033[93m'
ERRC = r'\033[91m'
ENDC = r'\033[0m'

# bad things to eval
UNSAFE_EVAL_STRS = [r'os\.(?![$\'\" ])', 'shutil', r'\.rm', r'\.mv']
UNSAFE_EVAL = re.compile(r'(%s)' % '|'.join(UNSAFE_EVAL_STRS))


# -- utilities ----------------------------------------------------------------

def elapsed_time():
    """Return the time (seconds) since this job started
    """
    import time
    return time.time() - globalv.START


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
    k = key and list(map(key, l)) or l

    def convert(text):
        if text.isdigit():
            return int(text)
        else:
            return text

    def alphanum_key(key):
        return [convert(c) for c in re.split('([0-9]+)', k[l.index(key)])]

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


_re_odc = re.compile('(OUTMON|OUT_DQ|LATCH)')


def get_odc_bitmask(odcchannel):
    return _re_odc.sub('BITMASK', str(odcchannel))


def safe_eval(val, strict=False, globals_=None, locals_=None):
    """Evaluate the given string as a line of python, if possible

    If the :meth:`eval` fails, a `str` is returned instead, unless
    `strict=True` is given.

    Parameters
    ----------
    val : `str`
        input text to evaluate

    strict : `bool`, optional, default: `False`
        raise an exception when the `eval` call fails (`True`) otherwise
        return the input as a `str` (`False`, default)

    globals_ : `dict`, optional
        dict of global variables to pass to `eval`, defaults to current
        `globals`

    locals_ : `dict`, optional
        dict of local variables to pass to `eval`, defaults to current
        `locals`

        .. note::

           Note the trailing underscore on the `globals_` and `locals_`
           kwargs, this is required to not clash with the builtin `globals`
           and `locals` methods`.

    Raises
    ------
    ValueError
        if the input string is considered unsafe to evaluate, normally
        meaning it contains something that might interact with the filesystem
        (e.g. `os.path` calls)
    NameError
    SyntaxError
        if the input cannot be evaluated, and `strict=True` is given

    See also
    --------
    eval
        for more documentation on the underlying evaluation method
    """
    # don't evaluate non-strings
    if not isinstance(val, string_types):
        return val
    # check that we aren't evaluating something dangerous
    try:
        match = UNSAFE_EVAL.search(val).group()
    except AttributeError:
        pass
    else:
        raise ValueError("Will not evaluate string containing %r: %r"
                         % (match, val))
    # format args for eval
    if globals_ is None and locals_ is None:
        args = ()
    elif globals_ is None and locals_ is not None:
        args = (globals(), locals_)
    elif locals_ is None and globals_ is not None:
        args = (globals_,)
    else:
        args = (globals_, locals_)
    # try and eval str
    try:
        return eval(val, *args)
    except (NameError, SyntaxError):
        return str(val)


# -- IFO parsing --------------------------------------------------------------

OBSERVATORY_MAP = {
    'G1': 'GEO',
    'H1': 'LIGO Hanford',
    'K1': 'KAGRA',
    'L1': 'LIGO Livingston',
    'V1': 'Virgo',
}


def get_default_ifo(fqdn=getfqdn()):
    """Find the default interferometer prefix (IFO) for the given host

    Parameters
    ----------
    fqdn : `str`
        the fully-qualified domain name (FQDN) of the host on which you
        wish to find the default IFO

    Returns
    -------
    IFO : `str`
        the upper-case X1-style prefix for the default IFO, if found, e.g. `L1`

    Raises
    ------
    ValueError
        if not default interferometer prefix can be parsed
    """
    if '.uni-hannover.' in fqdn or '.atlas.' in fqdn:
        return 'G1'
    elif '.ligo-wa.' in fqdn:
        return 'H1'
    elif '.ligo-la.' in fqdn:
        return 'L1'
    elif '.virgo.' in fqdn or '.ego-gw.' in fqdn:
        return 'V1'
    raise ValueError("Cannot determine default IFO for host %r" % fqdn)
