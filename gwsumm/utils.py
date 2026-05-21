# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013)
#               Evan Goetz (2026)
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

import argparse
import datetime
from dateutil.relativedelta import relativedelta
import os
import re
import sys
from socket import getfqdn

# import filter evals
from math import pi  # noqa: F401
import numpy  # noqa: F401

from gwpy.time import to_gps

from . import globalv, __version__

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

PROG = ('python -m gwsumm' if sys.argv[0].endswith('.py')
        else os.path.basename(sys.argv[0]))


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


def nat_sorted(iterable, key=None):
    """Sorted a list in the way that humans expect.

    Parameters
    ----------
    iterable : `iterable`
        iterable to sort
    key : `callable`
        sorting key

    Returns
    -------
    lsorted : `list`
        sorted() version of input ``l``
    """
    k = key and list(map(key, iterable)) or iterable

    def convert(text):
        if text.isdigit():
            return int(text)
        else:
            return text

    def alphanum_key(key):
        return [convert(c) for c in re.split(
            '([0-9]+)', k[iterable.index(key)])]

    return sorted(iterable, key=alphanum_key)


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
    if not isinstance(val, str):
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


# -- argparse utilities -------------------------------------------------------

class GWArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(GWArgumentParser, self).__init__(*args, **kwargs)
        self._positionals.title = 'Positional arguments'
        self._optionals.title = 'Optional arguments'


class GWHelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, *args, **kwargs):
        kwargs.setdefault('indent_increment', 4)
        super(GWHelpFormatter, self).__init__(*args, **kwargs)


class DateAction(argparse.Action):
    TIMESCALE = {'days': 1}

    @staticmethod
    def set_gps_times(namespace, startdate, enddate):
        setattr(namespace, 'gpsstart', to_gps(startdate))
        setattr(namespace, 'gpsend', to_gps(enddate))

    def __call__(self, parser, namespace, values, option_string=None):
        try:
            date = datetime.datetime.strptime(values, self.DATEFORMAT)
        except ValueError:
            raise parser.error("%s malformed: %r. Please format as %s"
                               % (self.dest.title(), values, self.METAVAR))
        else:
            self.set_gps_times(namespace, date,
                               date + relativedelta(**self.TIMESCALE))
            setattr(namespace, self.dest, date)
        return date


class DayAction(DateAction):
    TIMESCALE = {'days': 1}
    DATEFORMAT = '%Y%m%d'
    METAVAR = 'YYYYMMDD'


class WeekAction(DayAction):
    TIMESCALE = {'days': 7}
    DATEFORMAT = '%Y%m%d'
    METAVAR = 'YYYYMMDD'


class MonthAction(DateAction):
    TIMESCALE = {'months': 1}
    DATEFORMAT = '%Y%m'
    METAVAR = 'YYYYMM'


class YearAction(DateAction):
    TIMESCALE = {'years': 1}
    DATEFORMAT = '%Y'
    METAVAR = 'YYYY'


class GPSAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=False):
        try:
            values = float(values)
        except (TypeError, ValueError):
            pass
        setattr(namespace, self.dest, to_gps(values))


# -- parse command-line -------------------------------------------------------

def add_configuration_options(sharedopts_):
    try:
        default_ifo = get_default_ifo()
    except ValueError:
        default_ifo = None

    # configuration arguments
    copts = sharedopts_.add_argument_group(
        "Configuration options",
        "Provide a number of INI-format configuration files",
    )
    copts.add_argument(
        '-i',
        '--ifo',
        type=str,
        default=default_ifo,
        metavar='IFO',
        help="IFO prefix for interferometer to process. "
             "If this option is set in the [DEFAULT] of any of "
             "the INI files, giving it here is redundant.",
    )
    copts.add_argument(
        '-g',
        '--global-config',
        action='append',
        type=str,
        metavar='FILE',
        default=[],
        help="INI file for use in all workflow jobs, may be given "
             "multiple times",
    )
    copts.add_argument(
        '-f',
        '--config-file',
        action='append',
        type=str,
        metavar='FILE',
        default=[],
        help="INI file for analysis, may be given multiple times",
    )
    copts.add_argument(
        '-t',
        '--process-tab',
        action='append',
        type=str,
        help="process only this tab, can be given multiple times",
    )


def add_process_options(sharedopts_):
    popts = sharedopts_.add_argument_group(
        "Process options",
        "Configure how this summary will be processed.",
    )
    popts.add_argument(
        '--nds',
        action='store_true',
        help="use NDS as the data source",
    )
    popts.add_argument(
        '-j',
        '--multi-process',
        type=int,
        default=1,
        dest='multiprocess',
        metavar='N',
        help="use a maximum of N parallel processes at any time",
    )
    popts.add_argument(
        '-b',
        '--bulk-read',
        action='store_true',
        help="read all data up-front at the start of the job, "
             "rather than when it is needed for a tab",
    )
    popts.add_argument(
        '-S',
        '--on-segdb-error',
        type=str,
        default='raise',
        choices=['raise', 'ignore', 'warn'],
        help="action upon error fetching segments from SegDB",
    )
    popts.add_argument(
        '-G',
        '--on-datafind-error',
        type=str,
        default='raise',
        choices=['raise', 'ignore', 'warn'],
        help="action upon error querying for frames from the "
             "datafind server, default: %(default)s",
    )
    popts.add_argument(
        '--data-cache',
        action='append',
        default=[],
        help='path to LAL-format cache of TimeSeries data files',
    )
    popts.add_argument(
        '--event-cache',
        action='append',
        default=[],
        help='path to LAL-format cache of event trigger files',
    )
    popts.add_argument(
        '--segment-cache',
        action='append',
        default=[],
        help='path to LAL-format cache of state '
             'or data-quality segment files',
    )


def add_output_options(parser_):
    """Add output options to the subparser

    This is only needed because argparse can't handle mutually exclusive
    groups in a parent parser handed to a subparser for some reason.
    """
    outopts = parser_.add_argument_group("Output options")
    outopts.add_argument(
        '-o',
        '--output-dir',
        type=str,
        metavar='DIR',
        default=os.curdir,
        help="Output directory for summary information",
    )
    htmlopts = outopts.add_mutually_exclusive_group()
    htmlopts.add_argument(
        '-M',
        '--html-only',
        action='store_true',
        help="Generate container HTML and navigation only",
    )
    htmlopts.add_argument(
        '-n',
        '--no-html',
        action='store_true',
        help="Generate inner HTML and contents only, not supporting HTML",
    )
    outopts.add_argument(
        '-N',
        '--no-htaccess',
        action='store_true',
        help="don't create a .htaccess file to customise 404 errors",
    )


def add_archive_options(parser_, archive_dir_str):
    """Add archiving options to the subparser

    This is only needed because argparse can't handle mutually exclusive
    groups in a parent parser handed to a subparser for some reason.
    """
    hierarchopts = parser_.add_argument_group('Archive options')
    hierarchchoice = hierarchopts.add_mutually_exclusive_group()
    hierarchchoice.add_argument(
        '-a',
        '--archive',
        metavar='FILE_TAG',
        default=False,
        const='GW_SUMMARY_ARCHIVE',
        nargs='?',
        help="Read archived data from, and write processed data to "
             "an HDF archive file written with the FILE_TAG. If not "
             "given, no archive will be used, if given with no file "
             "tag, a default of '%(const)s' will be used.",
    )
    hierarchchoice.add_argument(
        '-d',
        '--daily-archive',
        metavar='FILE_TAG',
        default=False,
        const='GW_SUMMARY_ARCHIVE',
        nargs='?',
        help="Read data from the daily archives, with the given FILE_TAG."
             "If not given, daily archives will be used, if given with no "
             "file tag, a default of '%(const)s' will be used.",
    )
    hierarchopts.add_argument(
        '--archive-read-dir',
        metavar='DIR',
        type=str,
        default=archive_dir_str,
        help="Read archived data from this directory",
    )
    hierarchopts.add_argument(
        '--archive-write-dir',
        metavar='DIR',
        type=str,
        default=archive_dir_str,
        help="Write archived data to this directory",
    )


def add_condor_options(sharedopts_):
    """Add condor options to the subparser"""
    htcopts = sharedopts_.add_argument_group('Condor options')
    htcopts.add_argument(
        '--container-path',
        type=str,
        required=True,
        help="Path to the container image file (.sif)"
    )
    htcopts.add_argument(
        '-T',
        '--file-tag',
        type=str,
        default='gw_summary_pipe',
        help="file tag for pipeline files, default: %(default)s",
    )
    htcopts.add_argument(
        '-u',
        '--universe',
        type=str,
        default='container',
        choices=['container', 'local'],
        help="Universe for condor jobs, default: %(default)s",
    )
    htcopts.add_argument(
        '-l',
        '--log-dir',
        type=str,
        default=os.environ.get('LOCALDIR', None),
        help="Directory path for condor log files, default: %(default)s",
    )
    htcopts.add_argument(
        '-m',
        '--maxjobs',
        type=int,
        metavar='N',
        help="Restrict the DAG to submit only N jobs at any one "
             "time, default: %(default)s",
    )
    htcopts.add_argument(
        '-O',
        '--condor-timeout',
        type=float,
        metavar='T',
        help='Configure condor to terminate jobs after T hours '
             'to prevent idling, default: %(default)s',
    )
    htcopts.add_argument(
        '-c',
        '--condor-command',
        action='append',
        type=str,
        default=[],
        help="Extra condor submit commands to add to gw_summary submit file. "
             "Can be given multiple times in the form \"key=value\"",
    )
    htcopts.add_argument(
        '--condor-command-container',
        action='append',
        type=str,
        default=[],
        help="Extra condor submit commands to add to container universe "
             "gw_summary submit file. Can be given multiple times in the "
             "form \"key=value\"",
    )
    htcopts.add_argument(
        '-p',
        '--priority',
        action='append',
        type=str,
        default=[],
        help="priority for DAG node, should be given once "
             "for each --config-file in the same order",
    )


def create_parser(include_condor_opts=False):
    """Create a command-line parser for this entry point
    """
    # initialize top-level argument parser
    parser = GWArgumentParser(
        formatter_class=GWHelpFormatter,
        prog=PROG,
        description=__doc__,
        fromfile_prefix_chars="@",
        epilog="Arguments and options may be written into files and passed as "
               "positional arguments prepended with '@', e.g. '%(prog)s "
               "@args.txt'. In this format, options must be give as "
               "'--argument=value', and not '--argument value'.",
    )

    # global arguments
    parser.add_argument(
        '-V',
        '--version',
        action='version',
        version=__version__,
        help="show program's version number and exit",
    )

    # shared arguments
    sharedopts = GWArgumentParser(add_help=False)
    sharedopts.title = 'Progress arguments'
    sharedopts.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        default=False,
        help="show verbose logging output",
    )
    sharedopts.add_argument(
        '-D',
        '--debug',
        action='store_true',
        default=False,
        help="show information that could be useful in debugging",
    )

    # configuration options
    add_configuration_options(sharedopts)

    # process options
    add_process_options(sharedopts)

    # condor options, if needed
    if include_condor_opts:
        add_condor_options(sharedopts)

    # define sub-parser handler
    subparsers = parser.add_subparsers(
        dest='mode',
        title='Modes',
        description='Note: all dates are defined with UTC boundaries.\n'
                    'The valid modes are:',
    )
    subparser = dict()

    # DAY mode
    today = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%d')
    daydoc = """
    Run %s over a full UTC day, and link this day to others with a calendar
    built into the HTML navigation bar. In this mode you can also archive data
    in HDF-format files to allow progressive processing of live data without
    restarting from scratch every time.""" % parser.prog
    subparser['day'] = subparsers.add_parser(
        'day',
        description=daydoc,
        epilog=parser.epilog,
        parents=[sharedopts],
        formatter_class=GWHelpFormatter,
        help="Process one day of data",
    )
    subparser['day'].add_argument(
        'day',
        action=DayAction,
        type=str,
        nargs='?',
        metavar=DayAction.METAVAR,
        default=today,
        help="Day to process",
    )
    add_output_options(subparser['day'])

    darchopts = subparser['day'].add_argument_group(
        'Archive options',
        'Choose if, and how, to archive data from this run',
    )
    darchopts.add_argument(
        '-a',
        '--archive',
        metavar='FILE_TAG',
        default=False,
        const='GW_SUMMARY_ARCHIVE',
        nargs='?',
        help="Read archived data from, and write processed "
             "data to, an HDF archive file written with the "
             "FILE_TAG. If not given, no archive will be used, "
             "if given with no file tag, a default of "
             "'%(const)s' will be used.",
    )
    darchopts.add_argument(
        '--archive-read-dir',
        metavar='DIR',
        type=str,
        default=f'day/{today}/archive',
        help="Read archived data from this directory",
    )
    darchopts.add_argument(
        '--archive-write-dir',
        metavar='DIR',
        type=str,
        default=f'day/{today}/archive',
        help="Write archived data into this directory",
    )

    # WEEK mode
    subparser['week'] = subparsers.add_parser(
        'week',
        parents=[sharedopts],
        epilog=parser.epilog,
        formatter_class=GWHelpFormatter,
        help="Process one week of data",
    )
    subparser['week'].add_argument(
        'week',
        action=WeekAction,
        type=str,
        metavar=WeekAction.METAVAR,
        help="Week to process (given as starting day)",
    )
    add_output_options(subparser['week'])
    add_archive_options(subparser['week'], f'week/{today}/archive')

    # MONTH mode
    today = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m')
    subparser['month'] = subparsers.add_parser(
        'month',
        parents=[sharedopts],
        epilog=parser.epilog,
        formatter_class=GWHelpFormatter,
        help="Process one month of data",
    )
    subparser['month'].add_argument(
        'month',
        action=MonthAction,
        type=str,
        metavar=MonthAction.METAVAR,
        help="Month to process",
    )
    add_output_options(subparser['month'])
    add_archive_options(subparser['month'], f'month/{today}/archive')

    # GPS mode
    subparser['gps'] = subparsers.add_parser(
        'gps',
        parents=[sharedopts],
        epilog=parser.epilog,
        formatter_class=GWHelpFormatter,
        help="Process GPS interval",
    )
    subparser['gps'].add_argument(
        'gpsstart',
        action=GPSAction,
        type=str,
        metavar='GPSSTART',
        help='GPS start time',
    )
    subparser['gps'].add_argument(
        'gpsend',
        action=GPSAction,
        type=str,
        metavar='GPSEND',
        help='GPS end time.',
    )
    add_output_options(subparser['gps'])

    garchopts = subparser['gps'].add_argument_group(
        'Archive options',
        'Choose if, and how, to archive data from this run',
    )
    garchopts.add_argument(
        '-a',
        '--archive',
        metavar='FILE_TAG',
        default=False,
        const='GW_SUMMARY_ARCHIVE',
        nargs='?',
        help="Read archived data from, and write processed "
             "data to, an HDF archive file written with the "
             "FILE_TAG. If not given, no archive will be used, "
             "if given with no file tag, a default of "
             "'%(const)s' will be used.")
    garchopts.add_argument(
        '--archive-read-dir',
        metavar='DIR',
        type=str,
        default=None,
        help="Read archived data from this directory, default:"
             "gps/<start>-<stop>/archive",
    )
    garchopts.add_argument(
        '--archive-write-dir',
        metavar='DIR',
        type=str,
        default=None,
        help="Write archived data into this directory, default:"
             "gps/<start>-<stop>/archive",
    )

    # return the argument parser
    return parser
