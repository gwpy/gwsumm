# coding=utf-8
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

"""Command-line interface for the gravitational-wave interferometer summary
information system (`gwsumm`).

This module provides the command-line interface to the GWSumm package,
allowing generation of detector summary information.

Select a <mode> to run over a calendar amount of time ('day', 'week',
or 'month'), or an arbitrary GPS (semi-open) interval.

Run '%(prog)s <mode> --help' for details of the specific arguments and
options available for each mode.
"""

import argparse
import calendar
import datetime
import getpass
import os
import re
import sys
import warnings

from collections import OrderedDict
from configparser import (DEFAULTSECT, NoOptionError, NoSectionError)
from dateutil.relativedelta import relativedelta
from urllib.parse import urlparse

from glue.lal import (Cache, LIGOTimeGPS)

from gwpy import time
from gwpy.segments import (Segment, SegmentList)
from gwpy.signal.spectral import _lal as fft_lal
from gwpy.time import (Time, _tconvert, tconvert, to_gps)

from gwdetchar.cli import logger

from . import (
    __version__,
    archive,
    globalv,
    mode,
)
from .config import (
    GWSummConfigParser,
)
from .segments import get_segments
from .state import (
    ALLSTATE
)
from .tabs import (
    TabList,
    get_tab,
)
from .utils import (
    get_default_ifo,
    mkdir,
    re_flagdiv,
)
from .data import get_timeseries_dict

# set matplotlib backend
from matplotlib import use
use('Agg')

# XXX HACK XXX
# plots don't work with multiprocessing with lal.LIGOTimeGPS
time.LIGOTimeGPS = LIGOTimeGPS
_tconvert.LIGOTimeGPS = LIGOTimeGPS

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'

# set defaults
VERBOSE = False
PROFILE = False
try:
    DEFAULT_IFO = get_default_ifo()
except ValueError:
    DEFAULT_IFO = None

# initialize logger
PROG = ('python -m gwsumm' if sys.argv[0].endswith('.py')
        else os.path.basename(sys.argv[0]))
LOGGER = logger(name=PROG.split('python -m ').pop())

# find today's date
TODAY = datetime.datetime.utcnow().strftime('%Y%m%d')


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

def add_output_options(parser_):
    """Add outuput options to the subparser

    This is only needed because argparse can't handle mutually exclusive
    groups in a parent parser handed to a subparser for some reason.
    """
    outopts = parser_.add_argument_group("Output options")
    outopts.add_argument(
        '-o',
        '--output-dir',
        action='store',
        type=str,
        metavar='DIR',
        default=os.curdir,
        help="Output directory for summary information",
    )
    htmlopts = outopts.add_mutually_exclusive_group()
    htmlopts.add_argument(
        '-m',
        '--html-only',
        action='store_true',
        default=False,
        help="Generate container HTML and navigation only",
    )
    htmlopts.add_argument(
        '-n',
        '--no-html',
        action='store_true',
        default=False,
        help="Generate inner HTML and contents only, not supporting HTML",
    )
    outopts.add_argument(
        '-N',
        '--no-htaccess',
        action='store_true',
        default=False,
        help="don't create a .htaccess file to customise 404 errors",
    )


def add_archive_options(parser_):
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


def create_parser():
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

    # configuration arguments
    copts = sharedopts.add_argument_group(
        "Configuration options",
        "Provide a number of INI-format configuration files",
    )
    copts.add_argument(
        '-i',
        '--ifo',
        default=DEFAULT_IFO,
        metavar='IFO',
        help="IFO prefix for interferometer to process. "
             "If this option is set in the [DEFAULT] of any of "
             "the INI files, giving it here is redundant.",
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

    # process options
    popts = sharedopts.add_argument_group(
        "Process options",
        "Configure how this summary will be processed.",
    )
    popts.add_argument(
        '--nds',
        action='store_true',
        default=None,
        help="use NDS as the data source, default: 'guess'",
    )
    popts.add_argument(
        '-j',
        '--multi-process',
        action='store',
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
        default=False,
        help="read all data up-front at the start of the job, "
             "rather than when it is needed for a tab",
    )
    popts.add_argument(
        '-S',
        '--on-segdb-error',
        action='store',
        type=str,
        default='raise',
        choices=['raise', 'ignore', 'warn'],
        help="action upon error fetching segments from SegDB",
    )
    popts.add_argument(
        '-G',
        '--on-datafind-error',
        action='store',
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

    # define sub-parser handler
    subparsers = parser.add_subparsers(
        dest='mode',
        title='Modes',
        description='Note: all dates are defined with UTC boundaries.\n'
                    'The valid modes are:',
    )
    subparser = dict()

    # DAY mode
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
        default=TODAY,
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
    add_archive_options(subparser['week'])

    # MONTH mode
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
    add_archive_options(subparser['month'])

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

    add_output_options(subparser['gps'])

    # return the argument parser
    return parser


# -- main code block ----------------------------------------------------------

def main(args=None):
    """Run the GWSumm command-line interface
    """
    parser = create_parser()
    args = parser.parse_args(args=args)

    if args.debug:
        warnings.simplefilter('error', DeprecationWarning)

    # set verbose output options
    globalv.VERBOSE = args.verbose

    # find all config files
    args.config_file = [os.path.expanduser(fp) for csv in args.config_file for
                        fp in csv.split(',')]

    # check segdb option
    if args.on_segdb_error not in ['raise', 'warn', 'ignore']:
        parser.error("Invalid option --on-segdb-error='%s'" %
                     args.on_segdb_error)

    # read configuration file
    config = GWSummConfigParser()
    config.optionxform = str
    if args.ifo:
        config.set_ifo_options(args.ifo, section=DEFAULTSECT)
    config.set(DEFAULTSECT, 'user', getpass.getuser())
    config.read(args.config_file)

    try:
        ifo = config.get(DEFAULTSECT, 'IFO')
    except NoOptionError:
        ifo = None
    finally:
        globalv.IFO = ifo

    # interpolate section names
    interp = {}
    if ifo:
        interp['ifo'] = ifo.lower()
        interp['IFO'] = ifo.title()
    config.interpolate_section_names(**interp)

    # double-check week mode matches calendar setting
    if args.mode == 'week':
        if config.has_option("calendar", "start-of-week"):
            weekday = getattr(calendar,
                              config.get("calendar", "start-of-week").upper())
            if weekday != args.week.timetuple().tm_wday:
                msg = ("Cannot process week starting on %s. The "
                       "'start-of-week' option in the [calendar] section "
                       "of the INI file specifies weeks start on %ss."
                       % (args.week.strftime('%Y%m%d'),
                          config.get("calendar", "start-of-week")))
                raise parser.error(msg)

    # record times in ConfigParser
    config.set_date_options(args.gpsstart, args.gpsend, section=DEFAULTSECT)

    # convert times for convenience
    span = Segment(args.gpsstart, args.gpsend)
    utc = tconvert(args.gpsstart)
    starttime = Time(float(args.gpsstart), format='gps')
    endtime = Time(float(args.gpsend), format='gps')

    # set mode and output directory
    # This will set the output directory to be day/YYYYMMDD, week/YYYYMMDD,
    # month/YYYYMM, or gps/<gps start>-<gps end>.
    mode.set_mode(args.mode)
    try:
        path = mode.get_base(utc)
    except ValueError:
        path = os.path.join('gps', f'{args.gpsstart}-{args.gpsend}')

    # set LAL FFT plan wisdom level
    duration = min(globalv.NOW, args.gpsend) - args.gpsstart
    if duration > 200000:
        fft_lal.LAL_FFTPLAN_LEVEL = 3
    elif duration > 40000:
        fft_lal.LAL_FFTPLAN_LEVEL = 2
    else:
        fft_lal.LAL_FFTPLAN_LEVEL = 1

    # set global html only flag
    if args.html_only:
        globalv.HTMLONLY = True

    # build directories
    mkdir(args.output_dir)
    os.chdir(args.output_dir)
    plotdir = os.path.join(path, 'plots')
    mkdir(plotdir)

    # -- setup --------------------------------------

    LOGGER.info(" -- GW interferometer summary information system -- ")
    LOGGER.debug("This is process {}".format(os.getpid()))
    LOGGER.debug("You have selected {} mode".format(mode.get_mode().name))
    LOGGER.debug("Start time: {0} ({1})".format(starttime.utc.iso,
                                                starttime.gps))
    LOGGER.debug("End time: {0} ({1})".format(endtime.utc.iso,
                                              endtime.gps))
    LOGGER.debug("Output directory: {}".format(
        os.path.abspath(os.path.join(args.output_dir, path))))

    # -- Finalise configuration
    LOGGER.info("Loading configuration")
    plugins = config.load_plugins()
    if plugins:
        LOGGER.debug(" -- Loaded {} plugins:".format(len(plugins)))
        for mod in plugins:
            LOGGER.debug("        %s" % mod)
    units = config.load_units()
    LOGGER.debug("    Loaded %d units" % len(units))
    channels = config.load_channels()
    LOGGER.debug("    Loaded %d channels" % len(channels))
    states = config.load_states()
    LOGGER.debug("    Loaded %d states" % len(states))
    rcp = config.load_rcParams()
    LOGGER.debug("    Loaded %d rcParams" % len(rcp))

    # read list of tabs
    tablist = TabList.from_ini(config, match=args.process_tab,
                               path=path, plotdir=plotdir)
    tablist.sort(reverse=True)
    tabs = sorted(tablist.get_hierarchy(), key=tablist._sortkey)
    LOGGER.info("    Loaded %d tabs [%d parents overall]" % (
        len(tablist), len(tabs)))

    # read caches
    cache = {}
    for (key, var) in zip(['datacache', 'trigcache', 'segmentcache'],
                          [args.data_cache, args.event_cache,
                           args.segment_cache]):
        if var:
            LOGGER.info("Reading %s from %d files... " % (key, len(var)))
            cache[key] = Cache()
            for fp in var:
                with open(fp, 'r') as f:
                    cache[key].extend(Cache.fromfile(f))
            cache[key] = cache[key].sieve(segment=span)
            LOGGER.debug("Done [%d entries]" % len(cache[key]))

    # -- read archive -------------------------------

    if not hasattr(args, 'archive'):
        args.archive = False

    if args.html_only:
        args.archive = False
        args.daily_archive = False
    elif args.archive is True:
        args.archive = 'GW_SUMMARY_ARCHIVE'

    archives = []

    if args.archive:
        archivedir = os.path.join(path, 'archive')
        mkdir(archivedir)
        args.archive = os.path.join(archivedir, '%s-%s-%d-%d.h5'
                                    % (ifo, args.archive, args.gpsstart,
                                       args.gpsend - args.gpsstart))
        if os.path.isfile(args.archive):
            archives.append(args.archive)
        else:
            LOGGER.debug(
                "No archive found in %s, one will be created at the end"
                % args.archive)

    # read daily archive for week/month/... mode
    if hasattr(args, 'daily_archive') and args.daily_archive:
        # find daily archive files
        archives.extend(archive.find_daily_archives(
            args.gpsstart, args.gpsend, ifo, args.daily_archive, archivedir))
        # then don't read any actual data
        cache['datacache'] = Cache()

    for arch in archives:
        LOGGER.info("Reading archived data from %s" % arch)
        archive.read_data_archive(arch)
        LOGGER.debug("Archive data loaded")

    # -- read HTML configuration --------------------

    css = config.get_css(section='html')
    javascript = config.get_javascript(section='html')

    # enable comments
    try:
        globalv.HTML_COMMENTS_NAME = config.get('html', 'disqus-shortname')
    except (NoOptionError, NoSectionError):
        pass

    # find new ifo bases
    ifobases = {}
    try:
        bases_ = config.nditems('html')
    except NoSectionError:
        pass
    else:
        base_reg = re.compile(r'-base\Z')
        for key, val in bases_:
            if base_reg.search(key):
                ifobases[key.rsplit('-', 1)[0]] = val
    ifobases = OrderedDict(sorted(ifobases.items(), key=lambda x: x[0]))

    # -- write auxiliary pages ----------------------

    # get URL from output directory
    if 'public_html' in os.getcwd():
        urlbase = os.path.sep + os.path.join(
                      '~%s' % config.get(DEFAULTSECT, 'user'),
                      os.getcwd().split('public_html', 1)[1][1:])
        base = urlbase
    # otherwise get URL from html config
    elif ifo in ifobases:
        urlbase = urlparse(ifobases[ifo]).path
        base = urlbase
    # otherwise let the write_html processor work it out on-the-fly
    else:
        urlbase = None
        base = None

    # get link to issues report page
    try:
        issues = config.get('html', 'issues')
    except (NoSectionError, KeyError):
        issues = True

    # write 404 error page
    if not args.no_htaccess and not args.no_html and urlbase:
        top = os.path.join(urlbase, path)
        four0four = get_tab('404')(span=span, parent=None, path=path,
                                   index=os.path.join(path, '404.html'))
        four0four.write_html(css=css, js=javascript, tabs=tabs, ifo=ifo,
                             ifomap=ifobases, top=top, base=base,
                             writedata=not args.html_only,
                             writehtml=not args.no_html,
                             issues=issues)
        url404 = os.path.join(urlbase, four0four.index)
        with open(os.path.join(path, '.htaccess'), 'w') as htaccess:
            print('Options -Indexes', file=htaccess)
            print('ErrorDocument 404 %s' % url404, file=htaccess)
            print('ErrorDocument 403 %s' % url404, file=htaccess)

    # write config page
    about = get_tab('about')(span=span, parent=None, path=path)
    if not args.no_html:
        mkdir(about.path)
        about.write_html(
            css=css, js=javascript, tabs=tabs, config=config.files,
            prog=PROG, ifo=ifo, ifomap=ifobases, about=about.index,
            base=base, issues=issues, writedata=not args.html_only,
            writehtml=not args.no_html)

    # -- read bulk data -----------------------------

    # XXX: bulk data reading could optimise things
    # XXX: but has never been used, so should remove (DMM 18/01/16)
    if args.bulk_read and not args.html_only:
        LOGGER.info("Reading all data in BULK")
        allts = set()
        allsv = set()
        allflags = set()
        for tab in tablist:
            snames = []
            for state in tab.states:
                snames.append(state.name)
                if state.definition:
                    allflags.update(re_flagdiv.split(state.definition))
            # get all data defined for the 'all' state
            if ALLSTATE in snames:
                allts.update(tab.get_channels('timeseries', 'spectrogram',
                                              'spectrum', 'histogram'))
                allsv.update(tab.get_channels('statevector'))
                allflags.update(tab.get_flags('segments'))
            # or get data for plots defined over all states
            else:
                for plot in tab.plots:
                    if plot.state is not None:
                        continue
                    if plot.type in ['timeseries', 'spectrogram', 'spectrum',
                                     'histogram']:
                        allts.update(plot.channels)
                    elif plot.type in ['statevector']:
                        allsv.update(plot.channels)
                    elif plot.type in ['segments']:
                        allflags.update([f for cflag in plot.flags for f in
                                         re_flagdiv.split(cflag)[::2] if f])
        allseg = SegmentList([span])
        if len(allflags):
            LOGGER.info(
                "%d data-quality flags identified for segment query from all "
                "tabs" % len(allflags))
            get_segments(allflags, allseg, config=config, return_=False)
        if len(allts):
            LOGGER.info("%d channels identified for TimeSeries from all tabs"
                        % len(allts))
            get_timeseries_dict(allts, allseg,
                                config=config, nds=args.nds,
                                nproc=args.multiprocess, return_=False)
        if len(allsv):
            LOGGER.info("%d channels identified for StateVector from all tabs"
                        % len(allsv))
            get_timeseries_dict(allsv, allseg,
                                config=config, nds=args.nds, statevector=True,
                                nproc=args.multiprocess, return_=False)

    # -- process all tabs ---------------------------

    # TODO: consider re-working this loop as TabList.process_all

    for tab in tablist:
        if tab.parent:
            name = '%s/%s' % (tab.parent.name, tab.name)
        else:
            name = tab.name
        if not args.html_only and isinstance(tab, get_tab('_processed')):
            LOGGER.debug("Processing %s" % name)
            tab.process(config=config, nds=args.nds,
                        nproc=args.multiprocess,
                        segdb_error=args.on_segdb_error,
                        datafind_error=args.on_datafind_error, **cache)
        if not tab.hidden and not isinstance(tab, get_tab('link')):
            mkdir(tab.href)
            tab.write_html(
                css=css, js=javascript, tabs=tabs, ifo=ifo,
                ifomap=ifobases, about=about.index, base=base,
                issues=issues, writedata=not args.html_only,
                writehtml=not args.no_html)

        # archive this tab
        if args.archive:
            LOGGER.info("Writing data to archive")
            archive.write_data_archive(args.archive)
            LOGGER.debug("Archive written to {}".format(
                os.path.abspath(args.archive)))
        LOGGER.debug("%s complete" % (name))

    LOGGER.info("-- Data products written, all done --")


# -- run from command-line ----------------------------------------------------

if __name__ == "__main__":
    main()
