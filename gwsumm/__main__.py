# coding=utf-8
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

"""Command-line interface for the gravitational-wave interferometer summary
information system (`gwsumm`).

This module provides the command-line interface to the GWSumm package,
allowing generation of detector summary information.

Select a <mode> to run over a calendar amount of time ('day', 'week',
or 'month'), or an arbitrary GPS (semi-open) interval.

Run '%(prog)s <mode> --help' for details of the specific arguments and
options available for each mode.
"""

import calendar
import getpass
import os
import re
import sys
import warnings

from astropy.config import get_config_dir
from astropy.config.paths import get_cache_dir

from collections import OrderedDict
from configparser import (DEFAULTSECT, NoOptionError, NoSectionError)
from matplotlib import font_manager as fm
from urllib.parse import urlparse

from glue.lal import (Cache, LIGOTimeGPS)

from gwpy import time
from gwpy.segments import (Segment, SegmentList)
from gwpy.signal.spectral import _lal as fft_lal
from gwpy.time import (Time, _tconvert, tconvert)

from gwdetchar.utils.cli import logger

from . import (
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
    re_flagdiv,
    create_parser,
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
__credits__ = "Evan Goetz <evan.goetz@ligo.org>"

# set defaults
VERBOSE = False
PROFILE = False

# initialize logger
PROG = ('python -m gwsumm' if sys.argv[0].endswith('.py')
        else os.path.basename(sys.argv[0]))
LOGGER = logger(name=PROG.split('python -m ').pop())


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
        if args.archive_write_dir is None:
            args.archive_write_dir = f'{path}/archive'
        if args.archive_read_dir is None:
            args.archive_read_dir = f'{path}/archive'

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

    # build output directories
    os.makedirs(args.output_dir, exist_ok=True)
    os.chdir(args.output_dir)
    plotdir = os.path.join(path, 'plots')
    os.makedirs(plotdir, exist_ok=True)

    # build astropy XDG directories
    # The functions get_config_dir() and get_cache_dir() either return the
    # default ~/.astropy/config, ~/.astropy/cache directories or the
    # directories specified by the XDG_CONFIG_HOME and XDG_CACHE_HOME
    # environment variables. We only need to create directories if the
    # variables are set
    if os.environ.get('XDG_CONFIG_HOME', None):
        os.makedirs(os.path.join(get_config_dir(), 'astropy'), exist_ok=True)
    if os.environ.get('XDG_CACHE_HOME', None):
        os.makedirs(os.path.join(get_cache_dir(), 'astropy'), exist_ok=True)

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

    # Set fonts
    if ('XDG_DATA_HOME' in os.environ and
            os.path.exists(
                p := os.path.join(os.environ['XDG_DATA_HOME'], 'fonts')
            )):
        LOGGER.debug(f"Adding fonts for matplotlib from {p}")
        font_files = fm.findSystemFonts(fontpaths=p)
        for f in font_files:
            fm.fontManager.addfont(f)

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

    # EG: I don't know why this would be needed. Commenting out for now
    # if not hasattr(args, 'archive'):
    #     args.archive = False

    if args.html_only:
        args.archive = False
        args.daily_archive = False
    elif args.archive is True:
        args.archive = 'GW_SUMMARY_ARCHIVE'

    archives_read = []
    archives_write = []

    if args.archive:
        os.makedirs(args.archive_write_dir, exist_ok=True)
        archive_file_read = os.path.join(
            args.archive_read_dir,
            (f'{ifo}-{args.archive}-{args.gpsstart}-'
             f'{args.gpsend - args.gpsstart}.h5')
        )
        archive_file_write = os.path.join(
            args.archive_write_dir,
            (f'{ifo}-{args.archive}-{args.gpsstart}-'
             f'{args.gpsend - args.gpsstart}.h5')
        )
        if not os.path.isfile(archive_file_read):
            LOGGER.debug(f"No archive found in {archive_file_read}, one will"
                         f" be created at {archive_file_write}")
        else:
            archives_read.append(archive_file_read)
        archives_write.append(archive_file_write)

    # read daily archive for week/month/... mode
    if hasattr(args, 'daily_archive') and args.daily_archive:
        # find daily archive files
        archives_read.extend(archive.find_daily_archives(
            args.gpsstart,
            args.gpsend,
            ifo,
            args.daily_archive,
            args.archive_read_dir,
        ))
        # then don't read any actual data
        cache['datacache'] = Cache()

    for arch in archives_read:
        LOGGER.info(f"Reading archived data from {arch}")
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
        os.makedirs(about.path, exist_ok=True)
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
            os.makedirs(tab.href, exist_ok=True)
            tab.write_html(
                css=css, js=javascript, tabs=tabs, ifo=ifo,
                ifomap=ifobases, about=about.index, base=base,
                issues=issues, writedata=not args.html_only,
                writehtml=not args.no_html)

        # archive this tab
        if args.archive:
            LOGGER.info("Writing data to archive")
            try:
                archive.write_data_archive(archives_write[0])
                LOGGER.debug(
                    f"Archive written to {os.path.abspath(archives_write[0])}")
            except Exception:
                LOGGER.warning(
                    "New data archiving failed. Previous archive preserved.")
        LOGGER.debug("%s complete" % (name))

    LOGGER.info("-- Data products written, all done --")


# -- run from command-line ----------------------------------------------------

if __name__ == "__main__":
    main()
