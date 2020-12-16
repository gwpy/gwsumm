# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2018)
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

"""Plot the segments for a given Guardian node with a given definition
"""

import argparse
import os
import re
import shutil
import sys

from collections import OrderedDict
from configparser import DEFAULTSECT

from gwpy.time import to_gps

from gwdetchar.cli import logger

from ... import globalv
from ...archive import (write_data_archive, read_data_archive)
from ...config import GWSummConfigParser
from ...data import get_timeseries
from ...state import generate_all_state
from ...tabs import GuardianTab

# set matplotlib backend
from matplotlib import use
use('Agg')

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__credits__ = 'Alex Urban <alexander.urban@ligo.org>'

PROG = ('python -m gwsumm.plot.guardian' if sys.argv[0].endswith('.py')
        else os.path.basename(sys.argv[0]))
LOGGER = logger(name=PROG.split('python -m ').pop())

GWSummConfigParser.OPTCRE = re.compile(
    r'(?P<option>[^=\s][^=]*)\s*(?P<vi>[=])\s*(?P<value>.*)$')


# -- utilities ----------------------------------------------------------------

def safe_eval(val):
    """Evaluate the given string as a line of python, if possible

    If the :meth:`eval` fails, a `str` is returned instead.
    """
    try:
        return eval(val)
    except (NameError, SyntaxError):
        return str(val)


# -- parse command-line -------------------------------------------------------

def create_parser():
    """Create a command-line parser for this entry point
    """
    # initialize argument parser
    parser = argparse.ArgumentParser(
        prog=PROG,
        description=__doc__,
    )
    archopts = parser.add_mutually_exclusive_group()

    # positional arguments
    parser.add_argument('node')
    parser.add_argument('gpsstart', type=to_gps)
    parser.add_argument('gpsend', type=to_gps)
    parser.add_argument('config', help='config file defining Guardian node')

    # optional flags
    parser.add_argument(
        '-i',
        '--ifo',
        type=str,
        default="L1",
    )
    parser.add_argument(
        '-s',
        '--section',
        type=str,
        help="suffix of INI tab section to read, e.g. give "
             "--section='ISC_LOCK' to read [tab-ISC_LOCK] "
             "section, defaults to {node}",
    )
    parser.add_argument(
        '-t',
        '--epoch',
        type=to_gps,
        help="Zero-time for plot, defaults to GPSSTART",
    )
    parser.add_argument(
        '-p',
        '--plot-params',
        action='append',
        default=[],
        help="extra plotting keyword argument",
    )
    parser.add_argument(
        '-m',
        '--multi-process',
        type=int,
        default=1,
        dest='nproc',
        help="number of processes to use, default: %(default)s",
    )
    parser.add_argument(
        '-o',
        '--output-file',
        default="guardian.png",
        help="output file name, default: %(default)s",
    )
    parser.add_argument(
        '-v',
        '--verbose',
        action='store_true',
        help="print verbose output, default: False",
    )
    parser.add_argument(
        '-P',
        '--profile',
        action='store_true',
        help="print timing output, default: False",
    )

    # archive options
    archopts.add_argument(
        '-a',
        '--archive',
        help="full path of HDF archive for data",
    )
    archopts.add_argument(
        '-r',
        '--read-only-archive',
        help="full path of HDF archive for data, does not write",
    )

    # return the argument parser
    return parser


# -- main code block ----------------------------------------------------------

def main(args=None):
    """Run the online Guardian node visualization tool
    """
    parser = create_parser()
    args = parser.parse_args(args=args)

    globalv.VERBOSE = args.verbose
    globalv.PROFILE = args.profile
    args.epoch = args.epoch or args.gpsstart
    state = generate_all_state(args.gpsstart, args.gpsend)

    # format params
    params = {}
    for input_ in args.plot_params:
        key, val = input_.split('=', 1)
        params[key.strip('-')] = safe_eval(val)

    # read config
    config = GWSummConfigParser(dict_type=OrderedDict)
    config.read([args.config])
    config.set(DEFAULTSECT, 'gps-start-time', str(int(args.gpsstart)))
    config.set(DEFAULTSECT, 'gps-end-time', str(int(args.gpsend)))
    config.set(DEFAULTSECT, 'IFO', args.ifo)
    sec = 'tab-{}'.format(args.section or args.node)

    # read archive
    if args.archive and not args.read_only_archive:
        args.read_only_archive = args.archive
    if args.read_only_archive and os.path.isfile(args.read_only_archive):
        read_data_archive(args.read_only_archive)
        LOGGER.info(
            "Read data archive from {0.read_only_archive}".format(args))

    # make tab
    tab = GuardianTab.from_ini(config, sec, mode='gps', path='.', plotdir='.')
    tab.plots = tab.plots[:1]
    tab.plots[0].pargs.update(params)
    tab.plots[0].pargs['epoch'] = args.epoch

    # process
    LOGGER.info("Processing:")
    tab.process(nproc=args.nproc)
    plotfile = tab.plots[0].outputfile
    shutil.copy(plotfile, args.output_file)
    os.remove(plotfile)
    LOGGER.info("Plot saved to {0.output_file}".format(args))

    # crop and save archive
    if args.archive:
        for channel in globalv.DATA:
            globalv.DATA[channel] = get_timeseries(channel, state, query=False)
        write_data_archive(args.archive)
        LOGGER.info("Archive recorded as {0.archive}".format(args))


# -- run from command-line ----------------------------------------------------

if __name__ == "__main__":
    main()
