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

"""Plot the triggers for a given ETG and a given channel
"""

import argparse
import os
import sys

from collections import OrderedDict
from numpy import ndarray

from glue.lal import Cache

from gwpy.segments import Segment
from gwpy.time import to_gps

from gwdetchar.cli import logger

from ...segments import get_segments
from ...triggers import (get_triggers, keep_in_segments)
from ...utils import safe_eval

# set matplotlib backend
from matplotlib import use
use('Agg')

# backend-dependent imports
from matplotlib.artist import setp  # noqa: E402
from matplotlib.colors import LogNorm  # noqa: E402
from ...plot import get_plot  # noqa: E402

__author__ = 'Duncan Macleod <duncan.macleod@ligo.org>'
__credits__ = 'Alex Urban <alexander.urban@ligo.org>'

PROG = ('python -m gwsumm.plot.triggers' if sys.argv[0].endswith('.py')
        else os.path.basename(sys.argv[0]))
LOGGER = logger(name=PROG.split('python -m ').pop())

# default columns for plot
DEFAULT_COLUMNS = ['time', 'frequency', 'snr']


# -- parse command-line -------------------------------------------------------

def create_parser():
    """Create a command-line parser for this entry point
    """
    # initialize argument parser
    parser = argparse.ArgumentParser(
        prog=PROG,
        description=__doc__,
    )
    trigopts = parser.add_argument_group('trigger options')
    popts = parser.add_argument_group('plot options')

    # required arguments
    parser.add_argument('channel')
    parser.add_argument('gpsstart', type=to_gps)
    parser.add_argument('gpsend', type=to_gps)

    # options for reading triggers
    trigopts.add_argument(
        '-e',
        '--etg',
        default='omicron',
        help="name of ETG (default: %(default)s)",
    )
    trigopts.add_argument(
        '-l',
        '--cache-file',
        help="cache file containing event trigger file references",
    )
    trigopts.add_argument(
        '-C',
        '--columns',
        type=lambda x: x.split(','),
        help="comma-separated list of columns to read from "
             "files (default: {0})".format(','.join(DEFAULT_COLUMNS)),
    )
    trigopts.add_argument(
        '-s',
        '--snr',
        default=0,
        type=float,
        help="minimum SNR (default: %(default)s)",
    )
    trigopts.add_argument(
        '-a',
        '--state',
        metavar='FLAG',
        help="restrict triggers to active times for flag",
    )

    # options for output plot
    popts.add_argument(
        '-t',
        '--epoch',
        type=to_gps,
        help="Zero-time for plot (default: gpsstart)",
    )
    popts.add_argument(
        '-x',
        '--x-column',
        help="column for X-axis (default: first --column)",
    )
    popts.add_argument(
        '-f',
        '--y-column',
        help="column for Y-axis (default: second --column)",
    )
    popts.add_argument(
        '-c',
        '--color',
        help="column for colour (default: third --column)",
    )
    popts.add_argument(
        '-p',
        '--plot-params',
        action='append',
        metavar='"{arg}={param}"',
        help="extra plotting keyword arguments, "
             "can be given multiple times",
    )
    popts.add_argument(
        '--tiles',
        action='store_true',
        default=False,
        help="plot tiles instead of dots (default: %(default)s), "
             "if selected 'duration' and 'bandwidth' will be "
             "automatically added to --columns",
    )
    popts.add_argument(
        '--state-label',
        help="label for state segments (default: --state)",
    )
    popts.add_argument(
        '-o',
        '--output-file',
        default='triggers.png',
        help="output file name (default: %(default)s)",
    )

    # return the argument parser
    return parser


# -- main code block ----------------------------------------------------------

def main(args=None):
    """Run the online trigger visualization tool
    """
    parser = create_parser()
    args = parser.parse_args(args=args)

    args.epoch = args.epoch or args.gpsstart
    args.columns = args.columns or DEFAULT_COLUMNS
    if len(args.columns) < 2:
        parser.error("--columns must receive at least two columns, "
                     "got {0}".format(len(args.columns)))
    args.x_column = args.x_column or args.columns[0]
    args.y_column = args.y_column or args.columns[1]
    if not args.color and len(args.columns) >= 3:
        args.color = args.columns[2]

    # add columns for tile plot
    for c in ('duration', 'bandwidth'):
        if args.tiles and c not in args.columns:
            args.columns.append(c)

    span = Segment(args.gpsstart, args.gpsend)

    # format default params
    params = {
        'xscale': 'auto-gps',
        'epoch': args.epoch or args.gpsstart,
        'xlim': (args.gpsstart, args.gpsend),
        'yscale': 'log',
        'ylabel': 'Frequency [Hz]',
        'cmap': get_plot('triggers').defaults.get('cmap', 'YlGnBu'),
        'clim': (3, 50),
        'colorlabel': 'Signal-to-noise ratio (SNR)',
        'sortbycolor': True,
    }

    # update with user params
    for input_ in args.plot_params or []:
        key, val = input_.split('=', 1)
        params[key.strip('-')] = safe_eval(val)

    # -- load triggers ------------------------------

    # get segments
    if args.state:
        segs = get_segments(args.state, [span])

    # read cache
    if args.cache_file:
        with open(args.cache_file, 'r') as f:
            cache = Cache.fromfile(f).sieve(segment=span)
        LOGGER.info("Read cache of {0} files".format(len(cache)))
    else:
        cache = None

    # get triggers
    trigs = get_triggers(
        args.channel, args.etg, [span], cache=cache, columns=args.columns,
        verbose='Reading {0.channel} ({0.etg}) events'.format(args))
    LOGGER.info("Read {0} events".format(len(trigs)))
    if args.state:
        trigs = keep_in_segments(trigs, segs.active, etg=args.etg)
        LOGGER.info("{0} events in state {1!r}".format(len(trigs), args.state))
    if args.snr:
        trigs = trigs[trigs['snr'] >= args.snr]
        LOGGER.info("{0} events remaining with snr >= {1.snr}".format(
            len(trigs), args))

    # -- make plot ----------------------------------

    # format keywords for plot creation
    plot_kw = OrderedDict(
        (key, params.pop(key))
        for key in (
            'xscale',
            'xlim',
            'epoch',
            'yscale',
            'ylabel',
            'sortbycolor',
        )
    )

    # create plot
    if args.tiles:
        plot = trigs.tile(args.x_column, args.y_column,
                          'duration', 'bandwidth', color=args.color,
                          edgecolor=params.pop('edgecolor', 'face'),
                          linewidth=params.pop('linewidth', .8), **plot_kw)
    else:
        plot = trigs.scatter(args.x_column, args.y_column, color=args.color,
                             edgecolor=params.pop('edgecolor', 'none'),
                             s=params.pop('s', 12), **plot_kw)
    ax = plot.gca()
    mappable = ax.collections[0]

    # set mappable properties
    vmin, vmax = params.pop('clim', (3, 50))
    if params.pop('logcolor', True):
        mappable.set_norm(LogNorm(vmin=vmin, vmax=vmax))
    if mappable._A is None:
        mappable._A = ndarray((0,))

    # draw colorbar
    clabel = params.pop(
        'colorlabel',
        'Signal-to-noise ratio (SNR)',
    )
    cmap = params.pop(
        'cmap',
        get_plot('triggers').defaults.get('cmap',  'YlGnBu'),
    )
    ax.colorbar(mappable=mappable, label=clabel, cmap=cmap, clim=(vmin, vmax))

    for key, val in params.items():
        try:
            getattr(ax, 'set_%s' % key)(val)
        except AttributeError:
            setattr(ax, key, val)

    # add segments
    if args.state:
        sax = plot.add_segments_bar(segs, ax=ax,
                                    label=args.state_label or args.state)
        setp(sax.get_yticklabels(), fontsize=10)
        sax.set_epoch(args.epoch)

    # save and exit
    plot.save(args.output_file)
    LOGGER.info("Plot saved to {0.output_file}".format(args))


# -- run from command-line ----------------------------------------------------

if __name__ == "__main__":
    main()
