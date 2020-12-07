# -*- coding: utf-8 -*-
# Copyright (C) Alex Urban (2020)
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

"""Tests for the `gwsumm.plot.guardian` command-line interface
"""

import os
import shutil

from gwpy.timeseries import (
    TimeSeries,
    TimeSeriesList,
)

from .... import globalv
from ....archive import write_data_archive
from .. import __main__ as guardian_cli

__author__ = 'Alex Urban <alexander.urban@ligo.org>'

# -- test configuration

CONFIG = """
[tab-ISC_LOCK]
type = guardian
node = ISC_LOCK
name = %(node)s
; node states
600 = Low noise
"""

# -- test data

SUFFICES = [
    "STATE_N",
    "REQUEST_N",
    "NOMINAL_N",
    "OK",
    "MODE",
    "OP",
]
DATA = {
    key: TimeSeriesList(
        TimeSeries([600] * 3600 * 16, sample_rate=16, name=key, channel=key)
    ) for key in ["L1:GRD-ISC_LOCK_{}".format(suff) for suff in SUFFICES]
}


# -- utils --------------------------------------------------------------------

def _get_inputs(workdir):
    """Prepare, and return paths to, input data products
    """
    # set global timeseries data
    globalv.DATA = DATA
    # get path to data files
    ini = os.path.join(workdir, "config.ini")
    archive = os.path.abspath(os.path.join(workdir, "archive.h5"))
    # write to data files
    with open(ini, 'w') as f:
        f.write(CONFIG)
    write_data_archive(archive)
    # reset global data and return
    globalv.DATA = {}
    return (ini, archive)


# -- cli tests ----------------------------------------------------------------

def test_main(tmpdir, caplog):
    outdir = str(tmpdir)
    plot = os.path.join(outdir, "guardian.png")
    (ini, archive) = _get_inputs(outdir)
    args = [
        'ISC_LOCK',
        '0', '3600',
        ini,
        '--plot-params', 'title=Test figure',
        '--output-file', plot,
        '--verbose',
        '--archive', archive,
    ]
    # test output
    guardian_cli.main(args)
    assert os.path.exists(plot)
    assert len(os.listdir(outdir)) == 3  # 2 inputs, 1 output
    assert 'Read data archive from {}'.format(archive) in caplog.text
    assert 'Processing:' in caplog.text
    assert 'Plot saved to {}'.format(plot) in caplog.text
    assert 'Archive recorded as {}'.format(archive) in caplog.text
    # clean up
    shutil.rmtree(outdir, ignore_errors=True)
