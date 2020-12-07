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

"""Tests for the `gwsumm.plot.triggers` command-line interface
"""

import os
import pytest
import shutil

from unittest import mock

from gwpy.segments import (
    Segment,
    SegmentList,
    DataQualityFlag,
)

from .... import globalv
from .. import __main__ as triggers_cli

__author__ = 'Alex Urban <alexander.urban@ligo.org>'

# -- test configuration

CHANNEL = "H1:GDS-CALIB_STRAIN"

# -- test data

LOCK = DataQualityFlag(
    name="H1:DMT-GRD_ISC_LOCK_NOMINAL:1",
    active=SegmentList([Segment(2, 2048)]),
    known=SegmentList([Segment(0, 3600)]),
)


# -- cli tests ----------------------------------------------------------------

@mock.patch(
    'gwpy.segments.DataQualityFlag.query_dqsegdb',
    return_value=LOCK,
)
def test_main(dqflag, tmpdir, caplog):
    outdir = str(tmpdir)
    plot = os.path.join(outdir, "triggers.png")
    args = [
        CHANNEL,
        '0', '3600',
        '--snr', '1',
        '--state', LOCK.name,
        '--output-file', plot,
    ]
    # test output
    with pytest.warns(UserWarning) as record:
        triggers_cli.main(args)
    assert os.path.exists(plot)
    assert len(os.listdir(outdir)) == 1
    assert 'Read 0 events' in caplog.text
    assert "0 events in state '{}'".format(LOCK.name) in caplog.text
    assert '0 events remaining with snr >= 1.0' in caplog.text
    assert 'Plot saved to {}'.format(plot) in caplog.text
    # test the `UserWarning`
    # FIXME: once MatplotlibDeprecationWarning about colormaps is fixed,
    #        assert that this UserWarning is the **only** warning
    assert (record[0].message.args[0] ==
            "Caught ValueError: No channel-level directory found at "
            "/home/detchar/triggers/*/H1/GDS-CALIB_STRAIN_Omicron. Either "
            "the channel name or ETG names are wrong, or this channel is not "
            "configured for this ETG.")
    # clean up
    globalv.TRIGGERS = {}
    shutil.rmtree(outdir, ignore_errors=True)


def test_main_with_cache_and_tiles(tmpdir, caplog):
    outdir = str(tmpdir)
    cache = os.path.join(outdir, "empty.cache")
    plot = os.path.join(outdir, "triggers.png")
    args = [
        CHANNEL,
        '0', '3600',
        '--cache-file', cache,
        '--snr', '1',
        '--plot-params', 'legend-loc="upper right"',
        '--tiles',
        '--output-file', plot,
    ]
    # write an empty cache file
    with open(cache, 'w') as f:
        f.write("")
    # test output
    triggers_cli.main(args)
    assert os.path.exists(plot)
    assert len(os.listdir(outdir)) == 2  # 1 input, 1 output
    assert 'Read cache of 0 files' in caplog.text
    assert 'Read 0 events' in caplog.text
    assert '0 events remaining with snr >= 1.0' in caplog.text
    assert 'Plot saved to {}'.format(plot) in caplog.text
    # clean up
    globalv.TRIGGERS = {}
    shutil.rmtree(outdir, ignore_errors=True)


def test_main_invalid_columns(capsys):
    args = [
        CHANNEL,
        '0', '3600',
        '--columns', 'invalid',
    ]
    # test output
    with pytest.raises(SystemExit):
        triggers_cli.main(args)
    (_, err) = capsys.readouterr()
    assert err.endswith("--columns must receive at least two columns, got 1\n")
