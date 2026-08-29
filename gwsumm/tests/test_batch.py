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

"""Tests for the `gwsumm.batch` command-line interface
"""

import os
from pathlib import Path
import pytest

from .. import batch

__author__ = 'Alex Urban <alexander.urban@ligo.org>'
__credits__ = 'Evan Goetz <evan.goetz@ligo.org>'


# -- utilities ----------------------------------------------------------------

def _get_inputs():
    """Prepare and return paths to input data products
    """
    inputs = (
        "global.ini",
        "k1-test.ini",
        "K1-TEST-1.h5"
    )
    # write empty input files
    for filename in inputs:
        Path(filename).touch()
    return inputs


# -- cli tests ----------------------------------------------------------------

def test_main(tmpdir, caplog):
    outdir = str(tmpdir)
    (global_, k1test, archivefile) = _get_inputs()
    args = [
        'day',
        '--verbose',
        '--ifo', 'K1',
        '--maxjobs', '5',
        '--condor-timeout', '3',
        '--condor-command', 'notification=false',
        '--condor-command', 'environment=VAR1=VAL1',
        '--condor-command', 'environment=VAR2=VAL2',
        '--config-file', k1test,
        '--global-config', global_,
        '--nds',
        '--multi-process', '4',
        '--archive',
        '--archive-read-dir', str(Path('.').absolute()),
        '--event-cache', '/this/cache/is/not/real.cache',
        '--no-htaccess',
        '--output-dir', outdir,
        '--container-path', str(Path('./fake_container.sif').absolute()),
    ]
    # test log output
    batch.main(args)
    assert "Copied all INI configuration files to ./etc" in caplog.text
    assert "Appending VAR2=VAL2 to condor 'environment' value" in caplog.text
    assert " -- Configured HTML htmlnode job" in caplog.text
    assert " -- Configured job for config {}".format(
        os.path.join(outdir, "etc", os.path.basename(k1test))) in caplog.text
    assert "Setup complete, DAG written to: {}".format(
        os.path.join(outdir, "gw_summary_pipe.dag")) in caplog.text
    # test file output
    assert set(os.listdir(outdir)) == {
        "etc",
        "gw_summary_pipe_local.sub",
        "gw_summary_pipe.dag",
        "logs",
        "gw_summary_pipe.sub",
        "gw_summary_pipe.sh",
    }
    assert set(os.listdir(os.path.join(outdir, "etc"))) == {
        os.path.basename(k1test),
        os.path.basename(global_),
    }
    assert set(os.listdir(os.path.join(outdir, "logs"))) == set()
    # clean up
    for filename in (global_, k1test, archivefile):
        Path(filename).unlink()
    for filename in Path(outdir, 'etc').iterdir():
        Path(filename).unlink()
    Path(outdir, 'etc').rmdir()
    Path(outdir, 'logs').rmdir()
    for filename in Path(outdir).iterdir():
        Path(filename).unlink()
    Path(outdir).rmdir()


@pytest.mark.parametrize(
    'mode',
    (['day', '20170209'],
     ['week', '20170209'],
     ['month', '201702'],
     ['gps', '1170633618', '1170720018']),
)
def test_main_loop_over_modes(tmpdir, caplog, mode):
    outdir = str(tmpdir)
    (global_, k1test, archivefile) = _get_inputs()
    args = [
        '--verbose',
        '--ifo', 'K1',
        '--universe', 'local',
        '--config-file', k1test,
        '--global-config', global_,
        '--output-dir', outdir,
        '--container-path', str(Path('./fake_container.sif').absolute()),
    ]
    # test log output
    batch.main(mode + args)
    assert "Copied all INI configuration files to ./etc" in caplog.text
    assert " -- Configured HTML htmlnode job" in caplog.text
    assert " -- Configured job for config {}".format(
        os.path.join(outdir, "etc", os.path.basename(k1test))) in caplog.text
    assert "Setup complete, DAG written to: {}".format(
        os.path.join(outdir, "gw_summary_pipe.dag")) in caplog.text
    # clean up
    for filename in (global_, k1test, archivefile):
        Path(filename).unlink()
    for filename in Path(outdir, 'etc').iterdir():
        Path(filename).unlink()
    Path(outdir, 'etc').rmdir()
    Path(outdir, 'logs').rmdir()
    for filename in Path(outdir).iterdir():
        Path(filename).unlink()
    Path(outdir).rmdir()
