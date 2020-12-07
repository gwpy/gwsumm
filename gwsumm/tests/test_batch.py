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
import pytest
import shutil

from unittest import mock

from .. import batch

__author__ = 'Alex Urban <alexander.urban@ligo.org>'


# -- utilities ----------------------------------------------------------------

def _get_inputs():
    """Prepare and return paths to input data products
    """
    indir = os.getcwd()
    inputs = (
        os.path.join(indir, "global.ini"),
        os.path.join(indir, "k1-test.ini"),
        os.path.join(indir, "x509.cert"),
    )
    # write empty input files
    for filename in inputs:
        with open(filename, 'w') as f:
            f.write("")
    return inputs


# -- cli tests ----------------------------------------------------------------

@mock.patch(
    'gwsumm.batch.find_credential',
)
@mock.patch(
    'gwpy.io.kerberos.kinit',
    return_value=None,
)
def test_main(krb, x509, tmpdir, caplog):
    outdir = str(tmpdir)
    (global_, k1test, x509cert) = _get_inputs()
    x509.return_value = (x509cert, x509cert)
    args = [
        '--verbose',
        '--ifo', 'K1',
        '--maxjobs', '5',
        '--condor-timeout', '3',
        '--condor-command', 'notification=false',
        '--config-file', k1test,
        '--global-config', global_,
        '--nds',
        '--multi-process', '4',
        '--archive',
        '--event-cache', '/this/cache/is/not/real.cache',
        '--no-htaccess',
        '--output-dir', outdir,
    ]
    # test log output
    batch.main(args)
    assert "Copied all INI configuration files to ./etc" in caplog.text
    assert ("Configured Condor and Kerberos for NFS-shared credentials"
            in caplog.text)
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
        os.path.basename(x509cert),
        os.path.basename(global_),
    }
    assert set(os.listdir(os.path.join(outdir, "logs"))) == set()
    # clean up
    for filename in (global_, k1test, x509cert):
        os.remove(filename)
    shutil.rmtree(outdir, ignore_errors=True)


@pytest.mark.parametrize(
    'mode',
    (['--day', '20170209'],
     ['--week', '20170209'],
     ['--month', '201702'],
     ['--year', '2017'],
     ['--gps-start-time', '1170633618', '--gps-end-time', '1170720018']),
)
def test_main_loop_over_modes(tmpdir, caplog, mode):
    outdir = str(tmpdir)
    (global_, k1test, x509cert) = _get_inputs()
    args = [
        '--verbose',
        '--ifo', 'K1',
        '--universe', 'local',
        '--config-file', k1test,
        '--global-config', global_,
        '--single-process',
        '--output-dir', outdir,
    ]
    # test log output
    batch.main(args + mode)
    assert "Copied all INI configuration files to ./etc" in caplog.text
    assert " -- Configured HTML htmlnode job" in caplog.text
    assert " -- Configured job for config {}".format(
        os.path.join(outdir, "etc", os.path.basename(k1test))) in caplog.text
    assert "Setup complete, DAG written to: {}".format(
        os.path.join(outdir, "gw_summary_pipe.dag")) in caplog.text
    # clean up
    for filename in (global_, k1test, x509cert):
        os.remove(filename)
    shutil.rmtree(outdir, ignore_errors=True)


def test_main_invalid_modes(capsys):
    args = [
        '--ifo', 'V1',
        '--day', '20170209',
        '--month', '201702',
    ]
    # test output
    with pytest.raises(SystemExit):
        batch.main(args)
    (_, err) = capsys.readouterr()
    assert err.endswith("Please give only one of --day, --month, or "
                        "--gps-start-time and --gps-end-time.\n")
