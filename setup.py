#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright (C) Duncan Macleod (2013)
#
# This file is part of the GWSumm package.
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

"""Setup the GWSumm package
"""

import glob
import os.path
import versioneer

from setuptools import setup

version = versioneer.get_version()
cmdclass = versioneer.get_cmdclass()

try:
    from sphinx.setup_command import BuildDoc
except ImportError:
    pass
else:
    cmdclass["build_sphinx"] = BuildDoc

# configuration files
data_files = [
    ("etc/gwsumm/configuration",
     glob.glob(os.path.join("share", "*.ini"))),
]

# run setup
# NOTE: all other metadata and options come from setup.cfg
setup(
    version=version,
    project_urls={
         "Bug Tracker": "https://github.com/gwpy/gwsumm/issues",
         "Discussion Forum": "https://gwdetchar.slack.com",
         "Documentation": "https://gwsumm.readthedocs.io",
         "Source Code": "https://github.com/gwpy/gwsumm",
     },
    cmdclass=cmdclass,
    data_files=data_files,
)
