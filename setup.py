#!/usr/bin/env python
# -*- coding: utf-8 -*-
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

"""Setup the GWSumm package
"""

from __future__ import print_function

import sys
if sys.version < '2.6':
    raise ImportError("Python versions older than 2.6 are not supported.")

import glob
import hashlib
import os.path
import subprocess

try:
    import setuptools
except ImportError:
    import ez_setup
    ez_setup.use_setuptools()
finally:
    from setuptools import (setup, find_packages)
    from setuptools.command import (build_py, egg_info)

from distutils import log
from distutils.dist import Distribution
from distutils.cmd import Command
from distutils.command.clean import (clean, log, remove_tree)

# test for OrderedDict
extra_install_requires = []
try:
    from collections import OrderedDict
except ImportError:
    extra_install_requires.append('ordereddict>=1.1')

# import sphinx commands
try:
    from sphinx.setup_command import BuildDoc
except ImportError:
    cmdclass = {}
else:
    cmdclass = {'build_sphinx': BuildDoc}

# set basic metadata
PACKAGENAME = 'gwsumm'
AUTHOR = 'Duncan Macleod'
AUTHOR_EMAIL = 'duncan.macleod@ligo.org'
LICENSE = 'GPLv3'

VERSION_PY = os.path.join(PACKAGENAME, 'version.py')


# -----------------------------------------------------------------------------
# Clean up, including Sphinx, and setup_requires eggs

class GWSummClean(clean):
    def run(self):
        if self.all:
            # remove docs
            sphinx_dir = os.path.join(self.build_base, 'sphinx')
            if os.path.exists(sphinx_dir):
                remove_tree(sphinx_dir, dry_run=self.dry_run)
            else:
                log.warn("%r does not exist -- can't clean it", sphinx_dir)
            # remove version.py
            for vpy in [VERSION_PY, VERSION_PY + 'c']:
                if os.path.exists(vpy) and not self.dry_run:
                    log.info('removing %r' % vpy)
                    os.unlink(vpy)
                elif not os.path.exists(vpy):
                    log.warn("%r does not exist -- can't clean it", vpy)
            # remove setup eggs
            for egg in glob.glob('*.egg'):
                if os.path.isdir(egg):
                    remove_tree(egg, dry_run=self.dry_run)
                else:
                    log.info('removing %r' % egg)
                    os.unlink(egg)
            # remove Portfile
            portfile = 'Portfile'
            if os.path.exists(portfile) and not self.dry_run:
                log.info('removing %r' % portfile)
                os.unlink(portfile)
        clean.run(self)

cmdclass['clean'] = GWSummClean


# -----------------------------------------------------------------------------
# Custom builders to write version.py

class GitVersionMixin(object):
    """Mixin class to add methods to generate version information from git.
    """
    def write_version_py(self, pyfile):
        """Generate target file with versioning information from git VCS
        """
        log.info("generating %s" % pyfile)
        import vcs
        gitstatus = vcs.GitStatus()
        try:
            with open(pyfile, 'w') as fobj:
                gitstatus.write(fobj, author=AUTHOR, email=AUTHOR_EMAIL)
        except:
            if os.path.exists(pyfile):
                os.unlink(pyfile)
            raise
        return gitstatus

    def update_metadata(self):
        """Import package base and update distribution metadata
        """
        import gwsumm
        self.distribution.metadata.version = gwsumm.__version__
        desc, longdesc = gwsumm.__doc__.split('\n', 1)
        self.distribution.metadata.description = desc
        self.distribution.metadata.long_description = longdesc.strip('\n')


class GWSummBuildPy(build_py.build_py, GitVersionMixin):
    """Custom build_py command to deal with version generation
    """
    def __init__(self, *args, **kwargs):
        build_py.build_py.__init__(self, *args, **kwargs)

    def run(self):
        try:
            self.write_version_py(VERSION_PY)
        except ImportError:
            raise
        except:
            if not os.path.isfile(VERSION_PY):
                raise
        self.update_metadata()
        build_py.build_py.run(self)

cmdclass['build_py'] = GWSummBuildPy


class GWSummEggInfo(egg_info.egg_info, GitVersionMixin):
    """Custom egg_info command to deal with version generation
    """
    def finalize_options(self):
        try:
            self.write_version_py(VERSION_PY)
        except ImportError:
            raise
        except:
            if not os.path.isfile(VERSION_PY):
                raise
        if not self.distribution.metadata.version:
            self.update_metadata()
        egg_info.egg_info.finalize_options(self)

cmdclass['egg_info'] = GWSummEggInfo


# -----------------------------------------------------------------------------
# Process complicated dependencies

# XXX: this can be removed as soon as a stable release of glue can
#      handle pip/--user
try:
    from glue import git_version
except ImportError as e:
    e.args = ("GWSumm requires the GLUE package, which isn\'t available from "
              "PyPI.\nPlease visit\n"
              "https://www.lsc-group.phys.uwm.edu/daswg/projects/glue.html\n"
              "to download and install it manually.",)
    raise

# don't use setup_requires if just checking for information
# (credit: matplotlib/setup.py)
setup_requires = []
if not '--help' in sys.argv:
    dist_ = Distribution({'cmdclass': cmdclass})
    dist_.parse_config_files()
    dist_.parse_command_line()
    if not (any('--' + opt in sys.argv for opt in
            Distribution.display_option_names + ['help']) or
            dist_.commands == ['clean']):
        setup_requires = ['tornado', 'numpy >= 1.5', 'jinja2', 'gitpython']

# -----------------------------------------------------------------------------
# Find files

# Use the find_packages tool to locate all packages and modules
packagenames = find_packages()

# glob for all scripts
if os.path.isdir('bin'):
    scripts = glob.glob(os.path.join('bin', '*'))
else:
    scripts = []

if os.path.isdir('share'):
    configs = glob.glob(os.path.join('share', '*.ini'))
else:
    configs = []

# get all HTML files
if os.path.isdir('share'):
    htmlfiles = (glob.glob(os.path.join('share', '*.css')) +
                 glob.glob(os.path.join('share', '*.js')))
else:
    htmlfiles = []

# -----------------------------------------------------------------------------
# run setup

setup(name=PACKAGENAME,
      provides=[PACKAGENAME],
      version=None,
      description=None,
      long_description=None,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      url='https://gwpy.github.io/gwsumm',
      packages=packagenames,
      include_package_data=True,
      cmdclass=cmdclass,
      scripts=scripts,
      setup_requires=setup_requires,
      data_files=[('etc', configs),
                  (os.path.join('share', PACKAGENAME, 'html'), htmlfiles)],
      requires=[
          'glue',
          'numpy',
          'matplotlib',
          'astropy',
          'gwpy',
      ],
      install_requires=[
          'numpy >= 1.5',
          'matplotlib >= 1.3.0',
          'astropy >= 0.3',
          'gwpy >= 0.1a5',
      ] + extra_install_requires,
      dependency_links=[
          'https://www.lsc-group.phys.uwm.edu/daswg/download/'
              'software/source/glue-1.46.tar.gz#egg=glue-1.46',
      ],
      use_2to3=False,
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'Intended Audience :: End Users/Desktop',
          'Intended Audience :: Developers',
          'Natural Language :: English',
          'Topic :: Scientific/Engineering',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Physics',
          'Operating System :: POSIX',
          'Operating System :: Unix',
          'Operating System :: MacOS',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      ],
      )
