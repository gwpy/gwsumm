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

import sys
import glob
import os
from distutils import log
from shutil import copyfile

from setuptools import (Command, setup, find_packages)
from setuptools.command.egg_info import egg_info
from setuptools.command.build_py import build_py

# set basic metadata
PACKAGENAME = 'gwsumm'
DISTNAME = 'gwsumm'
AUTHOR = 'Alex Urban, Duncan Macleod'
AUTHOR_EMAIL = 'alexander.urban@ligo.org'
LICENSE = 'GPLv3'

cmdclass = {}

# -- versioning ---------------------------------------------------------------

import versioneer  # noqa: E402
__version__ = versioneer.get_version()
cmdclass.update(versioneer.get_cmdclass())

# -- documentation ------------------------------------------------------------

# import sphinx commands
try:
    from sphinx.setup_command import BuildDoc
except ImportError:
    pass
else:
    cmdclass['build_sphinx'] = BuildDoc

# -- dependencies -------------------------------------------------------------

# install requirements (module-level imported packages)
install_requires = [
    'python-dateutil',
    'lxml',
    'numpy>=1.10',
    'scipy>=1.2.0',
    'matplotlib>=2.2.0',
    'astropy>=1.2.1',
    'lalsuite',
    'lscsoft-glue>=1.60.0',
    'ligo-segments',
    'gwpy>=0.14.2',
    'gwtrigfind',
    'gwdatafind',
    'pygments',
    'MarkupPy',
    'markdown',
    'gwdetchar>=0.5.1',
    'configparser ; python_version < \'3.6\'',
]

# build and test requirements
setup_requires = ['pytest_runner'] if {
    'pytest', 'test'}.intersection(sys.argv) else []
tests_require = [
    'pytest>=2.8,<3.7',
    'pytest-cov',
    'flake8',
]
if sys.version < '3':
    tests_require.append('mock')

# extras
extras_require = {
    'docs': [
        'sphinx',
        'numpydoc',
        'sphinx-bootstrap-theme',
        'astropy_helpers',
    ],
}

# -- data files ---------------------------------------------------------------

SOURCE_FILES = [
     os.path.join('gwbootstrap', 'lib', 'gwbootstrap.min.css'),
     os.path.join('gwbootstrap', 'lib', 'gwbootstrap-extra.min.js'),
]

if not SOURCE_FILES:  # make sure submodule is not empty
    raise ValueError('gwbootstrap submodule is empty, please populate it '
                     'with `git submodule update --init`')


class BuildHtmlFiles(Command):
    """Grab compiled CSS and minified JavaScript
    """
    description = 'Grab compiled CSS and minified JS'

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    @property
    def staticdir(self):
        try:
            return self._static
        except AttributeError:
            self._static = os.path.join(PACKAGENAME, 'html', 'static')
            if not os.path.isdir(self._static):
                os.makedirs(self._static)
                log.info('created static dir in %s' % self._static)
        return self._static

    @property
    def staticpackage(self):
        return self.staticdir.replace(os.path.sep, '.')

    def build_utils(self):
        log.info('copying minified elements')
        for file_ in SOURCE_FILES:
            filename = os.path.basename(file_)
            target = os.path.join(self.staticdir, filename)
            copyfile(file_, target)
            log.info('minified CSS and JS written to %s' % target)

    def run(self):
        self.build_utils()
        if self.staticpackage not in self.distribution.packages:
            self.distribution.packages.append(self.staticpackage)
            log.info("added %s to package list" % self.staticpackage)


cmdclass['build_html_files'] = BuildHtmlFiles

old_build_py = cmdclass.pop('build_py', build_py)


class BuildPyWithHtmlFiles(old_build_py):
    """Custom build_py that grabs compiled CSS+JS sources as well
    """
    def run(self):
        self.run_command('build_html_files')
        old_build_py.run(self)


cmdclass['build_py'] = BuildPyWithHtmlFiles

old_egg_info = cmdclass.pop('egg_info', egg_info)


class EggInfoWithHtmlFiles(old_egg_info):
    """Custom egg_info that grabs compiled CSS+JS sources as well
    """
    def run(self):
        self.run_command('build_html_files')
        old_egg_info.run(self)


cmdclass['egg_info'] = EggInfoWithHtmlFiles

# configuration files
data_files = [
    (os.path.join('etc', PACKAGENAME, 'configuration'),
     glob.glob(os.path.join('share', '*.ini'))),
]

# -- run setup ----------------------------------------------------------------

packagenames = find_packages()
scripts = glob.glob(os.path.join('bin', '*'))

# read description
with open('README.rst', 'rb') as f:
    longdesc = f.read().decode().strip()

setup(name=DISTNAME,
      provides=[PACKAGENAME],
      version=__version__,
      description=("A python toolbox used by the LIGO Scientific "
                   "Collaboration for detector characterisation"),
      long_description=longdesc,
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      url='https://github.com/gwpy/gwsumm',
      packages=packagenames,
      include_package_data=True,
      cmdclass=cmdclass,
      scripts=scripts,
      setup_requires=setup_requires,
      install_requires=install_requires,
      tests_require=tests_require,
      extras_require=extras_require,
      data_files=data_files,
      use_2to3=False,
      classifiers=[
          'Programming Language :: Python',
          'Development Status :: 3 - Alpha',
          'Programming Language :: Python',
          'Programming Language :: Python :: 3.5',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
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
