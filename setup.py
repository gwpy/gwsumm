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

from __future__ import print_function

import sys
if sys.version < '2.6':
    raise ImportError("Python versions older than 2.6 are not supported.")

import glob
import os

from distutils import log
from setuptools import (Command, setup, find_packages)
from setuptools.command.egg_info import egg_info
from setuptools.command.build_py import build_py

# set basic metadata
PACKAGENAME = 'gwsumm'
DISTNAME = 'gwsumm'
AUTHOR = 'Duncan Macleod'
AUTHOR_EMAIL = 'duncan.macleod@ligo.org'
LICENSE = 'GPLv3'

cmdclass = {}

# -- versioning ---------------------------------------------------------------

import versioneer
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

setup_requires = [
    'pytest-runner',
    'libsass',
    'jsmin',
]
install_requires = [
    'numpy>=1.5',
    'matplotlib>=1.3.0',
    'astropy>=1.0',
    'gwpy>=0.1',
    'trigfind>=0.2.1',
]
requires = [
    'glue',
    'numpy',
    'matplotlib',
    'astropy',
    'gwpy',
    'h5py',
    'lxml',
]
tests_require = [
    'pytest>=2.8'
]
extras_require = {
    'doc': ['sphinx', 'numpydoc', 'sphinx-bootstrap-theme', 'astropy_helpers'],
}

# version-specific packages
if sys.version < '3':
    requires.append('enum34')
if sys.version < '2.7':
    tests_require.append('unittest2')

# -- data files ---------------------------------------------------------------

JAVASCRIPT_SOURCES = [
    f for f in glob.glob(os.path.join('share', 'js', '*.js')) + (
               glob.glob(os.path.join('bootstrap-ligo', 'js', '*.js')))
    if not f.endswith('.min.js')]
CSS_SOURCES = glob.glob(os.path.join('share', 'sass', '[!_]*.scss')) + (
    glob.glob(os.path.join('bootstrap-ligo', 'css', '[!_]*.scss')))


class BuildHtmlFiles(Command):
    """Compile SASS sources into CSS and minify javascript
    """
    description = 'Compile SASS into CSS and minify JS'
    user_options = [
        ('output-style=', None, 'CSS output style'),
    ]

    def initialize_options(self):
        self.output_style = 'compressed'

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

    def minify_js(self):
        jsdir = os.path.join('share', 'js')
        log.info('minifying js under %s' % jsdir)
        for jsfile in JAVASCRIPT_SOURCES:
            filename = os.path.basename(jsfile)
            target = os.path.join(
                self.staticdir, '%s.min.js' % os.path.splitext(filename)[0])
            self.minify(jsfile, target)

    def minify(self, source, target):
        import jsmin
        js = jsmin.jsmin(open(source).read())
        with open(target, 'w') as f:
            f.write(js)
        log.info('minified js written to %s' % target)

    def compile_sass(self):
        sassdir = os.path.join('share', 'sass')
        log.info('compiling SASS under %s to CSS' % sassdir)
        for sassfile in CSS_SOURCES:
            filename = os.path.basename(sassfile)
            target = os.path.join(
                self.staticdir, '%s.min.css' % os.path.splitext(filename)[0])
            self.compile(sassfile, target)

    def compile(self, source, target):
        import sass
        print(source, target)
        css = sass.compile(filename=source, output_style=self.output_style)
        with open(target, 'w') as f:
            f.write(css)
        log.info('%s CSS written to %s' % (self.output_style, target))

    def run(self):
        self.compile_sass()
        self.minify_js()
        if self.staticpackage not in self.distribution.packages:
            self.distribution.packages.append(self.staticpackage)
            log.info("added %s to package list" % self.staticpackage)

cmdclass['build_html_files'] = BuildHtmlFiles

old_build_py = cmdclass.pop('build_py', build_py)


class BuildPyWithHtmlFiles(old_build_py):
    """Custom build_py that compiles SASS+JS sources as well
    """
    def run(self):
        self.run_command('build_html_files')
        old_build_py.run(self)

cmdclass['build_py'] = BuildPyWithHtmlFiles

old_egg_info = cmdclass.pop('egg_info', egg_info)


class EggInfoWithHtmlFiles(old_egg_info):
    """Custom egg_info that compiles SASS+JS sources as well
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

setup(name=DISTNAME,
      provides=[PACKAGENAME],
      version=__version__,
      description=None,
      long_description=None,
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
      requires=requires,
      extras_require=extras_require,
      dependency_links=[
          'http://software.ligo.org/lscsoft/source/glue-1.49.1.tar.gz',
          'http://software.ligo.org/lscsoft/source/dqsegdb-1.2.2.tar.gz',
          'https://github.com/ligovirgo/trigfind/archive/v0.3.tar.gz#egg=trigfind-0.3',
      ],
      data_files=data_files,
      use_2to3=True,
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
