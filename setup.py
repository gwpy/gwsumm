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
import glob
import os
from distutils import log

from setuptools import (Command, setup, find_packages)
from setuptools.command.egg_info import egg_info
from setuptools.command.build_py import build_py

import versioneer

if sys.version < '2.6':
    raise ImportError("Python versions older than 2.6 are not supported.")

# set basic metadata
PACKAGENAME = 'gwsumm'
DISTNAME = 'gwsumm'
AUTHOR = 'Duncan Macleod, Alex Urban'
AUTHOR_EMAIL = 'alexander.urban@ligo.org'
LICENSE = 'GPLv3'

cmdclass = {}

# -- versioning ---------------------------------------------------------------

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

# build requirements
setup_requires = [
    'libsass',
    'jsmin',
]

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
    'gwdetchar>=0.5.1',
    'configparser ; python_version < \'3.6\'',
    'enum34 ; python_version < \'3.4\''
]

# testing requirements
if {'pytest', 'test'}.intersection(sys.argv):
    setup_requires.append('pytest_runner')  # python setup.py test
tests_require = [
    'pytest>=2.8'
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

JAVASCRIPT_SOURCES = [
    f for f in glob.glob(os.path.join('share', 'js', '*.js'))
    if not f.endswith('.min.js')]
CSS_SOURCES = glob.glob(os.path.join('share', 'sass', '[!_]*.scss'))


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
          'Programming Language :: Python :: 2.7',
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
