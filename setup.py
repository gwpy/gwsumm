#!/usr/bin/env python

import sys
import imp
try:
    # This incantation forces distribute to be used (over setuptools) if it is
    # available on the path; otherwise distribute will be downloaded.
    import pkg_resources
    distribute = pkg_resources.get_distribution('distribute')
    if pkg_resources.get_distribution('setuptools') != distribute:
        sys.path.insert(1, distribute.location)
        distribute.activate()
        imp.reload(pkg_resources)
except:  # There are several types of exceptions that can occur here
    from distribute_setup import use_setuptools
    use_setuptools()

import glob
import os
from setuptools import setup, find_packages

PACKAGENAME = 'gwsumm'
DESCRIPTION = 'Gravitational-wave interferometer summary information system'
AUTHOR = 'Duncan Macleod'
AUTHOR_EMAIL = 'duncan.macleod@ligo.org'
LICENSE = 'BSD'
VERSION = '0.0.dev'
RELEASE = 'dev' not in VERSION

from astropy.version_helpers import get_git_devstr, generate_version_py
generate_version_py(PACKAGENAME, VERSION, RELEASE)

setup(name=PACKAGENAME,
      version=VERSION,
      description=DESCRIPTION,
      packages=find_packages(),
      requires=['gwpy'],
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      licence=LICENSE)
