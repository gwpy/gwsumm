================================================
Gravitational-wave Summary Information Generator
================================================

GWSumm is a python toolbox used by the LIGO Scientific Collaboration to summarise and archive sundry facets of the performance of the LIGO instruments, and archive these data in a nested HTML structure.

|PyPI version| |Conda version|

|DOI| |License| |Supported Python versions|

|Build Status| |Coverage Status| |Code Climate|

https://gwsumm.readthedocs.io

------------
Installation
------------

GWSumm is best installed with `conda`_:

.. code:: bash

   conda install -c conda-forge gwsumm

but can also be installed with `pip`_:

.. code:: bash

   python -m pip install gwsumm

------------
Contributing
------------

All code should follow the Python Style Guide outlined in `PEP 0008`_;
users can use the `flake8`_ package to check their code for style issues
before submitting.

See `the contributions guide`_ for the recommended procedure for
proposing additions/changes.

.. _PEP 0008: https://www.python.org/dev/peps/pep-0008/
.. _flake8: http://flake8.pycqa.org
.. _the contributions guide: https://github.com/gwpy/gwsumm/blob/master/CONTRIBUTING.md
.. _conda: https://conda.io
.. _pip: https://pip.pypa.io/en/stable/

.. |PyPI version| image:: https://badge.fury.io/py/gwsumm.svg
   :target: http://badge.fury.io/py/gwsumm
.. |Conda version| image:: https://img.shields.io/conda/vn/conda-forge/gwsumm.svg
   :target: https://anaconda.org/conda-forge/gwsumm/
.. |DOI| image:: https://zenodo.org/badge/DOI/10.5281/zenodo.2647609.svg
   :target: https://doi.org/10.5281/zenodo.2647609
.. |License| image:: https://img.shields.io/pypi/l/gwsumm.svg
   :target: https://choosealicense.com/licenses/gpl-3.0/
.. |Supported Python versions| image:: https://img.shields.io/pypi/pyversions/gwsumm.svg
   :target: https://pypi.org/project/gwsumm/
.. |Build Status| image:: https://github.com/gwpy/gwsumm/actions/workflows/build.yml/badge.svg?branch=master
   :target: https://github.com/gwpy/gwsumm/actions/workflows/build.yml
.. |Coverage Status| image:: https://codecov.io/gh/gwpy/gwsumm/branch/master/graph/badge.svg
   :target: https://codecov.io/gh/gwpy/gwsumm
.. |Code Climate| image:: https://codeclimate.com/github/gwpy/gwsumm/badges/gpa.svg
   :target: https://codeclimate.com/github/gwpy/gwsumm
