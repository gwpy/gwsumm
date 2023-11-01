#########################################################
GWSumm: the gravitational-wave summary information system
#########################################################

|PyPI version| |Conda version|

|DOI| |License| |Supported Python versions|

The `gwsumm` package is a tool used by the
`Laser Interferometer Gravitational-wave Observatory (LIGO) <http://www.ligo.org>`_
to collect, aggregate, and visually summarise the sundry data produced
throughout the experiment in order to characterise instrument performance.

The output of this package, known internally as the 'summary pages', give an
archive of a number of figures or merit, including time-series amplitude
trends, frequency spectra and spectrograms, and transient event triggers.

This package includes a collection of command-line utilities and a python
module:

.. code:: python

   import gwsumm

Installation
============

GWSumm is best installed with `conda`_:

.. code:: bash

   conda install -c conda-forge gwsumm

but can also be installed with `pip`_:

.. code:: bash

   python -m pip install gwsumm

Note, users with `LIGO.ORG` credentials have access to a software
container with a regularly-updated build of GWSumm. For more
information please refer to the
`LSCSoft Conda <https://docs.ligo.org/lscsoft/conda/>`_ documentation.

Contributing
============

All code should follow the Python Style Guide outlined in `PEP 0008`_;
users can use the `flake8`_ package to check their code for style issues
before submitting.

See `the contributions guide`_ for the recommended procedure for
proposing additions/changes.

The GWSumm project is hosted on GitHub:

* Issue tickets: https://github.com/gwpy/gwsumm/issues
* Source code: https://github.com/gwpy/gwsumm


License
-------

GWSumm is distributed under the `GNU General Public License`_.

.. toctree::
   :maxdepth: 1
   :hidden:

   overview
   cli
   automation
   tabs/index
   tabs/websites
   tabs/modes
   tabs/api
   states
   plots
   modes
   api/index

.. _PEP 0008: https://www.python.org/dev/peps/pep-0008/
.. _flake8: http://flake8.pycqa.org
.. _the contributions guide: https://github.com/gwpy/gwsumm/blob/master/CONTRIBUTING.md
.. _conda: https://conda.io
.. _pip: https://pip.pypa.io/en/stable/
.. _GNU General Public License: https://github.com/gwpy/gwsumm/blob/master/LICENSE

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
