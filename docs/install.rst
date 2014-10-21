########################
How do I install GWSumm?
########################

============
Dependencies
============

**Build dependencies**

The GWSumm package has the following build-time dependencies (i.e. required for installation):

* `glue <https://www.lsc-group.phys.uwm.edu/daswg/projects/glue.html>`_
* `NumPy <http://www.numpy.org>`_ >= 1.5
* `Astropy <http://astropy.org>`_ >= 0.3
* `GWpy <https://gwpy.github.io>`_ >= 0.1a5

.. note::

   The `GLUE <https://www.lsc-group.phys.uwm.edu/daswg/projects/glue.html>`_ package isn't available through PyPI, meaning you will have to install it manually from the link.

**Runtime dependencies**

Additionally, in order for much of the code to import and run properly, users are suggested (but not required) to have the following packages:

* `lal <https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html>`_ and `lalframe <https://www.lsc-group.phys.uwm.edu/daswg/projects/lalsuite.html>`_ (same URL)
* `NDS2 <https://www.lsc-group.phys.uwm.edu/daswg/projects/nds-client.html>`_ (including SWIG-wrappings for python)

===================
Installing with pip
===================

GWSumm isn't mature enough to have a formal release, and so the best way to install is to point `pip <https://pip.pypa.io/en/latest/index.html>`_ at the GitHub repository and install from there:

.. code-block:: bash

   pip install --user git+https://github.com/gwpy/gwsumm

The ``--user`` option tells the installer to copy codes into the standard user library path, on linux machines this is

.. code-block:: bash

    ~/.local/lib

while on Mac OS this is

.. code-block:: bash

    ~/Library/Python/X.Y/lib

where ``X.Y`` is the python major and minor version numbers, e.g. ``2.7``.
For either operating system, python will automatically know about these directories, so you don't have to fiddle with any environment variables.

======================
Cloning the repository
======================

The source code for GWSumm is under ``git`` version control, hosted by http://github.com.
You can clone the repository from the Terminal as follows:

.. code-block:: bash

    git clone https://github.com/gwpy/gwsumm.git

You can then, if you wish, install the package by running the ``setup.py`` script as follows:

.. code-block:: bash

    cd gwsumm
    python setup.py install --user

=======================
Available installations
=======================

If you are a member of the LIGO Scientific Collaboration, both the GWSumm and astropy packages are installed on all shared computing centres.

If you use the ``bash`` shell, you can source the following script to set up the environment for the GWSumm package

.. code-block:: bash

    source /home/detchar/opt/gwpysoft/etc/gwpy-user-env.sh

If anyone wants to write an equivalent shell script for the ``csh`` shell, please e-mail it to `Duncan <duncan.macleod@ligo.org>`_.
