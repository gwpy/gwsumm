#############################
Configuring a simple HTML tab
#############################

.. currentmodule:: gwsumm.tabs

A :class:`Tab` is single HTML web-page containing some data. It can be as
simple as containing some text, or can be told to generated a number of plots
from GW interferometer data and display them in a specific format.

For full technical details on the `Tab` classes available, please `read the
tabs API page <../tabs>`_.

====
Role
====

A tab can take on one of two roles, depending on its configuration:

======  ======================================================================
Parent  The summary-page for a set of subordinate child tabs, displayed as the
        heading for a dropdown menu in the HTML navigation bar of the output
Child   A subordinate child tab of a given parent, linked under the relevant
        dropdown menu in the HTML navigation bar of the output
======  ======================================================================

=============
Configuration
=============

Each tab class provides a :meth:`Tab.from_ini` `classmethod`, allowing users to
create a new tab from a ``[section]`` in a configuration file.
Every tab should be configured with the following options:

===========  ===========================================================
`~Tab.type`  type of tab to configure [optional, default: ``'default'``]
`~Tab.name`  name of this tab
===========  ===========================================================

The ``type`` option, whilst optional, is recommended, mainly to make the
configuration more transparent to other users, who might not know which tab type
is the default.

Additionally, all tabs can be configured with the following keys:

.. autosummary::
   :nosignatures:

   ~Tab.shortname
   ~Tab.parent
   ~Tab.group

See the detailed definitions of each attribute for defaults.

-------------
`ExternalTab`
-------------

The `ExternalTab` allows users to embed any HTML page from the same domain
into a GWSumm page.
An `ExternalTab` can be configured as follows:

.. code-block:: ini

   [tab-external]
   type = external
   name = 'My results'
   url = '/~duncan.macleod/analysis/results/summary.html'

.. note::

   Only URLs on the same domain can be included by default on most servers.
   This is not a restriction of GWSumm, rather a safety measure of the Apache,
   and other, web server protocols.

---------
`PlotTab`
---------

The `PlotTab` allows users to embed any number of images into a new tab,
and choose the layout.
With this type of tab, users specify images to include by giving options of the
form ``X = /url/to/plot``, with ``X`` an integer, increasing for each plot,
and ``/url/to/plot`` the web URL at which the plot can be found.
Also, users can give the following extra options:

.. autosummary::
   :nosignatures:

   ~PlotTab.foreword
   ~PlotTab.afterword

An example `PlotTab` could be configured as follows:

.. literalinclude:: ../../share/examples/matplotlib.ini
   :language: ini

Here we have imported three examples plots from the
`matplotlib <http://matplotlib.org>`_ examples with a :attr:`~PlotTab.layout`
of one plot on the top row (full size), and two plots on the second row.

---------
`DataTab`
---------

The workhorse of the GWSumm package, at least as used by the LIGO Scientific
Collaboration is the `DataTab`. Please read this page for details on
configuring one of these

.. toctree::

   data
