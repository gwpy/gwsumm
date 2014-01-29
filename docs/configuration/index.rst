########################################
How to write your own configuration file
########################################

All input from the user, excepting a small number of command-line options,
are provided through one or more INI-format configuration files.
These files contain sets of ``key = value`` pairs grouped into named
``[sections]`` to define everything from the required data inputs for a given
page of output, to whether or not to display the gravitational-wave amplitude
spectrum with a logarithmic frequency scale.

.. toctree::

   format

The following pages outline how to define each of the cornerstones of the
output web-pages:

.. autosummary::
   :nosignatures:
   :toctree: ../_generated

   ~gwsumm.tabs.SummaryTab
   ~gwsumm.plot.TabPlot

Simply put, the output web-pages consist of a number of 'tabs' -- individual
pages containing some arbitrary data -- each containing plots -- graphs of
aforementioned arbitrary data.

The tab set can have up to two layers, meaning a tab can either act as the
parent for a number of sub-tabs, or can refer to its parent and itself be a
sub-tab of the parent.

.. toctree::
   :maxdepth: 1

   tabs
   plots
   default
