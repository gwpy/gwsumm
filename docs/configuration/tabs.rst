#######################
Configuring an HTML tab
#######################

.. currentmodule:: gwsumm.tabs.core

A :class:`tab <SummaryTab>` is single HTML web-page containing some data.
It is defined in the `gwsumm` package as follows:

.. autosummary::
   :toctree: ../_generated

   gwsumm.tabs.core.SummaryTab

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

Each tab is configured as a ``[section]`` in a configuration file.
Each tab should be configured with the following required set of keys:

.. autosummary::
   :nosignatures:
   :toctree: ../_generated

   ~SummaryTab.name

This configuration can be extend by defining any of the following optional
keys:

.. autosummary::
   :nosignatures:
   :toctree: ../_generated

   ~SummaryTab.parent
   ~SummaryTab.states
   ~SummaryTab.layout

See the detailed definitions of each attribute for defaults.

Crucially, each tab should be configured to contain any number of plots.
These should be given as ``X = plot definition`` pairs, where X is an integer
defining the ordering of each plot.
The ``plot definition`` can be given in a small number of ways:

Basic plots
-----------

The simplest plot definition is just a comma-separated list of data channels,
and a plot type indicator, e.g:

.. code-block:: ini

   [wind speed]
   1 = H1:PEM-LVEA_WIND,H1:PEM-EX_WIND,H1:PEM-EY_WIND timeseries

This configuration style is useful for simple plots that do not require any
high degree of customisation, and won't be repeated in any other tabs.

Channel groups
--------------

If a given plot contains a high number of channels, meaning the definition
line might become hard to read, or the list of channels will be used in a
number of other plots, meaning you don't want to repeat it over and over
again, the channels can be defined in a section and the section name given:

.. code-block:: ini

   [channels-wind]
   channels = H1:PEM-LVEA_WIND,H1:PEM-EX_WIND,H1:PEM-EY_WIND

   [wind speed]
   1 = channels-wind timeseries

This example is a bit more readable than that for the basic plot above, but
will give exactly the same output.

Plot sections
-------------

Similar to the channel groups above, if a particular plot style is likely to
be repeated in a number of tabs, or contains a number of very specific
configuration choices, it can be defined in its own section:

.. code-block:: ini

   [plot-wind]
   ; plot configuration goes here

   [wind speed]
   1 = plot-wind

where the ``[plot-wind]`` section contains the channel definition, and all
other necessary information.

Please see `plot configuration <plot>`_ for details on configuring the plots.

.. note::

   If the plot section does not contain a channel definition, for example if
   it will be repeated in different tabs with different channel sets, the
   channel lists or groups can be chained with a plot section:

   .. code-block:: ini

      [wind speed]
      1 = channels-wind plot-wind

   If both the channel section and plot sections contain channel lists, the
   channel group or list takes precedence.

Plot customisation
------------------

Any plot as defined in any of the above ways can be customised within the tab
section by giving any other plotting keys prepended with the plot number and
a dash, e.g:

.. code-block:: ini

   [wind speed]
   1 = channels-wind plot-wind
   1-ylim = 0,100
   1-title = "Wind speed [m/s]"

All customisations take precendence over values defined in plot sections, and
are parsed in the same way.

.. note::

   The tab configuration parser will attempt to evaluate customisation values
   into the relevant python types, e.g. above the ``1-ylim`` value would be
   parsed into a 2-tuple as ``(0, 100)``. However, if evaluation fails, the
   raw string value will be stored, meaning in the ``1-title`` key above, the
   quotation marks are not strictly required, but are a good idea.

A full example
--------------

The following code block, taken from the actual LIGO configuaration for the
Livingston Observatory seismic isolation for Basic Symmetric Chamber 1:

.. code-block:: ini

   [tab SEI BSC1]
   name = BSC1
   parent = Seismic
   states = All,BSC1-ODC
   layout = 1,2,1,2
   ; Platform motion spectrum (with ground reference)
   1 = %(ifo)s:HPI-BS_SENSCOR_Z_FIR_IN1_DQ,
       %(ifo)s:HPI-ITMY_BLND_L4C_Z_IN1_DQ,
       %(ifo)s:ISI-ITMY_ST1_BLND_Z_T240_CUR_IN1_DQ,
       %(ifo)s:SUS-ITMY_M0_ISIWIT_L_DQ plot-seismic-spectrum
   1-labels = 'Ground (Z)','HEPI (Z)','ISI Stage 1 (Z)','Platform motion (L)'
   1-ylim = 1e-4,1e5
   ; Platform motion spectrogram
   2 = %(ifo)s:SUS-ITMY_M0_ISIWIT_L_DQ plot-ground-spectrogram
   2-clim = 1e-4,1e5
   2-title = 'ITMY platform motion'3 = %(ifo)s:SUS-ITMY_M0_ISIWIT_L_DQ plot-ground-median-spectrogram3-title = %(2-title)s
   ; ODC
   4 = %(ifo)s:HPI-ITMY_ODC_CHANNEL_LATCH,%(ifo)s:ISI-ITMY_ODC_CHANNEL_LATCH statevector
   4-title = 'BSC1 ODCs (HPI/ISI)'

Here the 'BSC1' tab, under the 'Seismic' dropdown menu is defined to have a
small number of plots in a given layout (1 plot at the top, then 2 underneath,
then another 1, then 2 on all subsequent lines