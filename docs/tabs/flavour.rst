.. _flavours:

.. currentmodule:: gwsumm.tabs

##############
`Tab` flavours
##############

In its simplest form, the `Tab` is essentially a blank canvas on which to write whatever you want.
However, the original mandate for GWSumm was to provide a framework in which to generate automatic summaries of LIGO data, over a given timescale.

To handle data processing, rather than static HTML generation, each `Tab` has a flavour, based on its relation to any interval in time

.. currentmodule:: gwsumm.tabs.core

.. autosummary::
   :nosignatures:
   :toctree: ../api

   StaticTab
   IntervalTab
   EventTab

The 'flavour' of a `Tab` is set when it is created, based on the value of the :attr:`~Tab.mode` attribute, for full details see :ref:`modes`.
