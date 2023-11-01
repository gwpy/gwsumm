.. _modes:

.. currentmodule:: gwsumm.tabs

###########
`Tab` modes
###########

In its simplest form, the `Tab` is essentially a blank canvas on which to write whatever you want.
However, the original mandate for GWSumm was to provide a framework in which to generate automatic summaries of LIGO data, over a given timescale.

To handle data processing, rather than static HTML generation, each `Tab` has a type, based on its relation to any interval in time

.. currentmodule:: gwsumm.tabs.core

.. autosummary::
   :nosignatures:
   :toctree: ../api

   StaticTab
   IntervalTab
   EventTab

The type of a `Tab` is set automatically when it is created based on the value of the :attr:`~Tab.mode` attribute, so you don't need to remember the above objects.

=====
Modes
=====

GWSumm currently support seven different `Tab` modes:

==========  ====  =============================================================
Mode        Enum  Description
==========  ====  =============================================================
``STATIC``  0     No associated time interval
``EVENT``   1     Associated with a single GPS time, normally around an event
``GPS``     2     Simple (arbitrary) GPS ``[start, end)`` interval
``DAY``     10    One UTC 24-hour day
``WEEK``    11    One 7-day week
``MONTH``   12    One calendar month
``YEAR``    13    One calendar year
==========  ====  =============================================================

===============
Assigning modes
===============

Each `Tab` will be assigned a mode (unless specified as follows, the default mode
is ``STATIC``). The assignment can be done on a per-tab basis by
passing the `~Tab.mode` keyword argument when creating a `Tab`, or
globally, by using the :meth:`gwsumm.mode.set_mode`. The latter sets
the default mode for all subsequent tabs created in this session.

If a :attr:`~Tab.mode` is given that associates with a GPS time or
times, these must be given via the `~IntervalTab.span` or
`~EventTab.gpstime` keyword arguments, otherwise a `TypeError` will be
raised. The `span` tuple is the ``(GPS start time, GPS end time)``

.. code-block:: python

   >>> tab = Tab('My first tab', mode='day', span=(0, 100))
   >>> print(tab.mode, tab.span)
   (10, Segment(0, 100))
   >>> tab = Tab('My first tab', mode='EVENT', gpstime=101)
   >>> print(tab.mode, tab.gpstime)
   (1, LIGOTimeGPS(101,0))
