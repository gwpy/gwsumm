.. currentmodule:: gwsumm.tabs

#######################
Configuring a `DataTab`
#######################

While the `Tab`, `ExternalTab`, and `PlotTab` classes provide methods to
build an HTML page from some pre-defined content, the `DataTab` works in a
different model, allowing users to configure content that will be generated
on-the-fly.
The `DataTab` is configured by defining any number of :doc:`plots </plots>`,
with the raw data being read in and processed at runtime.

=================
Basic information
=================

Just a reminder that each `DataTab` should be configured with the following
basic options:

.. autosummary::

   ~DataTab.name
   ~DataTab.shortname
   ~DataTab.group
   ~DataTab.parent
   ~DataTab.layout

.. note::

   Of the above, only the ``name`` is required, all others are optional.

============================================
Configuring a :class:`~gwsumm.plot.DataPlot`
============================================

Plots in a `DataTab` are configured from an INI file in a similar manner to
those in a `PlotTab`, namely by giving options of the form

.. code-block:: ini

   X = <channels> <plot-type>

As for the `PlotTab`, ``X`` can be any integer, and is used to position the
plots in order in the HTML output (by default :mod:`ConfigParser` doesn't
preserve option ordering when INI files are read).

The ``<channels>`` specification can either be a comma-separate (no spaces)
list of
`data channels <http://gwpy.github.io/docs/stable/detector/channel.html>`_, or
the name of an INI ``[section]`` that defines the list of channels to be used.

Similarly, the ``<plot-type>`` specification can either be the
:attr:`type <gwsumm.plot.SummaryPlot.type>` attribute of the relevant plot
class, or the name of an INI ``[section]`` that defines a customised plot type.

For example:

.. code-block:: ini

   [tab-spectra]
   name = Spectra
   1 = L1:IMC-F_OUT_DQ spectrum
   2 = hoft-channels hoft-spectrum

would define a :class:`~gwsumm.plot.SpectrumDataPlot` of the
``L1:IMC-F_OUT_DQ`` channel (with no customisation), and a customised
:class:`~gwsumm.plot.SpectrumDataPlot` of the list of channels specified in
the ``[hoft-channels]`` section of the same INI file.
More details on configuring customised plots and channel sets are given below.

============================================
:class:`~gwsumm.plot.DataPlot` customisation
============================================

Plots for a `DataTab` can be customised in one of two ways: inline, or in a
``[section]``.

---------------------
Section configuration
---------------------

Plots can also be customised in a dedicated ``[section]`` of the INI file,
which defines at least the ``type`` of the plot.
Any other options given will be handled by the :meth:`process` of the specific
`~gwsumm.plot.DataPlot` class, the majority of which will be used to customise
the :class:`Axes` of the plot by executing the relevant :meth:`set_xxx` method.
For example:

.. code-block:: ini

   [hoft-spectrum]
   type = spectrum
   xscale = 'log'
   yscale = 'log'
   xlim = 10,4000
   ylim = 1e-24,1e-20
   xlabel = 'Frequency [Hz]'

This will define a custom plot type called ``hoft-spectrum``, a sub-type of
the built-in ``spectrum``, with logarithmic x- and y-axes, specific
x- and y-axes limits, and a fixes x-axis label.

Using a ``[section]`` to customise a specific plot is particularly useful when
the set of customisations will be used multiple times, as in this example:

.. code-block:: ini

   [tab-spectra]
   name = Spectra
   1 = L1:LDAS-STRAIN hoft-spectrum
   2 = H1:LDAS-STRAIN hoft-spectrum

   [hoft-spectrum]
   type = spectrum
   xscale = 'log'
   yscale = 'log'
   xlim = 10,4000
   ylim = 1e-24,1e-20
   xlabel = 'Frequency [Hz]'

Here the ``[hoft-spectrum]`` customisation has been applied to two different
plots in different tabs, since the ``h(t)`` gravitational-wave strain spectrum
for each of the LIGO observatories (LHO and LLO) should be very similar.

--------------------
Inline customisation
--------------------

In cases where plot customisations will only be used once, they can be given
as options of the form ``X-<key> = <value>`` where ``X`` is the plot to be
customised, and ``<key>`` and ``<value>`` are the option and value to be
passed to the plot processor.
For example:

.. code-block:: ini

   [tab-spectra]
   name = spectra
   1 = L1:IMC-F_OUT_DQ spectrum
   1-ylabel = 'Frequency noise [Hz]'

This will configure a :class:`~gwsumm.plot.SpectrumDataPlot` of the
``L1:IMC-F_OUT_DQ`` channel as previously, but will set a different y-axis
label than before.

As for options given in plot sections, the ``<key>`` can be given as anything
for which there is a :meth:`set_\<key\>` method in the underlying :class:`Axes`.

The two methods of customisation can be combined, as follows:

.. code-block:: ini

   [tab-spectra]
   name = Spectra
   1 = L1:LDAS-STRAIN hoft-spectrum
   1-title = 'LIGO Livingston Observatory'
   1-ylim = 5e-25,1e-20
   2 = H1:LDAS-STRAIN hoft-spectrum
   2-title = 'LIGO Hanford Observatory'

   [hoft-spectrum]
   type = spectrum
   xscale = 'log'
   yscale = 'log'
   xlim = 10,4000
   ylim = 1e-24,1e-20
   xlabel = 'Frequency [Hz]'

.. warning::

   Any inline customisations will override those given in a ``[section]``.
   In this example, the ``1-ylim`` option in ``[tab-spectra]`` will supercede
   the ``ylim`` option in ``[tab-spectrum]`` for plot ``1``, but will have no
   impact on plot ``2``.

.. currentmodule:: gwsumm.state

============================
Configuring a `SummaryState`
============================

Each `~gwsumm.tabs.DataTab` is configured to run over any number of
:doc:`states <../states>`.
By default a `~gwsumm.tabs.DataTab` is configured to run over the special 'all' state,
simply including all time from start to finish.
Custom states can be configured in their own ``[section]`` by providing the
following options

 .. autosummary::

   ~SummaryState.name
   ~SummaryState.definition
   ~SummaryState.key

For example:

.. code-block:: ini

   [state-science]
   name = 'Science'
   definition = L1:DMT-SCIENCE:1

Here the 'Science' state is defined by the single data-quality flag
``L1:DMT-SCIENCE:1``.
This state would then be used in the above example as follows:

.. code-block:: ini

   [tab-spectra]
   name = Spectra
   states = Science
   1 = L1:LDAS-STRAIN hoft-spectrum
   1-title = 'LIGO Livingston Observatory'
   1-ylim = 5e-25,1e-20
   2 = H1:LDAS-STRAIN hoft-spectrum
   2-title = 'LIGO Hanford Observatory'

The ``states`` option is always given in the plural form, and in general
can be a comma-separated list of states.
To process data for both the 'Science' state and for all times, we could give

.. code-block:: ini

   states = Science,All

Here plots and data summaries are generated for both states, with the first in
the list ('Science') being set as the default.

.. note::

   The ``key`` option can be given to assign a unique key to a `SummaryState`,
   useful in the event that you want multiple states to have the same name
   for HTML output, but with different definitions.
   For example:

   .. code-block:: ini

      [state-sei-odc-itmx]
      name = ODC
      key = SEI-ODC-ITMX
      definition = L1:ODC-ISI_ITMX_SUMMARY:1&L1:ODC-HPI_ITMX_SUMMARY:1

      [state-sei-odc-itmy]
      name = ODC
      key = SEI-ODC-ITMY
      definition = L1:ODC-ISI_ITMY_SUMMARY:1&L1:ODC-HPI_ITMY_SUMMARY:1

   These sections define the 'ODC' states for both of the input test masses
   (ITMs, X-arm and Y-arm) with the same name but with unique keys.

.. currentmodule:: gwsumm.tabs

=========================
Configuring data channels
=========================

For most `DataTab` output, the main source of data will be channels - streams
of data recorded from the myriad interferometer sensors.
Users can configure channels, or groups of channels, in order to control the
way these data are read in, manipulated, and displayed.

A single channel can be configured in a section named for that channel, for
example:

.. code-block:: ini

   [L1:LDAS-STRAIN]
   frametype = L1_LDAS_C00_L2
   unit = 'strain'
   resample = 4096

This section has configured the ``L1:LDAS-STRAIN`` channel to be read from
a specific frame file type and resampled to a given rate, with a specific unit.

Similary, groups of channels can be configured together in any section whose
name starts with ``channels-`` and that contains the ``channels`` option:

.. code-block:: ini

   [channels-hoft]
   channels = L1:PSL-ISS_PDA_OUT_DQ,
              L1:PSL-ISS_PDB_OUT_DQ
   frametype = L1_R
   resample = 8192
   unit = 'RIN'
   frequency-response = [2*pi*118.21261,2*pi*3.7887063,2*pi*2.792676],[2*pi*71.39292e-3,2*pi*71.535004e-3],-3.1972880804037495e-07

Here we have configured the in-loop and out-of-loop photodiodes for the
Pre-stabilized LASER (PSL) Intensity Stabilisation Servo (ISS) with a frametype
and re-sampling rate, a unit, and a ZPK-format filter (zero-pole-gain) to
transform the raw data into the correct units.

The following channel configuration options are understood for a `DataTab`:

==================  ===========================================================
channels            comma-separated list of channels to configure
frametype           GWF type code to query for data frames. Use a bar-separated
                    part to specify a string to match frame URLS for, e.g.
                    ``frametype = T|RDS9`` will find `T`-type frame files whose
                    full URLs contain the string `RDS9` (looking at you
                    GEO-600)
resample            number of samples per second at which to resample data
unit                physical unit of the data (after any filtering)
filter              time-domain filter to apply. Should be of the form \
                    ``zeros,poles,gain`` where ``zeros`` and ``poles`` should \
                    be a list of frequencies (``pi`` accepted) and ``gain`` \
                    should be a float
frequency-response  frequency-domain filter to apply. Should be used in \
                    favour of ``filter`` if only frequency-domain data are \
                    to be summarised. Takes the same format as ``filter``
fftlength           length of single Fourier transform (in seconds)
overlap             amount of overlap between successive Fourier transforms \
                    (in seconds)
stride              duration of single spectrogram bin (in seconds)
<xxx>-range         ``low,high`` axis range for plotting channel data. \
                    ``<xxx>`` should be ``psd``, ``asd``, ``frequency`` \
                    or any valid column of the source ``LIGO_LW`` table

==================  ===========================================================

======================
Variable interpolation
======================

The Python :mod:`ConfigParser` module supports variable interpolation inside
a given ``[section]``.
In addition to intra-section interpolation, all variables in the ``DEFAULT``
section of the configuration files are interpolated in other sections.
The `gw_summary` command-line executable sets the following variables in the
``DEFAULT`` section for configurations to use:

==============  ==============================================================
ifo             the two-character interferometer prefix, if given on the \
                command-line
user            the current username
gps-start-time  the integer GPS start time of the requested data stretch
gps-end-time    the integer GPS end time of the requested data stretch
yyyy            the four-integer year of the requested data stretch
yyyymm          the six-integer year and month of the requested data stretch
yyyydd          the eight-integer year, month, and day of the requested data \
                stretch
==============  ==============================================================

The full example below uses the ``%(ifo)s`` interpolation to make the
configuration file independent of interferometer.
Only when this file is passed to `gw_summary` along with the ``--ifo=L1``
argument does it actually become specific to the LIGO Livingston Observatory.

============
Full example
============

The following configuration file is one used by LIGO to monitor the Length
Sensing and Control (LSC) subsystem of the LIGO Livingston Observatory
interferometer.

**Accessed June 26th 2014**

.. literalinclude:: ../../share/examples/lsc.ini
   :language: ini
