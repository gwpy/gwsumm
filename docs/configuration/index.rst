########################################
How to write your own configuration file
########################################

GWSumm makes use of INI-format configuration files to customise which data
should be read, how it should be manipulated, and how it should be displayed.
Through the `gw_summary` executable, all input from the user, excepting a
small number of command-line options,
are provided through one or more INI-format configuration files.
These files contain sets of ``key = value`` pairs grouped into named
``[sections]`` to define everything from the required data inputs for a given
page of output, to whether or not to display the gravitational-wave amplitude
spectrum with a logarithmic frequency scale.

Each of the following pages will explain how to write an INI configuration
file to setup pretty much any data you could want:

.. toctree::
   :maxdepth: 2

   format
   tabs
   data
