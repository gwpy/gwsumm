.. _cli-page:

######################
Command-line interface
######################

GW Summary
==========

The primary interface for GWSumm is a command-line utility called `gw_summary`.
For a full explanation of the available command-line arguments and options, you
can run

.. command-output:: python ../bin/gw_summary --help

This tool can be run in four modes: daily, weekly, and monthly analyses, and
a specific range of GPS times.

Day mode
--------

To run in daily summary mode, the following command-line options are available:

.. command-output:: python gw_summary day --help

Week mode
---------

The arguments in weekly mode are as follows:

.. command-output:: python gw_summary week --help

Month mode
----------

In monthly mode:

.. command-output:: python gw_summary month --help

GPS mode
--------

To run within a specific (but arbitrary) range of GPS seconds:

.. command-output:: python gw_summary gps --help

Batch mode
==========

To stage a batch of analyses with a large collection of configuration files,
as is done in embarrassingly parallel fashion when the summary pages run
online, you can use the `gw_summary_pipe` command-line utility. This tool
uses `HT Condor <https://research.cs.wisc.edu/htcondor/>`_ to schedule
and run jobs in parallel.

To see all the available arguments and options for this tool, you can run
with `--help` as usual:

.. command-output:: python gw_summary_pipe --help
