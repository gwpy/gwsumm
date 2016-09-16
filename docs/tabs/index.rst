.. _tabs:

.. currentmodule:: gwsumm.tabs

####################
Introduction to Tabs
####################

A `Tab` is a single, configurable page of output, containing some data.
Each `Tab` is written in its own HTML page, and can be written to contain
any set of data, with any format.

The basic object provided by :mod:`gwsumm.tabs` is the `Tab`, which allows
embedding of arbitrary HTML text into a standard surrounding HTML framework.
The `Tab` also defines the API for other tabs.

----------------
Simple `Tab` use
----------------

A simple `Tab` can be created in only two steps

.. code-block:: python

   mytab = Tab('My first tab')
   mytab.write_html("This tab doesn't do very much")

This will create two directories under the current one,

- ``my_first_tab/`` containing the HTML for the new tab
- ``static/`` containing all of the stylesheets and javascript files that make the page functional

The content given to `Tab.write_html` is passed along untouched, so can contain any HTML you like.

-------------------
Generating websites
-------------------

The :ref:`next page <websites>` will guide you through created groups of tabs
and combining them to generate a fully-fledged website complete with
navigation.
