.. _websites:

.. currentmodule:: gwsumm.tabs

Generating Websites
===================

:ref:`As we have seen <tabs>`, generating standalone pages is trivial using GWSumm. What would be more useful would be to generate linked sets of pages, aka a website.

Navigation
----------

The key difference between standalone pages and a website is the ability to navigate between them. The `Tab.write_html` method will take care of that for you if you pass it all of the tabs:

.. code-block:: python

   from gwsumm.tabs import Tab
   tab1 = Tab('Tab 1')
   tab2 = Tab('Tab 2')
   tabs = [tab1, tab2]
   tab1.write_html('This is tab 1', tabs=tabs)
   tab2.write_html('This is tab 2', tabs=tabs)

This will write each tab into its own directory, as before, but the HTML will now contain an `<nav></nav>` block above the banner to allow navigation between the pages.

Tab parents
-----------

In the above example, each tab is included as a link in the navigation bar. However, in larger websites with many pages, the navigation can quickly become cluttered and will start to overflow the width of the page.
This can be avoided by declaring `~Tab.parent` for sets of tabs:

.. code-block:: python

   tab1 = Tab('Tab 1')
   tab2a = Tab('A', parent='Tab 2')
   tab2b = Tab('B', parent=tab2a.parent)
   tabs = [tab1, tab2a, tab2b]
   tab1.write_html('This is tab 1', tabs=tabs)
   tab2a.write_html('This is tab 2A', tabs=tabs)
   tab2b.write_html('This is tab 2B', tabs=tabs)

Here we have set a `parent` tab for 2A, and used the same for 2B, which creates a dropdown menu in the navigation bar linking to these tabs. 'Tab 2' is never created, but is used only for navigation.

Tab groups
----------

For even larger websites, sets of tabs under a single parent can be
further separated into `groups <Tab.group>`. For example, to put 2A
into group `1` and 2B into group `2`, we can write:

.. code-block:: python

   tab1 = Tab('Tab 1')
   tab2a = Tab('A', parent='Tab 2', group='1')
   tab2b = Tab('B', parent=tab2a.parent, group='2')
   tabs = [tab1, tab2a, tab2b]
   tab1.write_html('This is tab 1', tabs=tabs)
   tab2a.write_html('This is tab 2A', tabs=tabs)
   tab2b.write_html('This is tab 2B', tabs=tabs)
