#####################################
A whistle-stop tour of the INI format
#####################################

The INI file syntax is an intertionally-recognised method of provided basic
configuration options to any program, and is used extensively in, amongst
other projects, the red-hat linux distrobution.
Any `INI` file can be as simple as one line containing a `key` -- the name
of a given variable, perhaps -- and a corresponding `value`:

.. code-block:: ini

    name = John Doe

Most parsers would then return a variable ``name`` containg the string value
``John Doe``, or at least something very similar.

Moving beyond a simple set of variables and their values, the INI format
supports grouping the key-value pairs into named sections.
Any text enclosed in square brackets (``[]``) is interpreted as a section
name:

.. code-block:: ini

   [haggis]
   name = Haggis, Neeps, and Tatties
   ingredients = haggis, neeps, tatties
   time = 2.0

   [marsbar]
   name = Deep-fried Mars bar
   ingredients = Mars bar
   utensils = Deep-fat frier
   time = 0.01

The above example, an extract from the author's recipe book, contains two
sections with unique names, and a matchine set of keys describing the basic
information about each dish.
In this case, a standard file-parser would return a structured object with
two sub-structures, each representing one section of the file.

For more details and a good selections of examples, please see the
documentation of the python :mod:`ConfigParser` module.

.. warning::

   By default, the python :mod:`ConfigParser` module allows separation of
   keys and values using either the colon (``:``) or the equals sign (``=``),
   however, in order to increase functionality and user-friendliness,
   `gwsumm` configuration files must use only the equals sign (``=``) for
   this purpose. Do not use colons to separate key-value pairs.