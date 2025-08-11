Usage
=====

Once installed, import and initialize psico in PyMOL to expose its commands:

.. code-block:: python

   import psico.fullinit

You can then use :code:`help psico` on the PyMOL command line to discover
available commands, or search via the ``apropos`` helper:

.. code-block:: text

   PyMOL> apropos align

To generate a standalone HTML reference of all psico commands from within
PyMOL:

.. code-block:: python

   import psico.helping as H
   H.write_html_ref('psico-commands.html')

Many features require optional Python packages or external command line tools
(see README for details). Those are not mandatory for installation but enable
additional functionality when present.
