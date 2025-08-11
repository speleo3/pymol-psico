Installation
============

psico supports Python 3.7+ and PyMOL 2.3+.

Install from source (recommended for latest):

.. code-block:: bash

   pip install .

Install via a conda-compatible solver:

.. code-block:: bash

   micromamba install speleo3::pymol-psico

Enable all psico commands in PyMOL by adding this to ``~/.pymolrc.py`` or by
entering it into the PyMOL command line:

.. code-block:: python

   import psico.fullinit
