.. brutifus documentation master file

.. include:: ./substitutions.rst

brutifus |release|
==================

|pypi| |DOI_latest| |ascl| |astropy| |last-commit|

brutifus is a set of Python routines designed to post-process datacubes from Integral Field
Spectrographs, and in particular MUSE and WiFeS. The spirit of brutifus is to deal with
*generic* post-processing tasks with as little interactions as possible.

So far, brutifus allows to:

  - **correct a cube WCS solution** using Gaia,
  - **generate B&W and RGB images from the cube**,
  - **correct for the galactic reddening along the line of sight**,
  - **fit the stellar/nebular continuum** using the non-parametric Locally Weighted
    Scatterplot Smoothing (LOWESS) technique,
  - **perform an ad-hoc sky-subtraction** if required.

brutifus is modular: `execution` files allow the user to specify which recipe is to be
executed, and in which order. Wherever feasible/beneficial, brutifus routines can use up
to ``ny`` cpus at once (where ``nx*ny`` is the number of spaxels in the datacube;
``ny~300`` for a MUSE cube) to speed up the processing of large datasets.

.. note:: brutifus is very much a work in progress. Modules and recipes will get added as
          the need arises. The code is released publicly because I believe in an open
          Science policy.

          You're welcome to try it out ... so long as you remember that the code comes with
          no warranty!


.. toctree::
   :caption: Table of contents
   :maxdepth: 1

   Home <self>
   gallery
   installation
   running
   troubleshooting
   faq
   acknowledge
   changelog
   Contributing <https://github.com/brutifus/brutifus/blob/develop/CONTRIBUTING.md>
   Github <https://github.com/fpavogt/brutifus>
   modules


----

Copyright notice
****************

This file is part of the brutifus Python module.
The brutifus Python module is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software Foundation,
version 3 of the License.

The brutifus Python module is distributed in the hope that it will be useful, but without any
warranty; without even the implied warranty of merchantability or fitness for a particular
purpose.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the brutifus
Python module.  If not, see http://www.gnu.org/licenses/ .
