.. brutifus documentation master file

.. include:: ./substitutions.rst

brutifus |release| |watch| |stars|
==================================

|license| |github| |pypi| |DOI_latest| |ascl| |astropy| |last-commit|

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

brutifus is modular: `execution` files allow the user to specify which recipes are to be
executed, and in which order. Wherever feasible/beneficial, brutifus routines can use up
to ``ny`` cpus at once (where ``nx*ny`` is the number of spaxels in the datacube;
``ny~300`` for a MUSE cube) to speed up the processing of large datasets.

.. note:: brutifus is very much a work in progress. Modules and recipes will get added as
          the need arises. The code is released publicly because the devs believe in open
          Science.

          You're welcome to try it out ... so long as you remember that the code comes with
          no warranty!


.. toctree::
   :caption: Table of contents
   :maxdepth: 1

   Home <self>
   installation
   running
   troubleshooting
   faq
   acknowledge
   license
   changelog
   Contributing <https://github.com/brutifus/brutifus/blob/develop/CONTRIBUTING.md>
   Github <https://github.com/fpavogt/brutifus>
   modules
