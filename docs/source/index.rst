.. brutifus documentation master file

.. |DOI_latest| image:: https://zenodo.org/badge/157203434.svg
   :target: https://zenodo.org/badge/latestdoi/157203434

.. |ascl| image:: https://img.shields.io/badge/ascl-1903.004-blue.svg?colorB=262255
   :target: http://ascl.net/1903.004 

.. |pypi| image:: https://img.shields.io/pypi/v/brutifus.svg?colorB=<brightgreen>
    :target: https://pypi.python.org/pypi/brutifus/
    
.. |last-commit| image:: https://img.shields.io/github/last-commit/fpavogt/brutifus.svg?colorB=e6c000
   :target: https://github.com/fpavogt/brutifus

.. |issues| image:: https://img.shields.io/github/issues/fpavogt/brutifus.svg?colorB=b4001e   
   :target: https://github.com/fpavogt/brutifus/issues

.. |astropy| image:: http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat
    :target: http://www.astropy.org/
   
.. |stars| image:: https://img.shields.io/github/stars/fpavogt/brutifus.svg?style=social&label=Stars
   :target: https://github.com/fpavogt/brutifus/

.. |watch| image:: https://img.shields.io/github/watchers/fpavogt/brutifus.svg?style=social&label=Watch
   :target: https://github.com/fpavogt/brutifus/
   
   
.. |github| image:: https://img.shields.io/github/release/fpavogt/brutifus.svg
   :target: https://github.com/fpavogt/brutifus/releases   


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
          
Contents
------------------
.. toctree::
   :maxdepth: 1

   Home <self>
   gallery
   installation
   running
   faq
   changelog
   acknowledge
   modules/modules
   Github <https://github.com/fpavogt/brutifus>


----

Copyright notice
*******************
 
This file is part of the brutifus Python module.
The brutifus Python module is free software: you can redistribute it and/or modify it under 
the terms of the GNU General Public License as published by the Free Software Foundation, 
version 3 of the License.
 
The brutifus Python module is distributed in the hope that it will be useful, but without any 
warranty; without even the implied warranty of merchantability or fitness for a particular 
purpose.  See the GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License along with the brutifus 
Python module.  If not, see http://www.gnu.org/licenses/ .
 
