.. brutifus documentation master file, created by
   sphinx-quickstart on Fri Apr 22 11:23:38 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

brutifus |release|
==================

.. warning::
   These pages describe brutifus |release|. The module is currently under construction, and 
   NOT yet ready for deployment. For more information, send an email to 
   frederic.vogt@alumni.anu.edu.au.


brutifus is a set of Python routines designed to process datacubes from Integral Field 
Spectrographs, and in particular MUSE on the VLT/UT-4 at Paranal Observatory. **The focus 
of brutifus is primarily set on performing a detailed emission line analysis of supernova
remnants.**

Some of the features of brutifus include:

  - **the direct use of reduced datacubes** (e.g. produced from the official MUSE data reduction 
    pipeline) without requiring prior modifications,
  - **the fitting of the stellar/nebular continuum** using either a non-parametric Locally 
    Weighted Scatterplot Smoothing (LOWESS) technique,
  - **the fitting of emission lines** via ``mpfit`` that uses the Levenberg-Marquardt technique 
    to solve the least-squares problem, 
  - the ability to use **a varying number of gaussian components** for different emission 
    lines, tie parameters to one another, and set upper and lower bounds,
  - **an automated routine for identifying structures in the data** and 
    define associated extraction apertures, with the ability to then refine the selection 
    interactively,
  - **a modular structure** (inspired by 
    `pywifes <http://pywifes.github.io/pipeline/>`_; `Childress+, 2014 
    <http://adsabs.harvard.edu/abs/2014Ap%26SS.349..617C>`_) allowing users to choose 
    specific processing steps, or add their own with minimal efforts, and


brutifus can use up to ``ny`` cpus at once (where ``nx*ny`` is the number of spaxels in the
datacube; ``ny~300`` for one MUSE cube) to speed up the processing of large datasets. 

.. note::

    You can track the latest changes in the code via the `associated Github repository 
    <https://github.com/fpavogt/brutus>`_.
    
    See also the :ref:`changelog` for a global overview.
    

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

Indices and tables
--------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`


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
 
