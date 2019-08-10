
Installing brutifus
===================

brutifus is available on pypi, which makes its installation easier than ever. 
In a terminal, type:
::

   pip install brutifus

And that should take care of things.

The most recent release of brutifus is also available for download from its 
`Github repository <https://github.com/fpavogt/fcmaker/releases/latest/>`_. 
Interested users can fork the brutifus repository if they want to get access to the 
latest updates not yet released. *Push requests* for bug fixes and new features are 
welcome and will be examined in detail. 
      
Requirements
------------
brutifus is written in Python 3. The following packages are required for it to work 
properly:

* astropy (3.0 or above)
* astroquery (0.3.4 or above)
* matplotlib (3.0.0 or above)
* numpy (1.14.2 or or above)
* photutils (0.6.0 or above)
* scipy (1.1.0 or above)
* statsmodel (0.9.0 or above)

Testing the installation
------------------------

In a terminal shell, try calling the high-level brutifus help::
 
   brutifus --help
 
If that works, chances are, you will probably be fine ...

.. todo:: Implement an actual set of tests ... sigh!
