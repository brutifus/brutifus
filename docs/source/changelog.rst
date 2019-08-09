.. _changelog:

Changelog
=========

.. todo::
   - (!) allow to abort the continuum fitting with CTRL-C
   - `run_crude_snr_maps`: compute proper S/N for bright stars (not <>/std !)
   - enable minorticks in colorbars (see astropy bug `[8126] <https://github.com/astropy/astropy/issues/8126>`_)
   - set the color of NaN pixels in plots ... (see astropy bug `[8165] <https://github.com/astropy/astropy/issues/8165>`_)
   - check that white ticks appear in RGB images with imshow
   - when drawing contours in `run_skysub()`, show the actual pixel edges
   - allow an automated selection of the sky areas -> requires a better snr estimate ?
   - (?) add scalebar to plots ... tricky: what if pix scale not uniform, North not up, etc ...

v2019.08.1; F.P.A Vogt
- improved the WCS adjustment step, in case the reference pixel is not in the middle of the image
  (found by J. Suherli with OCam data).
- added automatic inclusion of code version in setup.py
- only load ``Gaia`` from astroquery if we need it, to avoid the need for the internet when it is not needed.

v2019.07.1: F.P.A. Vogt
- added courtesy function to create the execution files, via `python -m brutifus --setup`.
- added a proper entry point, so this can be run as `brutifus --setup`.
- carry over the full original headers when adjusting the WCS values.
- fixed a bug in the WCS adjustment function that led to an error of 0.5 pixel for odd array sizes

v2019.02.3: F.P.A. Vogt
 - following the suggestion from Nick Whyborn, use image convolution to find the best 
   spatial shift required to anchor the cube WCS on Gaia.
 - removed all references to aplpy

v2019.02.2: F.P.A. Vogt
 - added photutils to the list of required packages
 - fixed all docstrings
 - cleaned up code of all old brutus bits not (yet) ported to brutifus. 
 - started working on documentation 
 - first public release!

v2019.02.1: F.P.A. Vogt
 - added run_BWplot step
 - aded WCS correction step
 - worked a bit on refactoring the docstrings

v2019.01.1: F.P.A. Vogt
 - revived the reddening module

v2018.11.1: F.P.A. Vogt
 - added the 'raw_cube' name to fn_list for easy access

v2018.11.0: F.P.A. Vogt
 - ported brutus to python 3, renamed it to brutifus (available pip name)
 - removed ppxf from brutifus (focus on SNRs for now)
 - removed pyqz from brutifus (outdated)
 - added sky-subtraction step
 - initial upload to Github, initial upload to pip
