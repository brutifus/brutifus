.. _changelog:

Changelog
=========

TODO:
   - `run_crude_snr_maps`: compute proper S/N for bright stars (not <>/std !)
   - add scalebar to plots
   - enable minorticks to plots (see astropy bug [8126](https://github.com/astropy/astropy/issues/8126))
   - check that white ticks appear in RGB images with imshow
   - allow `run_plot_RGB` to plot RGB or single plane files (easy now without aplpy)
   - when drawing contours in `run_skysub()`, show the actual pixel edges
   - allow an automated selection of the sky areas -> requires a better snr estimate ?
   - allow to abort the continuum fitting with CTRL-C

v2018.11.1: F.P.A. Vogt
 - added the 'raw_cube' name to fn_list for easy access

v2018.11.0: F.P.A. Vogt
 - ported brutus to python 3, renamed it to brutifus (available pip name)
 - removed ppxf from brutifus (focus on SNRs for now)
 - removed pyqz from brutifus (outdated)
 - added sky-subtraction step
 - initial upload to Github, initial upload to pip
