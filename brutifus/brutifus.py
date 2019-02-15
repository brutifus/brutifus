# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018,  F.P.A. Vogt
 
-----------------------------------------------------------------------------------------
 
This file contains the master brutifus routines. Most of these routines call sub-routines,
after setting the scene/loading datasets/etc ...
 
Any processing step MUST have a dediacted routine in this file called 'run_XXX', which can
then refer to any existing/new brutifus/Python module.

Created November 2018, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
 
-----------------------------------------------------------------------------------------
  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
'''
# ----------------------------------------------------------------------------------------

import sys
import os
import warnings
from functools import partial

import multiprocessing
import pickle
import glob

import numpy as np


from matplotlib import MatplotlibDeprecationWarning
#warnings.filterwarnings("error",category=MatplotlibDeprecationWarning)
import matplotlib.gridspec as gridspec 
import matplotlib.pyplot as plt

from astropy.io import fits
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.io.fits.verify import VerifyWarning

from astropy import log
log.setLevel('WARNING') # Disable pesky INFO flags from aplpy

# For debugging purposes: turn all warnings into errors
#warnings.filterwarnings("error", category=AstropyDeprecationWarning)
warnings.filterwarnings("ignore", category=RuntimeWarning)

# Also filter stupid GAIA warnings
warnings.filterwarnings("ignore", module='astropy.io.votable.tree')

from astropy.visualization import (PercentileInterval, AsymmetricPercentileInterval,
                                   AsinhStretch, LinearStretch, ImageNormalize)
import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.units import Quantity
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.time import Time

from astroquery.gaia import Gaia

from photutils import DAOStarFinder
from photutils import CircularAperture

# Import brutifus-specific tools
from . import brutifus_tools as bifus_t
from . import brutifus_cof as bifus_cof
from . import brutifus_elf
from . import brutifus_plots as bifus_p
from . import brutifus_red as bifus_red
from . import brutifus_metadata as bifus_m
from .brutifus_metadata import __version__

# ---------------------------------------------------------------------------------------- 
def run_adjust_WCS(fn_list, params, suffix = None, name_in = None, name_out = None, 
                   pmin = 10., pmax = 99.0, max_offset = 5.):
   ''' Queries Gaia to check the WCS validity.
   
   :param fn_list: The dictionary containing all filenames created by brutifus.
   :type fn_list: dict
   :param params:  The dictionary containing all paramaters set by the user. 
   :type params: dict
   :param suffix: The tag of this step, to be used in all files generated for rapid id.
   :type suffix: string
   :param name_in: name tag to identify which cube to use to run the routine 
   :type name_in: string
   :param name_out: name tag to identify which cube comes out of the routine
   :type name_out: string
   :param pmin: lower percentile for the plot
   :type pmin: float
   :param pmax: upper percentile for the plot
   :type pmax: float
   :param max_offset: max offset between GAIA and the data ... in PIXELS!!!
   :type max_offset: float 
   
   
   :return: The updated dictionary of filenames.
   :rtype: dictionary
   
   .. note:: Proper motions for the GAIA entries are propagated to the DATE-OBS of the data.
             Source identification in white-light image based on :class:`photutils.DAOStarFinder`
             
   .. todo:: This routine only works if the correct match is the closest neighbor. 
             It needs to be reinforced in case of large offsets (and/or crowded fields).
   
   '''
   
   if params['verbose']:
      print('-> Adjusting the cube WCS using Gaia.')
   
   # Get the data
   [[lams, data, error], [header0, header_data, header_error]] = \
      bifus_t.extract_cube(fn_list[name_in],params['inst'])
   
   # Build a white-light image
   wl_im = np.nansum(data, axis=0)
   
   # Get some of the statistics on this image
   wl_mean, wl_median, wl_std = sigma_clipped_stats(wl_im, sigma=3.0) 
   
   if params['verbose']:
      print('   Looking for stars in the field ...')
      
   # Source: https://photutils.readthedocs.io/en/stable/detection.html   
   # Get ready to look for stars in there
   daofind = DAOStarFinder(fwhm = 3.0, threshold = 5. * wl_std) 
   
   # Start the search
   sources = daofind(wl_im - wl_median) 
   
   # Create a list of sources in RA-Dec ?
   sources_list =  np.array([sources['xcentroid'],sources['ycentroid']])
   
   # Reshape this properly as Nx2
   sources_list = sources_list.T
   
   # Save this 2D white-light image 
   hdu0 = fits.PrimaryHDU(None, header0)
   hdu1 = fits.ImageHDU(wl_im)
   # Make sure the WCS coordinates are included as well
   hdu1 = bifus_t.hdu_add_wcs(hdu1, header_data)
   # Also include a brief mention about which version of brutifus is being used
   hdu1 = bifus_t.hdu_add_brutifus(hdu1, suffix)
      
   # Write the file!
   hdu = fits.HDUList(hdus = [hdu0, hdu1])
   
   fn_list['white_light'] = os.path.join(params['prod_loc'], 
                                         suffix+'_'+params['target']+'_white_light.fits') 
   hdu.writeto(fn_list['white_light'], overwrite = True)
   
   # Store the 2D WCS info for later
   wcs2d = WCS(hdu1)
                                       
   # Query Gaia, and show the marker location on the plot. 
   # First get the field "center" and radius
   if params['inst'] == 'MUSE':
      # For MUSE, assume the data is always of the same format (because it is).
      core_coord = SkyCoord(ra=header_data['CRVAL1'], dec=header_data['CRVAL2'],
                            unit =(u.deg, u.deg), frame='icrs')
      search_radius = 0.75 * max([np.abs(header_data['NAXIS1'] * header_data['CD1_1']*u.deg),
                                np.abs(header_data['NAXIS2'] * header_data['CD2_2']*u.deg)])
   else:
      raise Exception('Ouch! ... instrument not yet supported?!')
   
   if params['verbose']:
      print('   Querying GAIA ...')
   
   # Make it a sync search, because I don't think I need more than 2000 targets ...
   j = Gaia.cone_search_async(core_coord, search_radius, verbose = False)
   r = j.get_results()
   
   if len(r)>= 2000:
      warnings.warn('   Maximum number of GAIA entries returned. Some stars might be missing.') 
   
   # Some of these proper motions are NaN's ... need to deal with this before feeding them to 
   # SkyCoord
   r['pmra'].fill_value = 0
   r['pmdec'].fill_value = 0

   # Turn this into a "catalogue" of SkyCoord
   gaia_cat = SkyCoord(ra=r['ra'], dec=r['dec'], unit=(u.deg, u.deg), frame='icrs',
                       pm_ra_cosdec = r['pmra'].filled(),
                       pm_dec = r['pmdec'].filled(),
                       obstime = [Time(item, format='decimalyear') for item in r['ref_epoch']],
                       equinox = 'J2000',
                       distance = 1.e3 * u.pc, # Assumes the default distance
                       )
   
   # Propagate Gaia to the observation date!
   gaia_cat = gaia_cat.apply_space_motion(new_obstime = Time(header0['DATE-OBS'])) 
   
   # Clean this up, to only show the points within the image
   gaia_cat = gaia_cat[wcs2d.footprint_contains(gaia_cat)]
   
   # Transform this in pixel coords
   gaia_cat_pix = wcs2d.wcs_world2pix(gaia_cat.ra,gaia_cat.dec, 0)
   
   # Reshape it properly
   gaia_cat_pix = np.array(gaia_cat_pix).T
   
   # Now, loop through each point, find the closest neighbor
   dxys = []
   for source in sources_list:
      (this_d, this_delta, this_neighbor) = bifus_t.nearest_2dpoint(source,gaia_cat_pix)
      
      if this_d <= max_offset:
         dxys += [this_delta.tolist()]
   
   dxys = np.array(dxys)
   
   # Ok, for finding the best possible offsets, I have a lot of choice. 
   # Assuming I have enough points, let's just take the median along x and along y.
   
   dx = np.median(dxys[:,0])
   dy = np.median(dxys[:,1])
   
   # Let's make a plot to show how good I'm doing with getting a reasonable offset
   plt.close(99)
   fig99 = plt.figure(99, figsize = (6.94,7))
   
   gs = gridspec.GridSpec(1,1, 
                          height_ratios = [1],
                          width_ratios = [1], 
                          left = 0.12, right = 0.98, 
                          bottom = 0.12, top = 0.9, 
                          wspace = 0.1, hspace = 0.1)
   
   ax99 = plt.subplot(gs[0,0])  
   
   ax99.scatter(dxys[:,0],dxys[:,1],marker='.', color='k')
   ax99.axvline(dx)
   ax99.axhline(dy)
   
   ax99.set_xlim((-1.2*max_offset, 1.2*max_offset))
   ax99.set_ylim((-1.2*max_offset, 1.2*max_offset))
   
   circle = plt.Circle((0, 0), max_offset, facecolor='none', edgecolor='k')
   ax99.add_artist(circle)
   ax99.set_aspect('equal')
   ax99.set_xlabel('Delta x [pix]')
   ax99.set_ylabel('Delta y [pix]')
   ax99.set_title('Max. search radius: %.1f pixels' % (max_offset))
   
   fig99.savefig(os.path.join(params['plot_loc'], suffix + '_' + params['target'] + 
                                                           '_get-WCS-offsets.pdf'))
   plt.close(99)
   
   # Very well, now apply the corrections!
   # First to the white-light image
   fits.setval(fn_list['white_light'], 'CRPIX1', value = hdu1.header['CRPIX1'] - dx, ext=1)
   fits.setval(fn_list['white_light'], 'CRPIX2', value = hdu1.header['CRPIX2'] - dy, ext=1)
   
   # Then also to the fixed datacube
   hdu0 = fits.PrimaryHDU(None, header0)
   hdu1 = fits.ImageHDU(data)
   hdu2 = fits.ImageHDU(error)
        
   for hdu in [hdu1,hdu2]:
      # Make sure the WCS coordinates are included as well
      hdu = bifus_t.hdu_add_wcs(hdu, header_data)
      hdu = bifus_t.hdu_add_lams(hdu, header_data)
      # Also include a brief mention about which version of brutifus is being used
      hdu = bifus_t.hdu_add_brutifus(hdu, suffix)
      
      # And update the WCS value
      hdu.header['CRPIX1'] = hdu.header['CRPIX1'] - dx
      hdu.header['CRPIX2'] = hdu.header['CRPIX2'] - dy
      
   # Add the filename to the dictionary of filenames
   fn_list[name_out] = os.path.join(params['prod_loc'], 
                                    suffix+'_'+params['target']+'_wcs-corr.fits')
   hdu = fits.HDUList(hdus=[hdu0, hdu1, hdu2])  
   hdu.writeto( fn_list[name_out], overwrite=True)
   
   # Plot the white-light image to show how good I am doing ...
   (fig, ax, ofn) = bifus_p.make_2Dplot(fn_list['white_light'], ext = 1, 
                                        ofn = os.path.join(params['plot_loc'],
                                                           suffix + '_' + params['target'] + 
                                                           '_check-WCS.pdf'),
                                       pmin = pmin, pmax = pmax,
                                       cmap = None, 
                                       cblabel = None,
                                       #scalebar = params['scalebar'],
                                       )   
   
   # And add them to the plot
   ax.scatter(gaia_cat.ra, gaia_cat.dec, marker=bifus_p.crosshair(pa=45), color='darkorange', 
              s=100, 
              transform=ax.get_transform('world'))
   # ax.scatter(gaia_cat_pix[0], gaia_cat_pix[1], marker='+', color='darkorange', s=50)
   
   # Show the stars I found in the data ... these are in pix coordinates ...
   # -> always "right" irrespective of the WCS!
   ax.scatter(sources['xcentroid'], sources['ycentroid'], color='crimson', s = 50, 
              marker = 'o', facecolor = 'none')
                    
   # Save the updated plot
   fig.savefig(ofn)                               
                                    
   return fn_list
   
# ---------------------------------------------------------------------------------------- 
def run_crude_snr_maps(fn_list, params, suffix = None, name_in = None, do_plot = False, 
                       zcorr_lams = False):
   '''
   This function computes the SNR maps for the continuum and Ha (or other line). 
   It also creates a map of spaxels with any signal at all.
   The resulting maps are saved to a fits file with full header and WCS coordinates.
    
   :Args:
      fn_list: dictionary
               The dictionary containing all filenames created by brutifus.
      params: dictionary
              The dictionary containing all paramaters set by the user. 
      suffix: string [default: None]
              The tag of this step, to be used in all files generated for rapid id.
      name_in: string
               name tag to identify which cube to use to run the routine
      do_plot: bool [default: True]
               Whether to make a plot of the SNR maps or not.
      zcorr_lams: bool
                  True -> correct the wavelength range for the target redshift.
        
   :Returns:
      fn_list: dictionary
               The updated dictionary of filenames. 
                 
   :Notes:
      This function is a "first guess" of the SNR for latter use in the code. A more
      reliable measure of the SNR for the emission line should be computed after they
      have been fitted.
   '''
    
   if params['verbose']:
      print('-> Computing the SNR maps.')

   # Get the data
   [[lams, data, error], [header0, header_data, header_error]] = \
      bifus_t.extract_cube(fn_list[name_in],params['inst'])
    
   # Take redshift into account ?
   if zcorr_lams: 
      lams /= params['z_target']+ 1.
   
   # Prepare a storage list 
   snrs = []
   
   # Here, hide those stupid warnings about all-NaNs slices
   #with warnings.catch_warnings():
      #warnings.simplefilter("ignore", category=RuntimeWarning)
   
   for r in params['snr_ranges']:
      
         # The signal
         if r[-1] == 'c':
            s = np.nanmedian(data[(lams>=r[0])*(lams<=r[1]),:,:], axis = 0) 
         
         elif r[-1] == 'e':
            s = np.nanmax(data[(lams>=r[0]) * (lams<=r[1]), :,:], axis = 0)
      
         else:
            raise Exception('S/N calculation type unknown: %s' % r[-1])
      
         # The noise = STD over the range -> NOT strictly correct for high S/N stars !!!
         n = np.nanstd(data[(lams>=r[0])*(lams<=r[1]),:,:], axis = 0) 
      
         # Compute S/N
         snr = s/n
      
         # Ensure this is always >= 0
         snr[snr<0] = 0
      
         # Store it for later
         snrs += [snr]
   
   # And create a map with just spaxels that have any data (i.e. have been observed).
   anything = np.ones_like(data[0,:,:])
   anything[np.all(np.isnan(data),axis=0)] = np.nan
       
   # Very well, now let's create a fits file to save this as required.
   hdu0 = fits.PrimaryHDU(None,header0)
   hdus = [hdu0]
   
   # Add the regions where I have data ... always useful !
   hdu = fits.ImageHDU(anything)
   # Make sure the WCS coordinates are included as well
   hdu = bifus_t.hdu_add_wcs(hdu, header_data)
   # Also include a brief mention about which version of brutifus is being used
   hdu = bifus_t.hdu_add_brutifus(hdu, suffix)
   # For reference, also include the line/region this maps are based on
   hdu.header['B_SNRANG'] = ('%.1f-%.1f' % (lams[0],lams[-1]), 'spectral range (A) used for SNR')
   hdu.header['B_SNTYPE'] = ('x', 'binary map: NaN= no data, 1 = valid spectra')
   
   hdus += [hdu]
   
   # Now also loop through all the maps I have created
   for (i,r) in enumerate(params['snr_ranges']):
   
      hdu = fits.ImageHDU(snrs[i])
      
      # Make sure the WCS coordinates are included as well
      hdu = bifus_t.hdu_add_wcs(hdu, header_data)
      
      # Also include a brief mention about which version of brutifus is being used
      hdu = bifus_t.hdu_add_brutifus(hdu, suffix)
      
      # For reference, also include the line/region this maps are based on
      hdu.header['B_SNRANG'] = ('%.1f-%.1f' % (r[0],r[1]), 'spectral range (A) used for SNR')
      hdu.header['B_SNTYPE'] = ('%s' % r[-1], '"c"ontinuum, or "e"mission')
      
      hdus += [hdu]
        
   # Write the file!
   hdu = fits.HDUList(hdus=hdus)
   
   fn_list['snr_maps'] = os.path.join(params['prod_loc'], 
                                      suffix+'_'+params['target']+'_snr-maps.fits')
   
   hdu.writeto(fn_list['snr_maps'], overwrite = True)
        
   # Make some plots if requested
   if do_plot:
      
      # First, plot the region with any signal at all
      bifus_p.make_2Dplot(fn_list['snr_maps'],
                             1, 
                             os.path.join(params['plot_loc'],
                                          suffix + '_' + params['target'] + 
                                          '_valid_spectra.pdf'),
                             #contour_fn = None, 
                             #contour_ext = None,
                             #contour_levels = [None], 
                             vmin = 0, vmax = 1.5,
                             cmap = None, 
                             cblabel = None,
                             scalebar = params['scalebar'],
                            )  
   
      # Alright, let's take out the big guns ...    
      for (i,r) in enumerate(params['snr_ranges']):
             
         bifus_p.make_2Dplot(fn_list['snr_maps'],
                             i+2, 
                             os.path.join(params['plot_loc'],
                                          suffix + '_' + params['target'] + 
                                          '_snr_%.1f-%.1f.pdf' % (r[0],r[1])),
                             vmin = 0, vmax = 50,
                             cmap = 'alligator_r',
                             cblabel = r"S/N ('%s') %.1f\AA-%.1f\AA" % (r[-1],r[0],r[1]),
                             scalebar = params['scalebar'],
                            )                                         
   return fn_list

# ---------------------------------------------------------------------------------------- 
def run_plot_BW(fn_list, params, suffix=None, 
                name_in = None,
                bands = [[7500.,9300.],], 
                conts = [None,],
                stretches = ['arcsinh',],
                stretch_plims = [[10.,99.5],],
                stretch_vlims = [[None,None],],
                gauss_blurs = [None,]
               ):   
   ''' 
   This function is designed to make some B&W images from slices in the cube.
    
   :Args:
      fn_list: dictionary
               The dictionary containing all filenames created by brutifus.
      params: dictionary
              The dictionary containing all paramaters set by the user. 
      suffix: string [default: None]
              The tag of this step, to be used in all files generated for rapid id.
      name_in: string
               name tag to identify which cube to use to run the routine
      bands: list of list of list [default: [[7500,9300],],
           A list of wavelength pairs.
      conts: list of list of list [default: [[None,None],],
             A list of wavelength pairs channels associated continuum regions.         
      stretches: list of string [default: ['arcsinh']]
                 The stretches to apply to the data, e.g. 'linear', 'log', 'arcsinh'.
      stretch_plims: list of list of floats [default: [[10.,99.5],]
                     The limiting percentiles for the plot
      stretch_vlims: list of list of floats [default: [[10.,99.5],]]
                     The limiting values for the plot (superseeds stretch_plims)
      gauss_blurs: list of scalars
                   Radius of gaussian blur, None for nothing
             
   :Returns:
      fn_list: dictionary
               The updated dictionary of filenames. 
   ''' 
    
   # First, I need to generate 3 individual fits files for each band
   if params['verbose']:
      print('-> Creating some B&W images.')
    
   # Get the data
   [[lams, data, error], [header0, header_data, header_error]] = \
      bifus_t.extract_cube(fn_list[name_in],params['inst'])
   
   # Step 1: Construct individual images for each band
   for (i,band) in enumerate(bands):
        
      # get the data out
      scidata = data[(lams>=band[0])*(lams<=band[1]),:,:]
      # Get the continuum
      if not(conts[i] is None):
         contdata = data[(lams>=conts[i][0])*(lams<=conts[i][1]),:,:]   
      else:
         contdata = 0
         
      # Subtract the continuum          
      scidata = np.sum(scidata,axis=0) - np.sum(contdata, axis=0)
         
      fn = os.path.join(params['plot_loc'],'BW_tmp.fits')
      hdu = fits.PrimaryHDU(scidata)
      # Add the wcs info
      hdu = bifus_t.hdu_add_wcs(hdu, header_data)
      outfits = fits.HDUList([hdu])
      outfits.writeto(fn,overwrite=True)
            
        
      ofn = os.path.join(params['plot_loc'], suffix+'_'+params['target']+'_'+name_in+
                         '_BW_%i-%i.pdf' % (band[0],band[0]))
                                                                             
      # Great, I am now ready to call the plotting function
      bifus_p.make_2Dplot(fn, ext = 0, ofn = ofn, 
                           stretch = stretches[i],
                           pmin = stretch_plims[i][0],
                           pmax = stretch_plims[i][1],
                           vmin = stretch_vlims[i][0],
                           vmax = stretch_vlims[i][1],
                           gauss_blur = gauss_blurs[i],
                           cmap = None,
                           #title = r'\smaller %s-%s\,\AA\' % (band[0],band[1]),
                           scalebar = params['scalebar'])
    
      # And remember to delete all the temporary files
      os.remove(fn)    
    
   return fn_list    
     
# ---------------------------------------------------------------------------------------- 
def run_plot_RGB(fn_list, params, suffix=None, 
                        name_in = None,
                        bands = [[[7500.,9300.],[6000.,7500.],[4800.,6000.]],], 
                        conts = [[None, None, None],],
                        stretches = [['arcsinh','arcsinh','arcsinh'],],
                        stretch_plims = [[10.,99.5,10.,99.5,10.,99.5],],
                        stretch_vlims = [[None,None,None,None,None,None],],
                        gauss_blurs = [[None,None,None,],]
                       ):   
   ''' 
   This function is designed to make some RGB images from slices in the raw cube.
    
   :Args:
      fn_list: dictionary
               The dictionary containing all filenames created by brutifus.
      params: dictionary
              The dictionary containing all paramaters set by the user. 
      suffix: string [default: None]
              The tag of this step, to be used in all files generated for rapid id.
      name_in: string
               name tag to identify which cube to use to run the routine
      bands: list of list of list [default: [[[7500,9300],[6000,7500],[4800,6000]],],
           A list of wavelength pairs limiting the R, G & B channels.
      conts: list of list of list [default: [[None,None,None],],
             A list of wavelength pairs limiting the R, G & B channels associated continuum
             regions.         
      stretches: list of string [default: ['log']]
                 The stretches to apply to the data, e.g. 'linear', 'log', 'arcsinh'.
      stretch_plims: list of list of floats [default: [[10.,99.5,10.,99.5,10.,99.5],]]
                     The limiting percentiles for the plot, as 
                     [pmin_r, pmax_r, pmin_g, pmax_g, pmin_b, pmax_b]
      stretch_vlims: list of list of floats [default: [[10.,99.5,10.,99.5,10.,99.5],]]
                     The limiting values for the plot (superseeds stretch_plims), as
                     [vmin_r, vmax_r, vmin_g, vmax_g, vmin_b, vmax_b]
      use_egal_dered: bool [default: False]
                      If available, whether to use the line fluxes corrected for 
                      extragalactic attentuation.
             
   :Returns:
      fn_list: dictionary
               The updated dictionary of filenames. 
   ''' 
    
   # First, I need to generate 3 individual fits files for each band
   if params['verbose']:
      print('-> Creating some RGB images.')
    
   # Get the data
   [[lams, data, error], [header0, header_data, header_error]] = \
      bifus_t.extract_cube(fn_list[name_in],params['inst'])
   
   # Step 1: Construct individual images for each band
   for (i,band) in enumerate(bands):
        
      fns = []
       
      for (j,b) in enumerate(band):
         # get the data out
         scidata = data[(lams>=b[0])*(lams<=b[1]),:,:]
         if not(conts[i][j] is None):
            contdata = data[(lams>=conts[i][j][0])*(lams<=conts[i][j][1]),:,:]   
         else:
            contdata = 0
                
         scidata = np.sum(scidata,axis=0) - np.sum(contdata, axis=0)
         
         fn = os.path.join(params['plot_loc'],'RGB_tmp_%i.fits' % j)
         fns.append(fn)
         hdu = fits.PrimaryHDU(scidata)
         # Add the wcs info
         hdu = bifus_t.hdu_add_wcs(hdu, header_data)
         outfits = fits.HDUList([hdu])
         outfits.writeto(fn,overwrite=True)
            
        
      ofn = os.path.join(params['plot_loc'], suffix+'_'+params['target']+'_'+name_in+
                         '_RGB_%i-%i_%i-%i_%i-%i.pdf' % (band[0][0],band[0][1],
                                                         band[1][0],band[1][1],
                                                         band[2][0],band[2][1]))
                                                                             
      # Great, I am now ready to call the plotting function
      bifus_p.make_RGBplot(fns, ofn, 
                           stretch = stretches[i],
                           plims = stretch_plims[i], 
                           vlims = stretch_vlims[i],
                           gauss_blur = gauss_blurs[i],
                           title = r'\smaller R: %s-%s\,\AA\  G: %s-%s\,\AA\  B: %s-%s\,\AA' %  
                                   (band[0][0],band[0][1],
                                    band[1][0], band[1][1],
                                    band[2][0], band[2][1]),
                           scalebar = params['scalebar'])
    
      # And remember to delete all the temporary files
      for f in fns:
         os.remove(f)    
    
   return fn_list    
   
# ----------------------------------------------------------------------------------------
def run_sky_sub(fn_list, params, suffix = None, name_in = None, name_out = None):
   '''Performs  a manual sky subtraction, when the observation strategy was flawed.
    
   :Args:
      fn_list: dictionary
               The dictionary containing all filenames created by brutifus.
      params: dictionary
              The dictionary containing all paramaters set by the user. 
      suffix: string [default: None]
              The tag of this step, to be used in all files generated for rapid id.
      name_in: string
               name tag to identify which cube to use to run the routine
      name_out: string
               name tag to identify which cube come out of the routine
    
   :Returns:
      fn_list: dictionary
               The updated dictionary of filenames.               
   '''
   
   if params['verbose']:
      print('-> Subtracting the sky ...')
    
   # Get the data
   [[lams, data, error], [header0, header_data, header_error]] = \
      bifus_t.extract_cube(fn_list[name_in],params['inst'])
  
   # Assemble a 2D map of the sky spaxels
   sky_spaxels = np.zeros_like(data[0,:,:])
   cys,cxs = np.mgrid[0:header_data['NAXIS2'],0:header_data['NAXIS1']]
   
   #  Flag all the (sky) spaxels within the user-defined aperture
   for sr in params['sky_regions']:
         
      if len(sr) == 3: # Circular aperture
         
         d = np.sqrt( (cxs-sr[0])**2 + (cys-sr[1])**2)
         sky_spaxels[d<=sr[2]] = 1
      
      elif len(sr) ==4: # Square aperture
         
         sky_spaxels[sr[1]:sr[1]+sr[3]+1, sr[0]:sr[0]+sr[2]+1] = 1
   
   # Very well, assemble the sky spectrum now
   sky_spec = np.array([np.nanmedian(plane[sky_spaxels==1]) for plane in data[:]])
   
   # Make a descent plot of the sky spectrum
   plt.close(1)
   plt.figure(1, figsize=(14.16,4))
   gs = gridspec.GridSpec(1,1, height_ratios=[1], width_ratios=[1],
                          left=0.08, right=0.98, bottom=0.18, top=0.98, wspace=0.02, hspace=0.03)
   
   ax1 = plt.subplot(gs[0,0])
   
   ax1.semilogy(lams,sky_spec, 'k-', drawstyle='steps-mid', lw=0.8)
   
   ax1.set_xlim((lams[0],lams[-1]))
   ax1.set_xlabel('Observed wavelength [\AA]')
   ax1.set_ylabel(bifus_m.ffmt[params['inst']]['funit'])
   
   plt.savefig(os.path.join(params['plot_loc'], suffix+'_'+params['target']+'_skyspec.pdf'))
   plt.close(1)
   
   # In addition, also make a plot showing all the sky locations
   
   # Start by assembling a white light image
   wl_im = np.nansum(data, axis = 0)
   
   # Very well, now let's create a fits file to save this as required.
   hdu0 = fits.PrimaryHDU(None, header0)
   # Add the white light image
   hdu1 = fits.ImageHDU(wl_im)
   # Make sure the WCS coordinates are included as well
   hdu1 = bifus_t.hdu_add_wcs(hdu1, header_data)
   # Also include a brief mention about which version of brutifus is being used
   hdu1 = bifus_t.hdu_add_brutifus(hdu1, suffix)
        
   # Write the file!
   hdu = fits.HDUList(hdus=[hdu0, hdu1])
   
   fn_list['wl_im'] = os.path.join(params['prod_loc'], 
                                   suffix+'_'+params['target']+'_wl-im.fits')
   hdu.writeto(fn_list['wl_im'], overwrite = True)
   
   # Image filename
   ofn = os.path.join(params['plot_loc'], 
                                   suffix+'_'+params['target']+'_sky-regions.pdf')
   
   # Make the base figure
   (fig, ax1, ofn) = bifus_p.make_2Dplot(fn_list['wl_im'], ext = 1, 
                                         ofn = ofn,
                                         stretch = 'linear',
                                         pmin = 5,
                                         pmax = 90,
                                         )
   
   # Show all the sky spaxels
   #sky_spaxels[sky_spaxels ==0] = np.nan
   #ax1.imshow(sky_spaxels, origin= 'lower', cmap='winter', vmin = 0, vmax = 1, alpha=0.5)
   ax1.contour(sky_spaxels, levels=[0.5], colors = ['w'], 
               linewidths = [0.75], origin = 'lower')
   # TODO: have the contours follow the pixel edges
   
   
   fig.savefig(ofn)
   
   # Now get started with the actual sky subtraction
   
   # Turn the spectrum into a cube 
   sky_cube = np.repeat(np.repeat(sky_spec[:,np.newaxis], 
                                  header_data['NAXIS2'], axis=1)[:,:,np.newaxis],
                        header_data['NAXIS1'], axis=2)
   
   # Perform the sky subtraction
   data -= sky_cube
   
   # Here, I assume that the sky correction is "perfect" (i.e. it comes with no errors)
   # Since it's just a subtraction, I don't need to alter the error cube
   
   # Export the new datacube
   hdu0 = fits.PrimaryHDU(None, header0)
   hdu1 = fits.ImageHDU(data)
   hdu2 = fits.ImageHDU(error)
        
   for hdu in [hdu1,hdu2]:
      # Make sure the WCS coordinates are included as well
      hdu = bifus_t.hdu_add_wcs(hdu, header_data)
      hdu = bifus_t.hdu_add_lams(hdu, header_data)
      # Also include a brief mention about which version of brutifus is being used
      hdu = bifus_t.hdu_add_brutifus(hdu, suffix)

   
   # Add the filename to the dictionary of filenames
   fn_list[name_out] = os.path.join(params['prod_loc'], 
                                    suffix+'_'+params['target']+'_skysub-cube.fits')
   hdu = fits.HDUList(hdus=[hdu0, hdu1, hdu2])  
   hdu.writeto( fn_list[name_out], overwrite=True)
        
   
   return fn_list
   
# ----------------------------------------------------------------------------------------
def run_gal_dered(fn_list, params, suffix = None, do_plot = False,
                  name_in = None, name_out = None):
   '''Corrects for Galactic extinction, given the Ab and Av extinction.
   
   This function derives the Alambda value for any wavelength, and corrects the data to
   correct for our "local" extinction. Intended for extragalactic sources.
   
   :Args:
      fn_list: dictionary
               The dictionary containing all filenames created by brutifus.
      params: dictionary
              The dictionary containing all paramaters set by the user. 
      suffix: string [default: None]
              The tag of this step, to be used in all files generated for rapid id.
      do_plot: bool [default: True]
               Whether to make a plot of the Alambda correction applied or not.
      name_in: string
               name tag to identify which cube to use to run the routine
      name_out: string
               name tag to identify which cube comes out of the routine
   :Returns:
      fn_list: dictionary
               The updated dictionary of filenames. 
                 
   :Notes:
        To reproduce the approach from NED, use the Ab and Av value for you object from 
        there, and set curve='f99', rv=3.1 in the params file.   
   '''
    
   if params['verbose']:
      print('-> Correcting for Galactic extinction.')
    
   # Get the data
   [[lams, data, error], [header0, header_data, header_error]] = \
      bifus_t.extract_cube(fn_list[name_in],params['inst'])
    
   # Compute Alambda
   alams = bifus_red.alam(lams, params['Ab'],params['Av'], curve=params['gal_curve'],
                           rv=params['gal_rv'])
    
   # Compute the flux correction factor
   etau = bifus_red.galactic_red(lams,params['Ab'],params['Av'], 
                                 curve=params['gal_curve'],
                                 rv=params['gal_rv'])
    
   # Correct the data and the propagate the errors
   data *= etau[:, np.newaxis, np.newaxis]
   error *= (etau[:, np.newaxis, np.newaxis])**2
    
   # Save the datacubes
   hdu0 = fits.PrimaryHDU(None,header0)
   hdu1 = fits.ImageHDU(data)
   hdu2 = fits.ImageHDU(error)
        
   for hdu in [hdu1,hdu2]:
      # Make sure the WCS coordinates are included as well
      hdu = bifus_t.hdu_add_wcs(hdu, header_data)
      hdu = bifus_t.hdu_add_lams(hdu, header_data)
      # Also include a brief mention about which version of brutifus is being used
      hdu = bifus_t.hdu_add_brutifus(hdu, suffix)
      # Add the line reference wavelength for future references
      hdu.header['BR_AB'] = (params['Ab'], 'Extinction in B')
      hdu.header['BR_AV'] = (params['Av'], 'Extinction in V')
            
   # Add the filename to the dictionary of filenames
   fn_list[name_out] = os.path.join(params['prod_loc'], 
                                    suffix+'_'+params['target']+'_gal-dered_cube.fits')
                                        
   hdu = fits.HDUList(hdus=[hdu0,hdu1,hdu2])
   hdu.writeto(fn_list[name_out], overwrite=True)
        
   # Do I want to save a plot ?
   if do_plot:
      ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+
                            '_gal_Alambda_corr.pdf')
      bifus_p.make_galred_plot(lams,alams,etau,params['Ab'],params['Av'], ofn)
    
   return fn_list
    
# ----------------------------------------------------------------------------------------
def run_fit_continuum(fn_list, params, suffix=None, name_in = None,
                      start_row = None, end_row = None,
                      method='lowess'):
   ''' 
   This function fits the continuum in the datacube. It is designed to use multiprocessing
   to speed things up on good computers. It deals with the data columns-per-columns, and
   can be restarted mid-course, in case of a crash. 
   
   :Args:
      fn_list: dictionary
               The dictionary containing all filenames created by brutifus.
      params: dictionary
              The dictionary containing all paramaters set by the user. 
      suffix: string [default: None]
              The tag of this step, to be used in all files generated for rapid id.
      name_in: string
               name tag to identify which cube to use to run the routine
   :Returns:
      fn_list: dictionary
               The updated dictionary of filenames. 
   '''
    
   # For the bad continuum: run a LOWESS filter from statsmodel.nonparametric:
   # sm.nonparametric.smoothers_lowess.lowess(spec,lams,frac=0.05, it=5)
    
   # Rather than launch it all at once, let's be smart in case of problems. I'll run
   # the fits row-by-row with multiprocessing (hence up to 300cpus can help!), and save
   # between each row.
   
    # Get the data
   [[lams, data, error], [header0, header_data, header_error]] = \
      bifus_t.extract_cube(fn_list[name_in],params['inst'])
   
   # I also need to load the SNR cube to know where I have data I want to fit
   '''
   ### DISABLED FOR NOW - ONLY MEANINGFUl WITH MORE THAN 1 CONT. FIT. TECHNIQUE ###
   try:
      print(os.path.join(params['prod_loc'],fn_list['snr_maps']))
      hdu = fits.open(os.path.join(params['prod_loc'],fn_list['snr_maps']))
      is_nonnan = hdu[1].data
      hdu.close()
   except:
      raise Exception('Failed to load the crude SNR cube. Did you run the step ?')
   '''
   
   # Get some info about the cube
   nrows = header_data['NAXIS1']
   
   # Where do I start and end the processing ?
   if start_row is None:
      start_row = 0
   if end_row is None:
      end_row = nrows-1

   # Ok, what do I want to do ?
   if method == 'lowess':
      if params['verbose']:
         print('-> Starting the continuum fitting using the LOWESS approach.')
            
      fit_func = partial(bifus_cof.lowess_fit, lams = lams, 
                         frac = params['lowess_frac'], it = params['lowess_it'])
      # Note here the clever use of the partial function, that turns the lowess_fit
      # function from something that takes 4 arguments into something that only takes 1
      # argument ... thus perfect for the upcoming "map" functions !
   
   else:
      raise Exception('Continuum fitting method "%s" not supported.' % (method))

   # Very well, let's start the loop on rows. If the code crashes/is interrupted, you'll
   # loose the current row. Just live with it.
   for row in np.linspace(start_row, end_row, end_row - start_row + 1).astype(int):   
      
      # Alright, now deal with the spaxels outside the user-chosen SNR range.
      # Replace them with nan's
      good_spaxels = np.ones((header_data['NAXIS2']))
      
      ''' 
      ### DISABLED FOR NOW - ONLY MEANINGFUl WITH MORE THAN 1 CONT. FIT. TECHNIQUE ###
      if params[method+'_snr_min']:
         good_spaxels[snr_cont[:,row] < params[method+'_snr_min']] = np.nan
      if params[method+'_snr_max']:
         good_spaxels[snr_cont[:,row] >= params[method+'_snr_max']] = np.nan
      '''
      
      # Build a list of spectra to be fitted
      specs = [data[:,i,row] * good_spaxels[i] for i in range(header_data['NAXIS2'])]
      errs = [error[:,i,row] * good_spaxels[i] for i in range(header_data['NAXIS2'])]
        
      # Set up the multiprocessing pool of workers
      if params['multiprocessing']:
         # Did the user specify a number of processes to use ?
         if type(params['multiprocessing']) == np.int:
            nproc = params['multiprocessing']
                
         else: # Ok, just use them all ...
            nproc = multiprocessing.cpu_count()
                        
         pool = multiprocessing.Pool(processes = nproc, 
                                     initializer = bifus_t.init_worker())    

         if params['verbose']:
            sys.stdout.write('\r   Fitting spectra in row %2.i, %i at a time ...' % 
                                 (row,nproc))
            sys.stdout.flush()

         # Launch the fitting ! Make sure I deal with KeyBoard Interrupt properly
         # Only a problem for multiprocessing. For the rest of the code, whatever.
         try:   
            conts = pool.map(fit_func, specs)
         except KeyboardInterrupt:
            print(' interrupted !')
            # Still close and join properly
            pool.close()
            pool.join()
            sys.exit('Multiprocessing continuum fitting interrupted at row %i'% row)
         else: # If all is fine
            pool.close()
            pool.join()  
              
              
      else: # just do things 1-by-1
         if params['verbose']:
            sys.stdout.write('\r   Fitting spectra in row %2.i, one at a time ...' % 
                                 row)
            sys.stdout.flush()
                                
         conts = map(fit_func, specs)

      # Here, I need to save these results. Pickle is fast and temporary,
      # until I re-build the entire cube later on. It also allows for better 
      # row-by-row flexibility.
      
      print(' ') # Need this to deal with the stdout mess
      
      if not(os.path.isdir(params['tmp_loc'])):
         print('   !Requested storage location does not exist! Creating it ...')
         print('    => %s' % params['tmp_loc'])
         os .mkdir(params['tmp_loc'])
      
      fn = os.path.join(params['tmp_loc'],
                        suffix+'_'+params['target']+'_'+method+'_row_'+
                        str(np.int(row)).zfill(4)+'.pkl')
      file = open(fn,'wb')
      pickle.dump(conts,file)
      file.close()
        
   # And add the generic pickle filename to the dictionary of filenames
   fn_list[method+'_pickle'] = suffix+'_'+params['target']+'_'+method+'_row_'
        
   print('   Fitting completed !')
    
   return fn_list
   
# ----------------------------------------------------------------------------------------
def run_make_continuum_cube(fn_list, params, suffix=None,
                            method='lowess'):   
   ''' 
   This function is designed to construct a "usable and decent" datacube out of the
   mess generated by the continuum fitting function, i.e. out of the many pickle 
   files generated.
   
   :Args:
      fn_list: dictionary
               The dictionary containing all filenames created by brutifus.
      params: dictionary
              The dictionary containing all paramaters set by the user. 
      suffix: string [default: None]
              The tag of this step, to be used in all files generated for rapid id.
      method: string
              which fits to gather ?
   
   '''
    
   if params['verbose']:
      print('-> Constructing the datacube for the continuum fitting (%s).' % method)
    
   # Get the raw data open, to know how big things are ...
   [[lams, data, error], [header0, header_data, header_error]] = \
      bifus_t.extract_cube(fn_list['raw_cube'],params['inst']) 
        
   # Get some info
   nrows = header_data['NAXIS1']
   
   # Prepare a continuum cube structure
   cont_cube = np.zeros_like(data) * np.nan
                                          
   # Loop through the rows, and extract the results. 
   # Try to loop through everything - in case this step was run in chunks.
   for row in range(0,nrows):
      progress = 100. * (row+1.)/nrows
      sys.stdout.write('\r   Building cube [%5.1f%s]' % 
                         (progress,'%'))
      sys.stdout.flush()

      fn = os.path.join(params['tmp_loc'],
                        fn_list[method+'_pickle']+str(np.int(row)).zfill(4)+'.pkl')

      if os.path.isfile(fn):
         # Very well, I have some fit here. Let's get them back
         myfile = open(fn,'rb')
         conts = pickle.load(myfile)  
         myfile.close() 
        
         # Mind the shape
         if method=='lowess':
            # I directly saved the continuum as an array just get it out.
            cont_cube[:,:,row] = np.array(conts).T 
         
         else:
            raise Exception(' Continuum fitting method "%s" unknown.' % (method))
                                     
   # Very well, now let's create a fits file to save this as required.
   hdu0 = fits.PrimaryHDU(None,header0)
   hdu1 = fits.ImageHDU(cont_cube)
   if method == 'lowess':
      hdu2 = fits.ImageHDU(np.zeros_like(cont_cube)) # Keep all the errors to 0 for now.
   else:
      raise Exception('What errors do you have in your fitted sky ???')
   
   # Make sure the WCS coordinates are included as well
   for hdu in [hdu1,hdu2]:
      hdu = bifus_t.hdu_add_wcs(hdu,header_data)
      hdu = bifus_t.hdu_add_lams(hdu,header_data)
   
      # Also include a brief mention about which version of brutifus is being used
      hdu = bifus_t.hdu_add_brutifus(hdu,suffix)
      
   hdu = fits.HDUList(hdus=[hdu0, hdu1, hdu2])
   fn_list[method+'_cube'] = os.path.join(params['prod_loc'],
                                          suffix+'_'+params['target']+'_'+method+'.fits')
   hdu.writeto(fn_list[method+'_cube'], overwrite=True)
    
   print(' ') # Required to clear the stdout stuff
    
   return fn_list

# ----------------------------------------------------------------------------------------
def run_subtract_continuum(fn_list, params, suffix = None, name_in = None,
                           name_out = None, method = None):
   ''' 
   This function will subtract the sky from a given cube.
   
   :Args:
      fn_list: dictionary
               The dictionary containing all filenames created by brutifus.
      params: dictionary
              The dictionary containing all paramaters set by the user. 
      suffix: string [default: None]
              The tag of this step, to be used in all files generated for rapid id.
      name_in: string
               name tag to identify which cube to use to run the routine
      name_out: string
               name tag to identify which cube comes out of the routine
   '''

   


   if params['verbose']:
      print('-> Subtracting the sky continuum (%s).' % method)
   
   # Get the data open
   [[lams, data, error], [header0, header_data, header_error]] = \
      bifus_t.extract_cube(fn_list[name_in],params['inst']) 
   
   # Get the skycube open
   [[skylams, skydata, skyerror], [skyheader0, skyheader_data, skyheader_error]] = \
      bifus_t.extract_cube(fn_list[method+'_cube'],params['inst']) 
   
   
   # Since this is a subtraction, and I assume no error on the sky, the errors remain unchanged
   data -= skydata 
     
   # And save this to a new fits file
   hdu0 = fits.PrimaryHDU(None,header0)
   hdu1 = fits.ImageHDU(data)
   hdu2 = fits.ImageHDU(error)
      
   for hdu in [hdu1,hdu2]:
      # Make sure the WCS coordinates are included as well
      hdu = bifus_t.hdu_add_wcs(hdu,header_data)
      hdu = bifus_t.hdu_add_lams(hdu,header_data)
      # Also include a brief mention about which version of brutifus is being used
      hdu = bifus_t.hdu_add_brutifus(hdu,suffix)
      
   hdus = fits.HDUList(hdus=[hdu0,hdu1,hdu2])
   
   fn_list[name_out] = os.path.join(params['prod_loc'], 
                                    suffix+'_'+params['target']+'_'+method+
                                    '-contsub-cube.fits')
   hdus.writeto(fn_list[name_out], overwrite=True)

   return fn_list 
   
# ----------------------------------------------------------------------------------------
def run_fit_elines(fn_list, params, suffix=None, start_row = None, end_row = None, 
                   ):
   ''' 
   This function fits the emission lines in the datacube, after subtracting the continuum
   derived using LOWESSS or PPXF. It is designed to use multiprocessing to speed
   things up on good computers. It deals with the data columns-per-columns, and 
   can be restarted mid-course, in case of a crash. 
   '''
   
   raise Exception('Check and update this function!!!')
   
   
   nlines = len(params['elines'].keys())
   if params['verbose']:
      print('-> Starting the emission line fitting for %2.i line(s).' % nlines)
    
   # Rather than launch it all at once, let's be smart in case of problems. I'll run
	# the fits row-by-row with multiprocessing (hence up to 300cpus can help!), and save
	# between each row.
	
	# First, load the datacube to be fitted. Use the galactic deredened one if it exists.
   if fn_list['galdered_cube'] is None:
      hdu = fits.open(os.path.join(params['data_loc'],params['data_fn']))
   else:
      hdu = fits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
   header0 = hdu[0].header
   data = hdu[1].data
   header1 = hdu[1].header
   error = hdu[2].data
   header2 = hdu[2].header
   hdu.close()
    
   # Build the wavelength array
   lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
   # Calculate the outer bin edges for the spectrum
   be = np.append(lams-header1['CD3_3']/2.,lams[-1]+header1['CD3_3']/2.)
    
   # I also need to load the SNR cube for the spaxel selection
   hdu = fits.open(os.path.join(params['prod_loc'],
                                   fn_list['snr_cube']))
   snr_cont = hdu[1].data                               
   snr_elines = hdu[2].data
   hdu.close()
    
   # I also need the continuum cubes
   if fn_list['cont_mix_cube']:
      fn = os.path.join(params['prod_loc'],fn_list['cont_mix_cube'])
      if os.path.isfile(fn):
         hdu = fits.open(fn)
         cont_mix_cube = hdu[1].data
         hdu.close()
        
      # Very well, now, perform the continuum subtraction.
      data -= cont_mix_cube
    
   # Get some info about the cube
   nrows = header1['NAXIS1']
   if start_row is None:
      start_row = 0
   if end_row is None:
	   end_row = nrows-1 
        
   fit_func = partial(brutifus_elf.els_mpfit, lams=lams, be=be, params=params)
	# Note here the clever use of the partial function, that turns the els_mpfit
	# function from something that takes 4 arguments into something that only takes 1
	# argument ... thus perfect for the upcoming "map" functions !
    
   # Very well, let's start the loop on rows. If the code crashes/is interrupted, you'll
	# loose the current row. Just live with it.
   for row in np.linspace(start_row,end_row, end_row-start_row+1):   
	    
		# Alright, now deal with the spaxels outside the user-chosen SNR range.
		# Replace them with nan's
      good_spaxels = np.ones((header1['NAXIS2']))
      if params['elines_snr_min']:
         good_spaxels[snr_elines[:,row]<params['elines_snr_min']] = np.nan
      if params['elines_snr_max']:
         good_spaxels[snr_elines[:,row]>params['elines_snr_max']] = np.nan
		
		# Build a list of spectra to be fitted
      specerrs = [[data[:,i,row] * good_spaxels[i],
                  error[:,i,row]* good_spaxels[i]] for i in range(header1['NAXIS2'])]
        
		# Set up the multiprocessing pool of workers
      if params['multiprocessing']:
            # Did the user specify a number of processes to use ?
            if type(params['multiprocessing']) == np.int:
                nproc = params['multiprocessing']
                pool = multiprocessing.Pool(processes = nproc, 
                                            initializer = brutifus_tools.init_worker())
                
            else: # Ok, just use them all ...
                nproc = multiprocessing.cpu_count()
                pool = multiprocessing.Pool(processes=None, 
                                            initializer = brutifus_tools.init_worker())
			
            if params['verbose']:
                sys.stdout.write('\r   Fitting spectra in row %2.i, %i at a time ...' % 
                                 (row,nproc))
                sys.stdout.flush()
			
		      # Launch the fitting ! Make sure to deal with KeyBoard Interrupt properly
		      # Only a problem for multiprocessing. For the rest of the code, whatever.
            try:   
                els = pool.map(fit_func, specerrs)
            except KeyboardInterrupt:
                print(' interrupted !')
                # Still close and join properly
                pool.close()
                pool.join()
                sys.exit('Multiprocessing line fitting interrupted at row %i'% row)
            else: # If all is fine
                pool.close()
                pool.join()  
              
      else: # just do things 1-by-1
            if params['verbose']:
                sys.stdout.write('\r   Fitting spectra in row %2.i, one at a time ...' % 
                                 row)
                sys.stdout.flush() 
                                
            els = map(fit_func, specerrs)
	    
	   # Here, I need to save these results. Pickle could be fast and temporary,
	   # Until I then re-build the entire cube later on ? Also allow for better
	   # row-by-row flexibility.
      fn = os.path.join(params['tmp_loc'],
                          suffix+'_'+params['target']+'_elines_row_'+
                          str(np.int(row)).zfill(4)+'.pkl')
      file = open(fn,'w')
      pickle.dump(els,file)
      file.close()
   	
   	# Add the generic filename to the dictionary of filenames
      fn_list['elines_pickle'] = suffix+'_'+params['target']+'_elines_row_'
   	     
   print(' done !')
    
   return fn_list
# ----------------------------------------------------------------------------------------

def run_make_elines_cube(fn_list, params, suffix=None):   
    ''' 
    This function is designed to construct a "usable and decent" datacube out of the
    mess generated by the emission line fitting function, i.e. out of the many pickle 
    files generated.
    '''
    
    raise Exception('Check and update this function!!!')
    
    # First, load the original datacube. I need to know how much stuff was fitted.
    if fn_list['galdered_cube'] is None:
        hdu = fits.open(os.path.join(params['data_loc'],params['data_fn']))
    else:
        hdu = fits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    # Build the wavelength array
    lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    # Calculate the outer bin edges for the spectrum
    be = np.append(lams-header1['CD3_3']/2.,lams[-1]+header1['CD3_3']/2.) 
    
    nrows = header1['NAXIS1']
       
    # How many emission lines were fitted ?
    nlines = len(params['elines'].keys()) 
       
    # Very well, what do I want to extract ?
    # 1) A "full emission line spectrum"
    # 2) For each line, a flux map, an intensity map, a velocity map and a dispersion map,
    # and h3 and h4 maps
    # Let's get to it.
    elines_fullspec_cube = data * np.nan # That way, I keep the nan's in the raw data
    # 7 params: Flux, I, v, sigma, h3, h4, SNR
    elines_params_cube = np.zeros((7*nlines,header1['NAXIS2'],header1['NAXIS1']))*np.nan
    # And the associated errors !
    elines_params_err = np.zeros((7*nlines,header1['NAXIS2'],header1['NAXIS1']))*np.nan
    # Also save the status from mpfit. Could be useful for sorting the good from the bad.
    elines_fit_status = np.zeros((header1['NAXIS2'],header1['NAXIS1']))*np.nan
    
    if params['verbose']:
        print('-> Constructing the datacube for the emission line fitting parameters.')
       
    # Loop through the rows, and extract the results. 
    # Try to loop through everything - in case this step was run in chunks.
    for row in range(0,nrows):
        progress = 100. * (row+1.)/nrows
        sys.stdout.write('\r   Building cubes [%5.1f%s]' % 
                         (progress,'%'))
        sys.stdout.flush()

        fn = os.path.join(params['tmp_loc'],
                          fn_list['elines_pickle']+str(np.int(row)).zfill(4)+'.pkl')

        if os.path.isfile(fn):
            # Very well, I have some fit here. Let's get them back
            myfile = open(fn,'r')
            ms = pickle.load(myfile)  
            myfile.close() 
        
            # Get all the parameters in a big array
            ps = [item.params for item in ms]
            errs = [item.perror for item in ms]
            stats = [item.status for item in ms]
            
            # Here, I need to make sure the ps and errs array have a decent shape, even 
            # when the fit failed. Also store the STD that comes out of mpfit (NOT the variance!) 
            ps = [np.zeros_like(ps[0])*np.nan if not(item.status in [1,2,3,4]) else ps[j] 
                                                            for (j,item) in enumerate(ms)]
            
            errs = [np.zeros_like(ps[0])*np.nan if not(item.status in [1,2,3,4]) else 
                                                 errs[j] for (j,item) in enumerate(ms)]
            
            
            # Fill the corresponding datacube
            # First, what are my indices ? I need to keep one slice empty for the SNR of 
            # each line. Get the indices where I will thus store the fit params.
            ind = np.array([np.arange(6)+7*j for j in range(nlines)]).reshape(nlines*6)
            
            elines_params_cube[ind,:,row] = np.array(ps).T
            elines_params_err[ind,:,row] = np.array(errs).T
            elines_fit_status[:,row] = np.array(stats).T
            
            # now, reconstruct the full emission line spectrum
            elines_specs = np.array([brutifus_elf.els_spec(lams,p,be=be, 
                                                        method=params['line_profile'],
                                                        inst=params['inst']) for p in ps])
       
            # Fill the corresponding cube 
            elines_fullspec_cube[:,:,row] = elines_specs.T
            
    
    # Now, for each line, the first of these slices is the reference wavelength. 
    # Replace this by the total flux instead ! And make some plots if requested.
    for (k,key) in enumerate(np.sort(params['elines'].keys())):
        
        # Calculate the sigma of the fit, in A (incl. instrument dispersion,etc ...)
        # as well as the associated error.
        # WARNING: unlike the data, these are the 1-sigma error coming out of mpfit.
        # Mind the **2 !
        zlams = elines_params_cube[7*k] * (1.+ elines_params_cube[7*k+2] / c )
        zlams_err = elines_params_cube[7*k] * elines_params_err[7*k+2]/c # The 1-std err
        sigma_obs_A = brutifus_elf.obs_sigma(elines_params_cube[7*k+3],zlams, 
                                          inst=params['inst'], 
                                          in_errs=[elines_params_err[7*k+3],zlams_err] )
        
        # Compute the line flux
        elines_params_cube[7*k,:,:] = np.sqrt(2*np.pi) * elines_params_cube[7*k+1,:,:] * \
                                      sigma_obs_A[0]                           
        # What about the error ?
        elines_params_err[7*k,:,:] = \
            np.abs(elines_params_cube[7*k,:,:]) * \
            np.sqrt(elines_params_err[7*k+1,:,:]**2/elines_params_cube[7*k+1,:,:]**2 +
                    sigma_obs_A[1]**2/sigma_obs_A[0]**2)           
    
    
    # Now for each line, calculate the SNR
    # For this, I need the continuum cubes
    if fn_list['cont_mix_cube']:
        fn = os.path.join(params['prod_loc'],fn_list['cont_mix_cube'])
        hdu = fits.open(fn)
        cont_mix_cube = hdu[1].data
        hdu.close()
    else:
        cont_mix_cube = np.zeros_like(data)
        
    # Compute the residual cube
    residuals = data - cont_mix_cube - elines_fullspec_cube
    
    cont_range = params['cont_range']
    lams0 = lams / (params['z_target'] + 1)
    
    # Compute the noise over the region of interest defined by the user
    noise = np.nanstd(residuals[(lams0>=cont_range[0])*(lams0<=cont_range[1]),:,:],axis=0) 
    
    for (k, key) in enumerate(np.sort(params['elines'].keys())):
        elines_params_cube[7*k+6,:,:] = elines_params_cube[7*k+1,:,:]/noise
       
    
    # Export the cube for each emission line parameters as a multi-extension fits file        
    # Do the same for the errors - i.e. params and errors are in two distinct cubes
    for (e,epc) in enumerate([elines_params_cube,elines_params_err]):
        hdu0 = fits.PrimaryHDU(None,header0)
    
        hdus = [hdu0]
        # Use the sorted keys, to ensure the same order as the fit parameters
        for (k,key) in enumerate(np.sort(params['elines'].keys())):
            hduk = fits.ImageHDU(epc[7*k:7*k+7,:,:])
            # Make sure the WCS coordinates are included as well
            hduk = brutifus_tools.hdu_add_wcs(hduk,header1)
            # Also include a brief mention about which version of brutifus is being used
            hduk = brutifus_tools.hdu_add_brutifus(hduk,suffix)
            # Add the line reference wavelength for future references
            hduk.header['BR_REFL'] = (params['elines'][key][0][0], 'reference wavelength')
            hduk.header['BR_CKEY'] = ('F,I,v,sigma,h3,h4,SNR','Content of the cube planes')
            
            hdus.append(hduk)
            
        hdu = fits.HDUList(hdus=hdus)
        fn_out = os.path.join(params['prod_loc'],
                              suffix+'_'+params['target']+'_elines_'+
                              ['params','perror'][e]+'.fits')
        hdu.writeto(fn_out, overwrite=True)
        
        # Add the filename to the dictionary of filenames
        fn_list['elines_'+['params','perror'][e]+'_cube'] = suffix+'_'+params['target']+\
                                                            '_elines_'+\
                                                            ['params','perror'][e]+'.fits'
                   
    # Very well, now let's also create a fits file to save the full emission line spectrum
    # as required.
    hdu0 = fits.PrimaryHDU(None,header0)
    hdu1 = fits.ImageHDU(elines_fullspec_cube)
    # Make sure the WCS coordinates are included as well
    hdu1 = brutifus_tools.hdu_add_wcs(hdu1,header1)
    hdu1 = brutifus_tools.hdu_add_lams(hdu1,header1)
    # Also include a brief mention about which version of brutifus is being used
    hdu1 = brutifus_tools.hdu_add_brutifus(hdu1,suffix)
    hdu = fits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_elines_fullspec.fits')
    hdu.writeto(fn_out, overwrite=True)
    
    # Add the filename to the dictionary of filenames
    fn_list['elines_spec_cube'] = suffix+'_'+params['target']+'_elines_fullspec.fits'
    
    # And finally, the plot with the fit status for each spaxel, to know the bad from the 
    # not-so-bad,
    hdu0 = fits.PrimaryHDU(None,header0)
    hdu1 = fits.ImageHDU(elines_fit_status)
    # Make sure the WCS coordinates are included as well
    hdu1 = brutifus_tools.hdu_add_wcs(hdu1,header1)
    # Also include a brief mention about which version of brutifus is being used
    hdu1 = brutifus_tools.hdu_add_brutifus(hdu1,suffix)
    hdu = fits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_elines_mpfit_status.fits')
    hdu.writeto(fn_out, overwrite=True)
    
    # Add the filename to the dictionary of filenames
    fn_list['elines_fit_status'] = suffix+'_'+params['target']+'_elines_mpfit_status.fits'
    
    print(' ')
    return fn_list
# ----------------------------------------------------------------------------------------

def run_plot_elines_cube(fn_list, params, suffix=None, vrange=None, 
                         sigrange=None):   
    '''Creates some plots for the emission lines: F, I, v and sigma.
        
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutifus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        vrange: list of int [default: None]
                If set, the range of the colorbar for the velocity plot.
        sigrangre: list of int [default: None]        
                If set, the range of the colorbar for the velocity dispersion plot.
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames.
    ''' 
    
    raise Exception('Check and update this function!!!')
       
    if params['verbose']:
        print('-> Making some nifty plots from the emission line fitting output.')
       
    fn = os.path.join(params['prod_loc'], fn_list['elines_params_cube'])
    
    # Open the file
    hdu = fits.open(fn)
    header0 = hdu[0].header

    # Create a temporary FITS file
    fn_tmp = os.path.join(params['tmp_loc'],suffix+'_tmp.fits')

    # Do this for each emission line fitted
    for (k,key) in enumerate(np.sort(params['elines'].keys())):
    
        # Because of the way aplpy works, I need to stored each "image" in its own fits
        # file. I don't want to keep them, so let's just use a temporary one, get the plot
        # done, and remove it. Not ideal, but who cares ?
        
        this_header = hdu[k+1].header
        this_data = hdu[k+1].data
        
        # Make single pretty plots for Flux, Intensity, velocity and velocity dispersion                                        
        for (t,typ) in enumerate(['F','I','v','sigma','h3', 'h4', 'SNR']):
            
            # Take into account the SNR from the user
            if typ != 'SNR':
                this_data[t][this_data[-1] <params['elines'][key][0][3]] = np.nan
            
            # Create a dedicated HDU    
            tmphdu = fits.PrimaryHDU(this_data[t])
            # Add the WCS information
            tmphdu = brutifus_tools.hdu_add_wcs(tmphdu,this_header)
            tmphdu.writeto(fn_tmp, overwrite=True)
            
            # Now plot it 
            this_ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+
                                                        '_eline-'+key+'-'+
                                                        np.str(this_header['BR_REFL'])+
                                                        '_'+
                                                        typ+'.pdf') 
            
            # For the velocity fields, set the vmin and vmax
            if t == 0:
                my_vmin = None
                my_vmax = None
                my_cmap = 'alligator'
                my_stretch = 'arcsinh'
                my_label = r'F [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$]'
                my_cbticks = [125,250,500,1000,2000,4000]
            elif t == 1:
                my_vmin = None
                my_vmax = None
                my_cmap = 'alligator'
                my_stretch = 'arcsinh'
                my_label = r'I [10$^{-20}$ erg s$^{-1}$ cm$^{-2}$]'
                my_cbticks = [50,100,200,400,800,1600]
            elif t == 2 and vrange:
                my_vmin = vrange[0]
                my_vmax = vrange[1]
                my_cmap = 'alligator'
                my_stretch = 'linear'
                my_label = r'$v$ [km s$^{-1}$]'
                my_cbticks = None
            elif t ==3 and sigrange:
                my_vmin = sigrange[0]
                my_vmax = sigrange[1]
                my_cmap = 'alligator'
                my_stretch = 'linear'
                my_label = r'$\sigma_{tot}$ [km s$^{-1}$]'
                my_cbticks = None
            else:
                my_vmin = None
                my_vmax = None
                my_cmap = None
                my_stretch = 'linear'
                my_label = ''
                my_cbticks = None
                                                            
            bifus_p.make_2Dplot(fn_tmp,ext=0, ofn=this_ofn, contours=False, 
                                    vmin=my_vmin, vmax = my_vmax, cmap=my_cmap,
                                    stretch = my_stretch, cblabel=my_label, 
                                    cbticks = my_cbticks,
                                    scalebar = params['scalebar'],)                                   
    
    # Delete the temporary fits file
    os.remove(fn_tmp)
    # Don't forget to close the initial hdu ...
    hdu.close()
    
    # And also create a plot of the fit status, to see if anything weird happened
    fn = os.path.join(params['prod_loc'], fn_list['elines_fit_status'])
    ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+
                                                        '_eline_mpfit_status.pdf') 
    bifus_p.make_2Dplot(fn,ext=1, ofn=ofn, contours=False, vmin=-16, vmax = 8,
                           cmap='alligator',stretch='linear', 
                           cbticks=[-16,0,1,2,3,4,5,6,7,8], cblabel='mpfit status',
                           scalebar = None)
    
    return fn_list
# ----------------------------------------------------------------------------------------

def run_inspect_fit(fn_list,params, suffix=None, irange=[None,None], vrange=[None,None]):   
   '''Setup the interactive inspection of the fit results.
    
   :Args:
      fn_list: dictionary
               The dictionary containing all filenames created by brutifus.
      params: dictionary
               The dictionary containing all paramaters set by the user. 
      suffix: string [default: None]
               The tag of this step, to be used in all files generated for rapid id.
      irange: list of int [default: [None,None]]
               The range of the colorbar for the intensity plot.
      vrange: list of int [default: [None,None]]    
               The range of the colorbar for the velocity dispersion plot.

   :Returns:
      fn_list: dictionary
               The updated dictionary of filenames.
   '''
   
   
   if params['verbose']:
      print('-> Starting the interactive inspection of the fitting.')
    
   # Load the different files I need  
   # Raw data:    
   if fn_list['galdered_cube'] is None:
      hdu = fits.open(os.path.join(params['data_loc'],params['data_fn']))
   else:
      hdu = fits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
   header0 = hdu[0].header
   data = hdu[1].data
   header1 = hdu[1].header
   error = hdu[2].data
   header2 = hdu[2].header
   hdu.close()
    
   lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    
   # The Lowess continuum fit
   if not(fn_list['lowess_cube'] is None):
      fn = os.path.join(params['prod_loc'],fn_list['lowess_cube'])
      hdu = fits.open(fn)
      lowess = hdu[1].data
      hdu.close()  
   else:
      lowess = np.zeros_like(data)*np.nan
        
   # The ppxf continuum fit
   if not(fn_list['ppxf_cube'] is None):
      fn = os.path.join(params['prod_loc'],fn_list['ppxf_cube'])
      hdu = fits.open(fn)
      ppxf = hdu[1].data
      hdu.close()  
   else:
      ppxf = np.zeros_like(data)*np.nan
    
   # The continuum mix datacube
   if not(fn_list['cont_mix_cube'] is None):
      fn = os.path.join(params['prod_loc'],fn_list['cont_mix_cube'])
      hdu = fits.open(fn)
      cont_mix = hdu[1].data
      hdu.close()  
   else:
      cont_mix = np.zeros_like(data)*np.nan
    
       
   # The elines fit
   if not(fn_list['elines_spec_cube'] is None):
      fn = os.path.join(params['prod_loc'],fn_list['elines_spec_cube'])
      hdu = fits.open(fn)
      elines = hdu[1].data
      hdu.close() 
   else:
      elines = np.zeros_like(data)*np.nan 
    
   # Also open the map and vmap to be displayed. 
   # TODO: give user the choice of image on the right
   if not(fn_list['elines_params_cube'] is None):
      fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
      hdu = fits.open(fn)
      map = hdu[1].data[0]
      vmap = hdu[1].data[2]
      hdu.close()
   else:
      map = np.zeros_like(data[0])*np.nan
      vmap = np.zeros_like(data[0])*np.nan    
    
   # A filename if the user decides to save anything.
   my_ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+'_fit_inspection_')

   # Launch the interactive plot
   bifus_p.inspect_spaxels(lams,data,lowess,ppxf,cont_mix, elines, map, vmap, 
                           irange, vrange,
                           ofn = my_ofn)    
        
   return fn_list
# ----------------------------------------------------------------------------------------

def run_find_structures(fn_list,params, suffix=None, interactive_mode=True, 
                        automatic_mode=True):   
    ''' 
    This function is designed to identify structures (e.g. HII regions) in the data from
    a 2D image (i.e. an line intensity map), and save them to a pickle file. When
    interactive_mode=True, the user can manually refine the selection. Set
    automatic_mode=False to skip the automatic detection.
    '''
    
    raise Exception('Check and update this function!!!')
    
    if params['verbose']:
        print('-> Starting the semi-automated procedure for structure identification.')
    
    # Where am I going to save the apertures information ?
    fn_ap = os.path.join(params['prod_loc'],suffix+'_'+params['target']+'_ap_list.pkl')
    
    # Add the filename to the dictionary of filenames
    fn_list['ap_list'] = suffix+'_'+params['target']+'_ap_list.pkl'
    
    # Do we want to build apertures based on multiple maps ? Loop, one after the other.
    for key in params['ap_map_lines']:
        print('   Loading eline %s' %key)
        # Is the file already there ? Do I want to overwrite it ? And/or use its content ?
        if os.path.isfile(fn_ap):
            print('    ')
            print('   Existing aperture list (%s)' % fn_ap.split('/')[-1])
            print('   Type [a] to load this aperture list and edit it manually')
            print('        [b] to start from scratch, (file will be overwritten!)')
            print('        [c] to do nothing, and continue to the next step')
            while True:
                letter = raw_input()
            
                if letter in ['a','b','c']:
                    break
                else:
                    print('        [%s] unrecognized. Try again.' % letter)
        else:
            letter = 'b'
            
        # Open the file, load the list of apertures
        if letter == 'a':
            f = open(fn_ap, 'r')
            start_aps = pickle.load(f)
            f.close()
        elif letter =='b':
            start_aps = None
        elif letter =='c':
            continue
    
        # Now, open the elines param datacube, and extract the Flux map I want to detect
        # stuctures from.

        fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
        hdu = fits.open(fn)
        header0 = hdu[0].header
        plane = np.sort(params['elines'].keys()).tolist().index(key)
        data = hdu[plane+1].data[0]
        hdu.close()
        
        # Launch the aperture finding routine
        apertures = bifus_p.build_ap_list(data, start_aps = start_aps, 
                                              radius = params['ap_radius'],
                                              automatic_mode = automatic_mode,
                                              interactive_mode = interactive_mode,
                                              lam = params['elines'][key][0][0],
                                              save_plot = os.path.join(params['plot_loc'],
                                                                       suffix+'_'+
                                                                       params['target']+
                                                                       '_ap_list_'),
                                             )    
        
        # Only if the user wants to save the apertures, do it
        if apertures:
            # Save the results for later use
            f = open(fn_ap, 'w')
            pickle.dump(apertures,f)
            f.close()
        
    return fn_list
# ----------------------------------------------------------------------------------------

def run_make_ap_cube(fn_list, params, suffix=None, do_plot=True):   
    ''' 
    This function is designed to make a cube from a series of apertures (x,y,rs).
    For compativilty with spaxels-by-spaxels analysis codes (incl.brutifus), make the cube
    the same size as the original, and repleace each spectra in a given aperture by the
    total aperture spectra. Spaxels outside any apertures are nan's. 
    Assigned spaxels to one aperture only, in order of decreasing flux peak. This makes 
    the data redondant, but will allow for a rapid and direct processing of the resulting
    cube by brutifus.
    '''        
    
    raise Exception('Check and update this function!!!')
    
    if params['verbose']:
        print('-> Constructing the cube with the integrated '+\
              'aperture spectra.')
    
    # Very well, where is the aperture file ?
    fn_ap = os.path.join(params['prod_loc'],fn_list['ap_list'])
    f = open(fn_ap, 'r')
    start_aps = pickle.load(f)
    f.close()
    
    xs,ys,rs = zip(*start_aps)
    
    # I will also need the Flux map from the strongest line - use that set for the 
    # velocity reference of the line fitting
    for key in params['elines'].keys():
        if params['elines'][key][0][0] == params['ref_dv_line']:
            ref_key = key
    
    # Very well, now load the corresponding flux map
    fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
    hdu = fits.open(fn)
    fheader0 = hdu[0].header
    plane = np.sort(params['elines'].keys()).tolist().index(ref_key)
    flux = hdu[plane+1].data[0]
    fheader1 = hdu[plane+1].header
    hdu.close()
    
    # I also need to load the raw data cube. Here I ONLY use the raw cube - and NOT the
    # one possibly deredened for galactic extinction.
    hdu = fits.open(os.path.join(params['data_loc'],params['data_fn']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    # Get all the peak intensities associated with each aperture. Needed for sorting them.
    fs = flux[ys,xs]
    
    # Sort them in decreasing order of peak intensity
    sort_index = np.argsort(fs)[::-1]
    
    # Now, construct a map where each pixel contains the number of the region it belongs 
    # to. Starting from 0.
    ap_map = np.zeros_like(flux) * np.nan
    (ny,nx) = np.shape(flux)
    indices = np.mgrid[0:ny,0:nx]
    
    # Also construct an aperture spectra cube
    ap_spec_cube = np.zeros_like(data) * np.nan
    ap_spec_err = np.zeros_like(data) * np.nan
    
    # Loop through each ap. Hopefully not too slow ...
    for (i,ind) in enumerate(sort_index):
        progress = 100. * (i+1.)/len(sort_index)
        sys.stdout.write('\r   Dealing with aperture %i [%5.1f%s]' % 
                         (i,progress,'%'))
        sys.stdout.flush()
        x = xs[ind]
        y = ys[ind]
        r = rs[ind]
        # Find all spaxels with the ap radius 
        # Don't do anything fancy, just measure the distance to each spaxel center.
        in_ap = (indices[1]-x)**2+(indices[0]-y)**2 <= r**2
        
        # Assign each spaxel (not yet assigned to another brighter feature) the 
        # corresponding ap number.
        ap_map[in_ap * np.isnan(ap_map)] = i

        #For this aperture, sum all spaxels into master aperture spectra, and fill the 
        # cube. Avoid the nans (e.g. mosaic edges, etc ...)

        spec = np.nansum(data[:,indices[0][ap_map==i],indices[1][ap_map==i]],axis=1)
        err = np.nansum(error[:,indices[0][ap_map==i],indices[1][ap_map==i]],axis=1)

        # Loop through each spaxel in the aperture. There MUST be a smarter way, but I
        # can't figure it out. Should not be too costly time-wise anyway ...
        for k in range(len(indices[1][ap_map==i])):
            xi = indices[1][ap_map==i][k]
            yi = indices[0][ap_map==i][k]
            ap_spec_cube[:,yi,xi] = spec
            ap_spec_err[:,yi,xi]  = err
        
    # All done. Save the aperture map to a fits file.
    hdu0 = fits.PrimaryHDU(None,fheader0)
    hdu1 = fits.ImageHDU(ap_map)
    # Make sure the WCS coordinates are included as well
    hdu1 = brutifus_tools.hdu_add_wcs(hdu1,fheader1)
    # Also include a brief mention about which version of brutifus is being used
    hdu1 = brutifus_tools.hdu_add_brutifus(hdu1,suffix)
    hdu = fits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_ap_map.fits')
    hdu.writeto(fn_out, overwrite=True)    
    
    # Add this filename to the dictionary of filenames
    fn_list['ap_map'] = suffix+'_'+params['target']+'_ap_map.fits'
    
    # Make a plot of the apertures ?
    if do_plot:
        this_ofn = os.path.join(params['plot_loc'],
                          suffix+'_'+params['target']+'_ap_map.pdf')
        bifus_p.make_2Dplot(fn_out,ext=1, ofn=this_ofn, contours=False, 
                                    vmin=0, vmax = len(xs), cmap='alligator',
                                    stretch = 'linear', 
                                    cblabel=r'Aperture idendification number', 
                                    cbticks = None,
                                    scalebar = None,
                                )  
    
    # And also save the aperture spectral cube to a fits file
    hdu0 = fits.PrimaryHDU(None,header0)
    hdu1 = fits.ImageHDU(ap_spec_cube)
    hdu2 = fits.ImageHDU(ap_spec_err)
    # Make sure the WCS coordinates are included as well
    for hdu in [hdu1,hdu2]:
        hdu = brutifus_tools.hdu_add_wcs(hdu,header1)
        hdu = brutifus_tools.hdu_add_lams(hdu,header1)
        # Also include a brief mention about which version of brutifus is being used
        hdu = brutifus_tools.hdu_add_brutifus(hdu,suffix)
        
    hdu = fits.HDUList(hdus=[hdu0,hdu1,hdu2])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_ap_spec_cube.fits')
    hdu.writeto(fn_out, overwrite=True)   
      
    # Add this filename to the dictionary of filenames
    fn_list['ap_spec_cube'] = suffix+'_'+params['target']+'_ap_spec_cube.fits'
                                     
    print(' ')
    
    return fn_list           
# ----------------------------------------------------------------------------------------

def run_extragal_dered(fn_list, params, suffix=None, do_plot=True):   
    '''Corrects the line fluxes for the extragalactic reddening using Ha and Hb.
    
    This function returns a corrected set of line fluxes, as well as the associated Av 
    map.
     
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutifus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        do_plot: bool [default: True]
                 Whether to make a plot of the Av map or not.
        
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
                 
    :Notes:
        Should I add some info about each curve here ? 
    '''  
    
    raise Exception('Check and update this function!!!')
    
    if params['verbose']:
        print('-> Correcting for the extragalactic attenuation.')
    
    # Open the raw data cube
    if fn_list['galdered_cube'] is None:
        hdu = fits.open(os.path.join(params['data_loc'],params['data_fn']))
    else:
        hdu = fits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    nlines = len(params['elines'].keys())

    # Import the Halpha and Hbeta maps, and mask any spaxel not within the proper SNR range
    fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
    hdu = fits.open(fn)
    for (k,key) in enumerate(np.sort(params['elines'].keys())):
        if params['elines'][key][0][0] == ha:
            ha_map = hdu[k+1].data[0]
            ha_map[hdu[k+1].data[-1] < params['elines'][key][0][3]] = np.nan
        elif params['elines'][key][0][0] == hb:
            hb_map = hdu[k+1].data[0]
            hb_map[hdu[k+1].data[-1] < params['elines'][key][0][3]] = np.nan
    hdu.close()


    hahb = ha_map / hb_map
    # Construct the Av map, just because I can
    av = brutifus_red.hahb_to_av(hahb, params['hahb_0'], curve = params['egal_curve'], 
                              rv = params['egal_rv'], rva = params['egal_rva'])

    # Now, for each emission line fitted, let's correct the flux:
    # A storage structure for the emission lines. Keep the un-reddened flux as well.
    drelines_params_cube = np.zeros((8*nlines,header1['NAXIS2'],header1['NAXIS1']))*np.nan
    drelines_perror_cube = np.zeros((8*nlines,header1['NAXIS2'],header1['NAXIS1']))*np.nan

    # Load the current line parameters, loop through, and save to a new - bigger - array.
    # Also take care of the error
    fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])
    fn_err = os.path.join(params['prod_loc'],fn_list['elines_perror_cube'])
    hdu = fits.open(fn)
    hdu_err = fits.open(fn_err)
    fheader0 = hdu[0].header
    fheader0_err = hdu_err[0].header
    
    # Loop through each emission line
    for (k,key) in enumerate(np.sort(params['elines'].keys())):
        drelines_params_cube[8*k+1:8*(k+1)] = hdu[k+1].data
        drelines_perror_cube[8*k+1:8*(k+1)] = hdu_err[k+1].data
        # And add the de-reddened line flux in the first layer.
        # Note: the reddening correction is based on the de-redshifted line wavelength,
        # because this happens in the rest-frame of the target !
        this_lam = params['elines'][key][0][0]
        drelines_params_cube[8*k] = hdu[k+1].data[0] * \
                                     brutifus_red.extragalactic_red(this_lam, hahb, 
                                                                 params['hahb_0'],
                                                                 curve = params['egal_curve'], 
                                                                 rv = params['egal_rv'],
                                                                 rva = params['egal_rva'])
        # Assume the reddening correction is error free ... sigh...
        drelines_perror_cube[8*k] = hdu_err[k+1].data[0] * \
                                     (brutifus_red.extragalactic_red(this_lam, hahb, 
                                                                 params['hahb_0'],
                                                                 curve = params['egal_curve'], 
                                                                 rv = params['egal_rv'],
                                                                 rva = params['egal_rva']))**2                                                        
    fheader1 = hdu[1].header
    fheader1_err = hdu_err[1].header
    hdu.close()
    hdu_err.close()
    
    
    for (e,epc) in enumerate([drelines_params_cube, drelines_perror_cube]):
        # Ok, now save the data 
        hdu0 = fits.PrimaryHDU(None,fheader0)
    
        hdus = [hdu0]
        # Use the sorted keys, to ensure the same order as the fit parameters
        for (k,key) in enumerate(np.sort(params['elines'].keys())):
            hduk = fits.ImageHDU(epc[8*k:8*(k+1),:,:])
            # Make sure the WCS coordinates are included as well
            hduk = brutifus_tools.hdu_add_wcs(hduk,header1)
            # Also include a brief mention about which version of brutifus is being used
            hduk = brutifus_tools.hdu_add_brutifus(hduk,suffix)
            # Add the line reference wavelength for future references
            hduk.header['BR_REFL'] = (params['elines'][key][0][0], 'reference wavelength')
            hduk.header['BR_CKEY'] = ('dredF,F,I,v,sigma,h3,h4','Content of the cube planes')
            
            hdus.append(hduk)
            
        hdu = fits.HDUList(hdus=hdus)
        fn_out = os.path.join(params['prod_loc'],
                              suffix+'_'+params['target']+'_elines_'+
                              ['params','perror'][e]+'_dered.fits')
        hdu.writeto(fn_out, overwrite=True)
        
        # Add the filename to the dictionary of filenames
        fn_list['dered_elines_'+['params','perror'][e]] = suffix+'_'+\
                                                                  params['target']+\
                                                                  '_elines_'+\
                                                                  ['params','perror'][e]+\
                                                                  '_dered.fits'
                   
    # Very well, now let's also create a fits file to save the Av map
    # as required.
    hdu0 = fits.PrimaryHDU(None,header0)
    hdu1 = fits.ImageHDU(av)
    # Make sure the WCS coordinates are included as well
    hdu1 = brutifus_tools.hdu_add_wcs(hdu1,header1)
    # Also include a brief mention about which version of brutifus is being used
    hdu1 = brutifus_tools.hdu_add_brutifus(hdu1,suffix)
    hdu = fits.HDUList(hdus=[hdu0,hdu1])
    fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_Av.fits')
    hdu.writeto(fn_out, overwrite=True)
    
    # Add the filename to the dictionary of filenames
    fn_list['Av_map'] = suffix+'_'+params['target']+'_Av.fits'
    
    # If requested, save a plot for the Av map.
    if do_plot:
        
        this_ofn = os.path.join(params['plot_loc'],
                          suffix+'_'+params['target']+'_Av_map.pdf')
        bifus_p.make_2Dplot(fn_out,ext=1, ofn=this_ofn, contours=False, 
                                    vmin=0, vmax = 2, cmap='alligator',
                                    stretch = 'linear', 
                                    cblabel=r'A$_V$ [mag]', 
                                    cbticks = None,
                                    scalebar = params['scalebar'],
                                )
                                    
    return fn_list
# ----------------------------------------------------------------------------------------

def run_plot_elines_RGB(fn_list, params, suffix=None, 
                        mixes = [['[NII]','Ha','[OIII]'],], 
                        stretches = ['log',],
                        stretch_plims = [[10.,99.5,10.,99.5,10.,99.5],],
                        stretch_vlims = [[None,None,None,None,None,None],],
                        use_egal_dered = False,
                       ):   
    ''' 
    This function is designed to make some RGB images from emission lines.
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutifus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        mixes: list of list [default: [['[NII]','Ha','[OIII]']]]
             A list of line triplets indicatin the R, G & B channels.
        stretches: list of string [default: ['log']]
                   The stretches to apply to the data, e.g. 'linear', 'log', 'arcsinh'.
        stretch_plims: list of list of floats [default: [[10.,99.5,10.,99.5,10.,99.5],]]
                       The limiting percentiles for the plot, as 
                       [pmin_r, pmax_r, pmin_g, pmax_g, pmin_b, pmax_b]
        stretch_vlims: list of list of floats [default: [[10.,99.5,10.,99.5,10.,99.5],]]
                       The limtiing values for the plot (superseeds stretch_plims), as
                       [vmin_r, vmax_r, vmin_g, vmax_g, vmin_b, vmax_b]
        use_egal_dered: bool [default: False]
                        If available, whether to use the line fluxes corrected for 
                        extragalactic attentuation.
             
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
    ''' 
    
    raise Exception('Check and update this function!!!')
    
    # First, I need to generate 3 individual fits files with the line fluxes
    if params['verbose']:
        print('-> Creating some RGB emission line images.')
    
    # Very well, now I need to open the line fluxes.
    if use_egal_dered:
        fn = os.path.join(params['prod_loc'],fn_list['dered_elines_params'])
    else:
        fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])

    # Alright, now I need to extract the line fluxes
    # Open the file
    hdu = fits.open(fn)
    header0 = hdu[0].header
    header1 = hdu[1].header

    # All elines definition are inside elines_pyqz inside brutifus_metadata. 
    # For the ones the user requested, construct a flux map.
    elines_fluxes = {}
    for (m, mix) in enumerate(mixes):
        elines_fluxes[m] = {0:None, 1:None, 2:None}    

    # Look at each emission line fitted. Do I need it ?
    for (k,akey) in enumerate(np.sort(params['elines'].keys())):
    
        this_header = hdu[k+1].header
        this_data = hdu[k+1].data
                               
        lam = this_header['BR_REFL']
        
        for (m,mix) in enumerate(mixes):
            
            for (cc,col) in enumerate(mix): # Search all possible lines
                if lam in elines_pyqz[col]: # Wavelength match ? I found a line.
            
                    flux = this_data[0]
                
                    # Take into account the SNR
                    flux[this_data[-1]<params['elines'][akey][0][3]] = np.nan
            
                    # Already exists ? - then I found a line with multiple components
                    if not(elines_fluxes[m][cc] is None): 
                        elines_fluxes[m][cc] += flux # append the flux
                    else:
                        elines_fluxes[m][cc] = flux # Save the data for later
                       
    # Close the hdus
    hdu.close()
    
    # Very well, I now have all my line fluxes. Let's deal with each RGB set one at a time.
    for (m,mix) in enumerate(mixes):

        fns = []
        
        ofn = os.path.join(params['plot_loc'],
                 suffix+'_'+params['target']+'_RGB_%s-%s-%s.pdf' % (mix[0],mix[1],mix[2]))
                 
        # Now, let' create 3 independant temporary fits file with WCS coordinates
        for (cc,col) in enumerate(mix):
            fn = os.path.join(params['plot_loc'],'RGB_tmp_%i.fits' % cc)
            fns.append(fn)
            
            hdu = fits.PrimaryHDU(elines_fluxes[m][cc], header = header1)
            # Add the wcs info
            hdu = brutifus_tools.hdu_add_wcs(hdu, header1)
            outfits = fits.HDUList([hdu])
            outfits.writeto(fn,overwrite=True)
    
        # Great, I am now ready to call the plotting function
        bifus_p.make_RGBplot(fns, ofn,  stretch = stretches[m],
                                  plims = stretch_plims[m], vlims = stretch_vlims[m],
                                  title = r'RGB: %s %s %s' %  (mix[0],mix[1],mix[2]),
                                  scalebar = params['scalebar'])
    
        # And remember to delete all the temporary files
        for f in fns:
            os.remove(f)    
    
    return fn_list    
# ----------------------------------------------------------------------------------------

def run_plot_flux_ratio(fn_list, params, suffix=None, 
                        ratios = ['[NII]/[OIII]',],
                        vrange = [[None,None],],
                        use_egal_dered = False,
                       ):   
    ''' 
    This function is designed to make some images from emission lines flux ratios.
    
    :Args:
        fn_list: dictionary
                 The dictionary containing all filenames created by brutifus.
        params: dictionary
                The dictionary containing all paramaters set by the user. 
        suffix: string [default: None]
                The tag of this step, to be used in all files generated for rapid id.
        ratios: list of string [default: ['[NII]/[OIII]']]
             A list of line ratios to plot.
        vrange: list of lists [default: [[Mone,None]]]
                The plot range.
        use_egal_dered: bool [default: False]
                        If available, whether to use the line fluxes corrected for 
                        extragalactic attentuation.
             
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
    ''' 
    
    raise Exception('Check and update this function!!!')
    
       # First, I need to generate 3 individual fits files with the line fluxes
    if params['verbose']:
        print('-> Creating some line ratio images.')
    
    # Very well, now I need to open the line fluxes.
    if use_egal_dered:
        fn = os.path.join(params['prod_loc'],fn_list['dered_elines_params'])
    else:
        fn = os.path.join(params['prod_loc'],fn_list['elines_params_cube'])

    # Alright, now I need to extract the line fluxes
    # Open the file
    hdu = fits.open(fn)
    header0 = hdu[0].header
    header1 = hdu[1].header

    # All elines definition are inside elines_pyqz inside brutifus_metadata. 
    # For the ones the user requested, construct a flux map.
    elines_fluxes = {}
    for (m, mix) in enumerate(ratios):
        elines_fluxes[m] = {0:None, 1:None}    

    # Look at each emission line fitted. Do I need it ?
    for (k,akey) in enumerate(np.sort(params['elines'].keys())):
    
        this_header = hdu[k+1].header
        this_data = hdu[k+1].data
                               
        lam = this_header['BR_REFL']
        
        for (m,mix) in enumerate(ratios):
            
            for (cc,col) in enumerate(mix.split('/')): # Search all possible lines
                if lam in elines_pyqz[col]: # Wavelength match ? I found a line.
            
                    flux = this_data[0]
                
                    # Take into account the SNR
                    flux[this_data[-1]<params['elines'][akey][0][3]] = np.nan
            
                    # Already exists ? - then I found a line with multiple components
                    if not(elines_fluxes[m][cc] is None): 
                        elines_fluxes[m][cc] += flux # append the flux
                    else:
                        elines_fluxes[m][cc] = flux # Save the data for later
                       
    # Close the hdus
    hdu.close()
    
    # Very well, I now have all my line fluxes. Let's deal with each set one at a time.
    for (m,mix) in enumerate(ratios):
        
        ofn = os.path.join(params['plot_loc'],
                 suffix+'_'+params['target']+'_ratio_%s_vs_%s.pdf' % (mix.split('/')[0],
                                                                      mix.split('/')[1]))
                 
        # Now, let' create a fits file with WCS coordinates
        fn = os.path.join(params['prod_loc'],
                 suffix+'_'+params['target']+'_ratio_%s_vs_%s.fits' % (mix.split('/')[0],
                                                                       mix.split('/')[1]))
            
        hdu0 = fits.PrimaryHDU(header = header0)
        hdu1 = fits.ImageHDU(np.log10(elines_fluxes[m][0]/elines_fluxes[m][1]))
        # Make sure the WCS coordinates are included as well
        hdu1 = brutifus_tools.hdu_add_wcs(hdu1,header1)
        # Also include a brief mention about which version of brutifus is being used
        hdu1 = brutifus_tools.hdu_add_brutifus(hdu1,suffix)
        # Add the line reference wavelength for future references
        hdu1.header['BR_RATIO'] = (mix, 'line ratio')
            
        hdu = fits.HDUList(hdus=[hdu0,hdu1])
        hdu.writeto(fn, overwrite=True)
        
        # Great, I am now ready to call the plotting function
        bifus_p.make_2Dplot(fn, ext = 1, ofn = ofn, contours = False, 
                                 stretch = 'linear',
                                 vmin = vrange[m][0], 
                                 vmax = vrange[m][1],
                                 cmap = alligator,
                                 cblabel = r'$\log$ %s' %mix,
                                 cbticks = None,
                                 scalebar = params['scalebar'],
                                )
                                
    
    return fn_list  