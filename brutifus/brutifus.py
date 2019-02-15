# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018-2019,  F.P.A. Vogt
 
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
# from . import brutifus_elf
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
   :type suffix: str
   :param name_in: name tag to identify which cube to use to run the routine 
   :type name_in: str
   :param name_out: name tag to identify which cube comes out of the routine
   :type name_out: str
   :param pmin: lower percentile for the plot
   :type pmin: float
   :param pmax: upper percentile for the plot
   :type pmax: float
   :param max_offset: max offset between GAIA and the data ... in PIXELS!!!
   :type max_offset: float 
   
   
   :return: The updated dictionary of filenames.
   :rtype: dict
   
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
def run_crude_snr_maps(fn_list, params, suffix = None, name_in = None, 
                       zcorr_lams = False):
   ''' Computes a crude S/N ratio for a continuum range or emission line.
   
   This function computes the SNR maps for the continuum or emission lines. 
   It also creates a map of spaxels with any signal at all.
   The resulting maps are saved to a fits file with full header and WCS coordinates.
    
   :param fn_list: The dictionary containing all filenames created by brutifus.
   :type fn_list: dict
   :param params:  The dictionary containing all paramaters set by the user. 
   :type params: dict
   :param suffix: The tag of this step, to be used in all files generated for rapid id.
   :type suffix: str
   :param name_in: name tag to identify which cube to use to run the routine 
   :type name_in: str
   :param zcorr_lams: True -> correct the wavelength range for the target redshift.
   :type zcorr_lams: bool
        
   :return: The updated dictionary of filenames.
   :rtype: dict
                 
   .. note:: This function is a "first guess" of the SNR for latter use in the code. A more
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
        
   # Make some plots
   # First, plot the region with any signal at all
   bifus_p.make_2Dplot(fn_list['snr_maps'], 1, 
                       os.path.join(params['plot_loc'],
                                    suffix + '_' + params['target'] + 
                                    '_valid_spectra.pdf'),
                       #contour_fn = None, 
                       #contour_ext = None,
                       #contour_levels = [None], 
                       vmin = 0, vmax = 1.5,
                       cmap = None, cblabel = None, scalebar = params['scalebar'])  
   
   # Alright, let's take out the big guns ...    
   for (i,r) in enumerate(params['snr_ranges']):
             
      bifus_p.make_2Dplot(fn_list['snr_maps'], i+2, 
                          os.path.join(params['plot_loc'],
                                       suffix + '_' + params['target'] + 
                                       '_snr_%.1f-%.1f.pdf' % (r[0],r[1])),
                          vmin = 0, vmax = 50, cmap = 'alligator_r',
                          cblabel = r"S/N ('%s') %.1f\AA-%.1f\AA" % (r[-1],r[0],r[1]),
                          scalebar = params['scalebar'])     
                                                              
   return fn_list

# ---------------------------------------------------------------------------------------- 
def run_plot_BW(fn_list, params, suffix=None, name_in = None,
                bands = [[7500.,9300.],], 
                conts = [None,],
                stretches = ['arcsinh',],
                stretch_plims = [[10.,99.5],],
                stretch_vlims = [[None,None],],
                gauss_blurs = [None,]
               ):   
   ''' Make some B&W images from slices in the cube.
    
   :param fn_list: The dictionary containing all filenames created by brutifus.
   :type fn_list: dict
   :param params:  The dictionary containing all paramaters set by the user. 
   :type params: dict
   :param suffix: The tag of this step, to be used in all files generated for rapid id.
   :type suffix: string
   :param name_in: name tag to identify which cube to use to run the routine 
   :type name_in: string
   :param bands: A list of wavelength pairs.
   :type bands: list
   :param conts: A list of wavelength pairs channels associated continuum regions.     
   :type conts: list    
   :param stretches: The stretches to apply to the data, e.g. 'linear', 'log', 'arcsinh'.
   :type stretches: list
   :param stretch_plims: The limiting percentiles for the plot
   :type stretch_plims: list
   :param stretch_vlims: The limiting values for the plot (superseeds stretch_plims)
   :type stretch_vlims: list
   :param gauss_blurs: Radius of gaussian blur, None for nothing
   :type gauss_blurs: list
             
   :return: The updated dictionary of filenames.
   :rtype: dictionary
   
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
                           #scalebar = params['scalebar']
                           )
    
      # And remember to delete all the temporary files
      os.remove(fn)    
    
   return fn_list    
     
# ---------------------------------------------------------------------------------------- 
def run_plot_RGB(fn_list, params, suffix=None, name_in = None,
                 bands = [[[7500.,9300.],[6000.,7500.],[4800.,6000.]],], 
                 conts = [[None, None, None],],
                 stretches = [['arcsinh','arcsinh','arcsinh'],],
                 stretch_plims = [[10.,99.5,10.,99.5,10.,99.5],],
                 stretch_vlims = [[None,None,None,None,None,None],],
                 gauss_blurs = [[None,None,None,],]
                 ):   
   ''' Make some RGB images from slices in the raw cube.
    
   :param fn_list: The dictionary containing all filenames created by brutifus.
   :type fn_list: dict
   :param params:  The dictionary containing all paramaters set by the user. 
   :type params: dict
   :param suffix: The tag of this step, to be used in all files generated for rapid id.
   :type suffix: string
   :param name_in: name tag to identify which cube to use to run the routine 
   :type name_in: string
   :param bands: A list of 3 wavelength pairs.
   :type bands: list
   :param conts: A list of 3 wavelength pairs: channels associated continuum regions.     
   :type conts: list    
   :param stretches: The stretches to apply to each slice: 'linear', 'log', 'arcsinh', ...
   :type stretches: list
   :param stretch_plims: The limiting percentiles for the plot
   :type stretch_plims: list
   :param stretch_vlims: The limiting values for the plot (superseeds stretch_plims)
   :type stretch_vlims: list
   :param gauss_blurs: Radius of gaussian blur, None for nothing
   :type gauss_blurs: list
             
   :return: The updated dictionary of filenames.
   :rtype: dictionary
   
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
                           # scalebar = params['scalebar']
                           )
    
      # And remember to delete all the temporary files
      for f in fns:
         os.remove(f)    
    
   return fn_list    
   
# ----------------------------------------------------------------------------------------
def run_sky_sub(fn_list, params, suffix = None, name_in = None, name_out = None):
   ''' Performs  a manual sky subtraction, when the observation strategy was flawed.
    
   :param fn_list: The dictionary containing all filenames created by brutifus.
   :type fn_list: dict
   :param params: The dictionary containing all paramaters set by the user. 
   :type params: dict
   :param suffix: The tag of this step, to be used in all files generated for rapid id.
   :type suffix: str
   :param name_in: name tag to identify which cube to use to run the routine
   :type name_in: str
   :param name_out: name tag to identify which cube comes out of the routine
   :type name_out: str
    
   :return: The updated dictionary of filenames.
   :rtype: dict  
             
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
def run_gal_dered(fn_list, params, suffix = None, name_in = None, name_out = None):
   ''' Corrects for Galactic extinction, given the Ab and Av extinction.
   
   This function derives the Alambda value for any wavelength, and corrects the data to
   correct for our "local" extinction. Intended for extragalactic sources.
   
   :param fn_list: The dictionary containing all filenames created by brutifus.
   :type fn_list: dict
   :param params: The dictionary containing all paramaters set by the user. 
   :type params: dict
   :param suffix: The tag of this step, to be used in all files generated for rapid id.
   :type suffix: str
   :param name_in: name tag to identify which cube to use to run the routine
   :type name_in: str
   :param name_out: name tag to identify which cube comes out of the routine
   :type name_out: str

   :return: The updated dictionary of filenames.
   :rtype: dict
                 
   .. note:: To reproduce the approach from NED, use the Ab and Av value for you object
             from there, and set ``curve='f99'``, ``rv=3.1`` in the params file.   
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
        
   # Make a plot of this.
   ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+
                            '_gal_Alambda_corr.pdf')
   bifus_p.make_galred_plot(lams,alams,etau,params['Ab'],params['Av'], ofn)
    
   return fn_list
    
# ----------------------------------------------------------------------------------------
def run_fit_continuum(fn_list, params, suffix=None, name_in = None,
                      start_row = None, end_row = None,
                      method='lowess'):
   ''' Fits the stellar/nebular continuum using a specific method.
   
   This function fits the continuum in the datacube. It is designed to use multiprocessing
   to speed things up on good computers. It deals with the data columns-per-columns, and
   can be restarted mid-course, in case of a crash. 
   
   :param fn_list: The dictionary containing all filenames created by brutifus.
   :type fn_list: dict
   :param params: The dictionary containing all paramaters set by the user. 
   :type params: dict
   :param suffix: The tag of this step, to be used in all files generated for rapid id.
   :type suffix: str
   :param name_in: name tag to identify which cube to use to run the routine
   :type name_in: str
   :param start_row: the row from which to sart the fitting. None for min.
   :type start_row: int
   :param end_row: the row on which to end the fitting (inclusive). None for max.
   :type end_row: int
   :param method: which method to use for the fitting ?
   :type method: str

   :return: The updated dictionary of filenames.
   :rtype: dict
   
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
def run_make_continuum_cube(fn_list, params, suffix = None, method = 'lowess'):   
   ''' Assemble a continuum cube from the individually-pickled fits.
   
   This function is designed to construct a "usable and decent" datacube out of the
   mess generated by the continuum fitting function, i.e. out of the many pickle 
   files generated.
   
   :param fn_list: The dictionary containing all filenames created by brutifus.
   :type fn_list: dict
   :param params: The dictionary containing all paramaters set by the user. 
   :type params: dict
   :param suffix: The tag of this step, to be used in all files generated for rapid id.
   :type suffix: str
   :param method: which fits to gather ?
   :type method: str
   
   :return: The updated dictionary of filenames.
   :rtype: dict
   
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
   ''' Subtracts the sky from a given cube.
   
   :param fn_list: The dictionary containing all filenames created by brutifus.
   :type fn_list: dict
   :param params: The dictionary containing all paramaters set by the user. 
   :type params: dict
   :param suffix: The tag of this step, to be used in all files generated for rapid id.
   :type suffix: str
   :param name_in: name tag to identify which cube to use to run the routine
   :type name_in: str
   :param name_out: name tag to identify which cube comes out of the routine
   :type name_out: str
               
   :return: The updated dictionary of filenames.
   :rtype: dict
   
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