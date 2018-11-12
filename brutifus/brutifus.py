# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018,  F.P.A. Vogt
 
-----------------------------------------------------------------------------------------
 
This file contains the master brutifus routines to fit the stellar continuum and the 
emission lines in an IFU data cube (i.e. MUSE). Most of these routines call sub-routines,
after setting the scene/loading datasets/etc ...
 
Any processing step MUST have a dediacted routine in this file call 'run_XXX', which can
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

from astropy.io import fits
from astropy.utils.exceptions import AstropyDeprecationWarning
from astropy.io.fits.verify import VerifyWarning

from astropy import log
log.setLevel('WARNING') # Disable pesky INFO flags from aplpy

# For debugging purposes: turn all warnings into errors
#warnings.filterwarnings("error", category=AstropyDeprecationWarning)
#warnings.filterwarnings("error", category=VerifyWarning)


from matplotlib import MatplotlibDeprecationWarning
#warnings.filterwarnings("error",category=MatplotlibDeprecationWarning)
import matplotlib.gridspec as gridspec 
import matplotlib.pyplot as plt
from astropy.visualization import (PercentileInterval, AsymmetricPercentileInterval,
                                   AsinhStretch, LinearStretch, ImageNormalize)

# Import brutifus-specific tools
from . import brutifus_tools
from . import brutifus_cof
from . import brutifus_elf
from . import brutifus_plots as bifus_p
from . import brutifus_red
from .brutifus_metadata import *
from .brutifus_metadata import __version__

# Import other important tools, but don't fail, in case the user doesn't need them.
try:
    import brutifus_ppxf
    import ppxf_util as util
except:
    print(" ... failed to load brutifus_ppxf ...")

try:
    import fit_kinematic_pa as fkp
except:
   print(" ... failed to load fit_kinematic_pa ...")


# ---------------------------------------------------------------------------------------- 
  
def run_snr_maps(fn_list, params, suffix = None, do_plot = False):
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
      do_plot: bool [default: True]
               Whether to make a plot of the SNR maps or not.
        
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

   # Import the raw fits file 
   hdu = fits.open(os.path.join(params['data_loc'],params['data_fn']))
   header0 = hdu[0].header
   data = hdu[ffmt[params['inst']]['data']].data
   header1 = hdu[ffmt[params['inst']]['data']].header
   error = hdu[ffmt[params['inst']]['var']].data
   header2 = hdu[ffmt[params['inst']]['var']].header
   hdu.close()
    
   # Build the wavelength array - REST frame !
   lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
   lams /= params['z_target']+ 1.
   
   # Prepare a storage list 
   snrs=[]
   
   # Here, hide those stupid warnings about all-NaNs slices
   with warnings.catch_warnings():
      warnings.simplefilter("ignore", category=RuntimeWarning)
   
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
   hdu = brutifus_tools.hdu_add_wcs(hdu,header1)
   # Also include a brief mention about which version of brutifus is being used
   hdu = brutifus_tools.hdu_add_brutifus(hdu,suffix)
   # For reference, also include the line/region this maps are based on
   hdu.header['B_SNRANG'] = ('%.1f-%.1f' % (lams[0],lams[-1]), 'spectral range (A) used for SNR')
   hdu.header['B_SNTYPE'] = ('x', 'binary map, 0 = no data, 1 = valid spectra')
   
   hdus += [hdu]
   
   # Now also loop through all the maps I have created
   for (i,r) in enumerate(params['snr_ranges']):
   
      hdu = fits.ImageHDU(snrs[i])
      
      # Make sure the WCS coordinates are included as well
      hdu = brutifus_tools.hdu_add_wcs(hdu,header1)
      
      # Also include a brief mention about which version of brutifus is being used
      hdu = brutifus_tools.hdu_add_brutifus(hdu,suffix)
      
      # For reference, also include the line/region this maps are based on
      hdu.header['B_SNRANG'] = ('%.1f-%.1f' % (r[0],r[1]), 'spectral range (A) used for SNR')
      hdu.header['B_SNTYPE'] = ('%s' % r[-1], '"c"ontinuum, or "e"mission')
      
      hdus += [hdu]
        
   # Write the file!
   hdu = fits.HDUList(hdus=hdus)
   fn_out = os.path.join(params['prod_loc'],suffix+'_'+params['target']+'_snr.fits')
   hdu.writeto(fn_out, overwrite = True)
    
   # And add the filename to the dictionary of filenames
   fn_list['snr_cube'] = suffix+'_'+params['target']+'_snr.fits'
        
   if do_plot:
      # Alright, let's take out the big guns ...    
      for (i,r) in enumerate(params['snr_ranges']):
             
         bifus_p.make_2Dplot(fn_out,i+2, 
                             os.path.join(params['plot_loc'],
                                          suffix + '_' + params['target'] + 
                                          '_snr_%.1f-%.1f.pdf' % (r[0],r[1])),
                             contour_fn = fn_out, 
                             contour_ext = i+2,
                             contour_levels = [5], 
                             vmin=0, vmax=50,
                             cbticks=[0, 3, 5, 10, 50], 
                             cblabel = r"S/N ('%s') %.1f\AA-%.1f\AA" % (r[-1],r[0],r[1]),
                             scalebar = params['scalebar'],
                            )  
                                            
   return fn_list
# ---------------------------------------------------------------------------------------- 
def run_plot_RGB(fn_list, params, suffix=None, 
                        bands = [[[7500.,9300.],[6000.,7500.],[4800.,6000.]],], 
                        conts = [[None, None, None],],
                        stretches = ['log',],
                        stretch_plims = [[10.,99.5,10.,99.5,10.,99.5],],
                        stretch_vlims = [[None,None,None,None,None,None],],
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
                     The limtiing values for the plot (superseeds stretch_plims), as
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
    
   # Very well, now I need to open the datacube.
   fn = os.path.join(params['data_loc'],params['data_fn'])

   # Alright, now I need to extract the line fluxes
   # Open the file
   hdu = fits.open(fn)
   header0 = hdu[0].header
   header1 = hdu[ffmt[params['inst']]['data']].header
   data = hdu[ffmt[params['inst']]['data']].data
   hdu.close()
    
   # Build the wavelength array in Angstroem
   lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']

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
         hdu = brutifus_tools.hdu_add_wcs(hdu, header1)
         outfits = fits.HDUList([hdu])
         outfits.writeto(fn,overwrite=True)
            
        
      ofn = os.path.join(params['plot_loc'],
               suffix+'_'+params['target']+'_RGB_%i-%i_%i-%i_%i-%i.pdf' % (band[0][0],
                                                                           band[0][1],
                                                                           band[1][0],
                                                                           band[1][1],
                                                                           band[2][0],
                                                                           band[2][1]))
                                                                             
      # Great, I am now ready to call the plotting function
      bifus_p.make_RGBplot(fns, ofn,  stretch = stretches[i],
                           plims = stretch_plims[i], vlims = stretch_vlims[i],
                           title = r'\smaller R: %s-%s\,\AA | G: %s-%s\,\AA | B: %s-%s\,\AA' %  
                                   (band[0][0],band[0][1],
                                    band[1][0], band[1][1],
                                    band[2][0], band[2][1]),
                           scalebar = params['scalebar'])
    
      # And remember to delete all the temporary files
      for f in fns:
         os.remove(f)    
    
   return fn_list    
# ----------------------------------------------------------------------------------------

def run_sky_sub(fn_list, params, suffix = None):
   '''Performs  a manual sky subtraction, when the observation strategy was flawed.
    
   :Args:
      fn_list: dictionary
               The dictionary containing all filenames created by brutifus.
      params: dictionary
              The dictionary containing all paramaters set by the user. 
      suffix: string [default: None]
              The tag of this step, to be used in all files generated for rapid id.
        
   :Returns:
      fn_list: dictionary
               The updated dictionary of filenames. 
                 
   '''
   
   if params['verbose']:
      print('-> Subtracting the sky ...')
    
   # Very well, now I need to open the datacube.
   fn = os.path.join(params['data_loc'],params['data_fn'])

   # Open the file
   hdu = fits.open(fn)
   header0 = hdu[0].header
   header1 = hdu[ffmt[params['inst']]['data']].header
   data = hdu[ffmt[params['inst']]['data']].data
   header2 = hdu[ffmt[params['inst']]['var']].header
   error = hdu[ffmt[params['inst']]['var']].data
   hdu.close()
    
   # Build the wavelength array in Angstroem
   lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
   
   # Assemble a 2D map of my sky spaxels
   sky_spaxels = np.zeros_like(data[0,:,:])
   cys,cxs = np.mgrid[0:header1['NAXIS2'],0:header1['NAXIS1']]
    
   for sr in params['sky_regions']:
         
      if len(sr) == 3:
         
         d = np.sqrt( (cxs-sr[0])**2 + (cys-sr[1])**2)
         sky_spaxels[d<=sr[2]] = 1
      
      elif len(sr) ==4:
         
         sky_spaxels[sr[1]:sr[1]+sr[3]+1, sr[0]:sr[0]+sr[2]+1] = 1
   
   # Very well, assemble the sky spectrum now
   sky_spec = np.array([np.nanmedian(plane[sky_spaxels==1]) for plane in data[:]])
   
   # Make a descent plot of the sky spectrum
   plt.close(1)
   plt.figure(1, figsize=(14.16,4))
   gs = gridspec.GridSpec(1,1, height_ratios=[1], width_ratios=[1],
                          left=0.08, right=0.98, bottom=0.18, top=0.98, wspace=0.02, hspace=0.03)
   
   ax1 = plt.subplot(gs[0,0])
   
   ax1.semilogy(lams,sky_spec, 'k-', 
            drawstyle='steps-mid', lw=0.8)
   
   ax1.set_xlim((lams[0],lams[-1]))
   ax1.set_xlabel('Observed wavelength [\AA]')
   ax1.set_ylabel(r'F$_\lambda$ [$10^{-20}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]')
   
   plt.savefig(os.path.join(params['plot_loc'], suffix+'_'+params['target']+'_skyspec.pdf'))
   plt.close(1)
   
   # In addition, also make a plot showing all the sky locations
   plt.figure(1, figsize=(6.94,6.5))
   gs = gridspec.GridSpec(1,1, height_ratios=[1], width_ratios=[1],
                          left=0.08, right=0.98, bottom=0.1, top=0.98, wspace=0.02, hspace=0.03)
   
   ax1 = plt.subplot(gs[0,0])
   
   # Build a white-light image
   norm = ImageNormalize(np.sum(data,axis=0), interval=AsymmetricPercentileInterval(0,90),
                         stretch=AsinhStretch())
   ax1.imshow(np.sum(data, axis = 0), origin = 'lower', cmap = 'gray_r', norm=norm)
   
   # Show all the sky spaxels
   sky_spaxels[sky_spaxels ==0] = np.nan
   ax1.imshow(sky_spaxels, origin= 'lower', cmap='Reds', vmin = 0, vmax = 1, alpha=0.7)
   
   # Fine tune things a bit
   ax1.set_xlabel('x [pixel]')
   ax1.set_ylabel('y [pixel]')
   
   plt.savefig(os.path.join(params['plot_loc'], suffix+'_'+params['target']+'_skyspaxels.pdf'))
   plt.close(1)
   
   # Turn the spectrum into a cube 
   sky_cube = np.repeat(np.repeat(sky_spec[:,np.newaxis],header1['NAXIS2'],axis=1)[:,:,np.newaxis],header1['NAXIS1'],axis=2)
   
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
      hdu = brutifus_tools.hdu_add_wcs(hdu,header1)
      hdu = brutifus_tools.hdu_add_lams(hdu,header1)
      # Also include a brief mention about which version of brutifus is being used
      hdu = brutifus_tools.hdu_add_brutifus(hdu,suffix)
      # Add the line reference wavelength for future references
                     
   hdu = fits.HDUList(hdus=[hdu0,hdu1,hdu2])
   fn_out = os.path.join(params['prod_loc'],
                         suffix+'_'+params['target']+'_skysub.fits')
   hdu.writeto(fn_out, overwrite=True)
        
   # Add the filename to the dictionary of filenames
   fn_list['skysub_cube'] = suffix+'_'+params['target']+'_skysub.fits'
   
   return fn_list
# ----------------------------------------------------------------------------------------
  
def run_gal_dered(fn_list, params, suffix = None, do_plot = False):
    '''Corrects for Galactic extinction, given the Ab and Av extinction.
    
    This function erives the Alambda value for any wavelength, and corrects the data to
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
        
    :Returns:
        fn_list: dictionary
                 The updated dictionary of filenames. 
                 
    :Notes:
        To reproduce the approach from NED, use the Ab and Av value for you object from 
        there, and set curve='f99', rv=3.1.   
    '''
    
    raise Exception('Check and update this function!!!')
    
    
    if params['verbose']:
        print('-> Correcting for Galactic extinction.')
    
    # Import the raw data 
    hdu = fits.open(os.path.join(params['data_loc'],params['data_fn']))
    header0 = hdu[0].header
    data = hdu[1].data
    header1 = hdu[1].header
    error = hdu[2].data
    header2 = hdu[2].header
    hdu.close()
    
    # Build the wavelength array
    lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    
    # Compute Alambda
    alams = brutifus_red.alam(lams, params['Ab'],params['Av'], curve=params['gal_curve'],
                           rv=params['gal_rv'])
    
    # Compute the flux correction factor
    etau = brutifus_red.galactic_red(lams,params['Ab'],params['Av'], 
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
        hdu = brutifus_tools.hdu_add_wcs(hdu,header1)
        hdu = brutifus_tools.hdu_add_lams(hdu,header1)
        # Also include a brief mention about which version of brutifus is being used
        hdu = brutifus_tools.hdu_add_brutifus(hdu,suffix)
        # Add the line reference wavelength for future references
        hdu.header['BR_AB'] = (params['Ab'], 'Extinction in B')
        hdu.header['BR_AV'] = (params['Av'], 'Extinction in V')
            
            
    hdu = fits.HDUList(hdus=[hdu0,hdu1,hdu2])
    fn_out = os.path.join(params['prod_loc'],
                              suffix+'_'+params['target']+'_galdered.fits')
    hdu.writeto(fn_out, overwrite=True)
        
    # Add the filename to the dictionary of filenames
    fn_list['galdered_cube'] = suffix+'_'+params['target']+'_galdered.fits'
    
    # Do I want to save a plot ?
    if do_plot:
        ofn = os.path.join(params['plot_loc'],suffix+'_'+params['target']+
                            '_gal_Alambda_corr.pdf')
        bifus_p.make_galred_plot(lams,alams,etau,params['Ab'],params['Av'], ofn)
    
    return fn_list
    
# ----------------------------------------------------------------------------------------

def run_fit_continuum(fn_list, params, suffix=None, start_row = None, end_row = None, 
                      method='lowess'):
   ''' 
   This function fits the continuum in the datacube, either using ppxf (if SNR is decent)
   or using a simple polynomial value. It is designed to use multiprocessing to speed
   things up on good computers. It deals with the data columns-per-columns, and 
   can be restarted mid-course, in case of a crash. 
   '''
    
   # For the bad continuum: run a LOWESS filter from statsmodel.nonparametric
   #sm.nonparametric.smoothers_lowess.lowess(spec,lams,frac=0.05, it=5)
   # For the good fits, run ppxf
    
   # Rather than launch it all at once, let's be smart in case of problems. I'll run
   # the fits row-by-row with multiprocessing (hence up to 300cpus can help!), and save
   # between each row.
   
   # First, load the datacube to be fitted. Use the galactic deredened one, if it exists.
   #if os.path.isfile(os.path.join(params['prod_loc'],fn_list['galdered_cube'])):
   #   hdu = fits.open(os.path.join(params['prod_loc'],fn_list['galdered_cube']))
   #else: 
   #hdu = fits.open(os.path.join(params['data_loc'],params['data_fn']))
   hdu = fits.open(os.path.join(params['prod_loc'],fn_list['skysub_cube']))
   
   header0 = hdu[0].header
   data = hdu[ffmt[params['inst']]['data']].data
   header1 = hdu[ffmt[params['inst']]['data']].header
   error = hdu[ffmt[params['inst']]['var']].data
   header2 = hdu[ffmt[params['inst']]['var']].header
   hdu.close()
    
   # Build the wavelength array
   lams = np.arange(0, header1['NAXIS3'],1) * header1['CD3_3'] + header1['CRVAL3']
    
   # I also need to load the SNR cube for the spaxel selection
   hdu = fits.open(os.path.join(params['prod_loc'],fn_list['snr_cube']))
   snr_cont = hdu[1+params['lowess_snr_id']].data
   hdu.close()

   # Get some info about the cube
   nrows = header1['NAXIS1']
   if start_row is None:
      start_row = 0
   if end_row is None:
      end_row = nrows-1

   # Ok, what do I want to do ?
   if method == 'lowess':
      if params['verbose']:
         print('-> Starting the continuum fitting using the LOWESS approach.')
            
      fit_func = partial(brutifus_cof.lowess_fit, lams=lams, 
                         frac=params['lowess_frac'], it=params['lowess_it'])
      # Note here the clever use of the partial function, that turns the lowess_fit
      # function from something that takes 4 arguments into something that only takes 1
      # argument ... thus perfect for the upcoming "map" functions !
        
   elif method == 'ppxf':
      if params['verbose']:
         print('-> Starting the continuum fitting using PPXF.')
      
      raise Exception('Check me!')
      # I need to do some special preparation for 'ppxf'. Namely, log-rebin the 
      # templates. I could do that for each fit, bit it would cost me time. At the price 
      # of a larger code here, do it only once, and save "tons" of time (~0.5s per fit).
      # The code below is based on the examples provided by M. Cappellari within ppxf 
      # istelf (ppxf_kinematics_example_sauron.py & ppxf_population_example_sdss.py), 
      # but has been modified to account a wavelength dependant spectral resolution of 
      # the data. And avoid doing the same thing multiple times via multiprocessing. 
      
      # The full spectral range
      lam_range = np.array([lams[0],lams[-1]])
      
      # De-redshift it using the user's best guess
      # Here's the original text from M. Cappellari:
      #
      # << If the galaxy is at a significant redshift (z > 0.03), one would need to 
      # apply a large velocity shift in PPXF to match the template to the galaxy 
      # spectrum. This would require a large initial value for the velocity (V>1e4 km/s)
      # in the input parameter START = [V,sig]. This can cause PPXF to stop!
      # The solution consists of bringing the galaxy spectrum roughly to the
      # rest-frame wavelength, before calling PPXF. In practice there is no
      # need to modify the spectrum before the usual LOG_REBIN, given that a
      # red shift corresponds to a linear shift of the log-rebinned spectrum.
      # One just needs to compute the wavelength range in the rest-frame
      # and adjust the instrumental resolution of the galaxy observations. >>
      lams0 = lams/(params['z_target']+1)
      lam_range0 = np.array([lams0[0],lams0[-1]])
        
      # We can only fit the spectra where they overlap with the spectral library.
      # Get one of the templates, figure out its range, and derive the fit limits.
      sl_fns = glob.glob(os.path.join(sl_models[params['ppxf_sl_name']]['sl_loc'],'*'))
      hdu_sl = fits.open(sl_fns[0])
      header_sl = hdu_sl[0].header
      hdu_sl.close()
      
      lams_sl = np.arange(0, header_sl['NAXIS1'],1)*header_sl['CDELT1'] + \
                                                                       header_sl['CRVAL1']
      lam_range_sl = [lams_sl[0],lams_sl[-1]]
        
      # What are my fit limits, then ?
      fit_lims =[np.max([lam_range0[0],lam_range_sl[0]]), 
                 np.min([lam_range0[1],lam_range_sl[1]])]  
      mask = (lams0 > fit_lims[0]) * (lams0 < fit_lims[1])
      lams0c = lams0[mask]
      lam_range0c = np.array([lams0c[0], lams0c[-1]])
        
      # Soon, I will be log-rebining the spectra. What are the new bins going to be ?
      log_lams0c, emptyspec, velscale = brutifus_tools.log_rebin(lams0c,np.zeros_like(lams0c), 
                                                     sampling = params['ppxf_sampling'])
        
      # Deal with the template spectra: disperse them according to the instrument, and
      # log-rebin them
      if params['verbose']:
         sys.stdout.write('\r   Preparing the templates ... ')
         sys.stdout.flush()
            
      templates, lam_range_temp, logAge_grid, metal_grid, log_lams_temp = \
      brutifus_ppxf.setup_spectral_library(velscale, params['inst'], 
                                          params['ppxf_sl_name'], fit_lims, 
                                          params['z_target'])
      if params['verbose']:   
         sys.stdout.write('done.')
         print(' ')
         sys.stdout.flush()                                      
                                          
                                          
      # For the fit reconstruction later on, save the various wavelength ranges, etc ...
      fn = os.path.join(params['tmp_loc'],
                          suffix+'_'+params['target']+'_ppxf_lams.pkl')
      file = open(fn,'w')
      pickle.dump([lams,lams0,lams0c,log_lams0c],file)
      file.close()
      # And add the generic pickle filename to the dictionary of filenames
      fn_list['ppxf_lams'] = suffix+'_'+params['target']+'_ppxf_lams.pkl'              
                                          
      # From M. Cappellari:
      # << The galaxy and the template spectra do not have the same starting wavelength.
      # For this reason an extra velocity shift DV has to be applied to the template
      # to fit the galaxy spectrum. We remove this artificial shift by using the
      # keyword VSYST in the call to PPXF below, so that all velocities are
      # measured with respect to DV.>>
      #
      dv = c*np.log(lam_range_temp[0]/lam_range0c[0]) 
        
      # Initial estimate of the galaxy velocity in km/s:        
      vel = 0.   # It's been de-redshifted (using user's guess)
      start = [vel, 100.]  # (km/s), starting guess for [V,sigma]
      
      # Now create a list of "good pixels", i.e. mask emission lines and stuff.
      # List is the default one from ppxf.
      goodpixels = util.determine_goodpixels(np.log(log_lams0c), lam_range_temp, vel/c)

      # See the pPXF documentation for the keyword REGUL, and an explanation of the 
      # following two lines
      templates /= np.median(templates) # Normalizes templates by a scalar
      regul_err = params['ppxf_regul_err']  # Desired regularization error
      
      fit_func = partial(brutifus_ppxf.ppxf_population, templates=templates, 
                         velscale=velscale, start=start, goodpixels=goodpixels, 
                         plot=False, moments=params['ppxf_moments'], 
                         degree=params['ppxf_degree'], vsyst=dv, clean=False, 
                         mdegree=params['ppxf_mdegree'], regul=1./regul_err)    

   # Very well, let's start the loop on rows. If the code crashes/is interrupted, you'll
   # loose the current row. Just live with it.
   for row in np.linspace(start_row,end_row, end_row-start_row+1).astype(int):   
      
      # Alright, now deal with the spaxels outside the user-chosen SNR range.
      # Replace them with nan's
      good_spaxels = np.ones((header1['NAXIS2']))
      if params[method+'_snr_min']:
         good_spaxels[snr_cont[:,row] < params[method+'_snr_min']] = np.nan
      if params[method+'_snr_max']:
         good_spaxels[snr_cont[:,row] >= params[method+'_snr_max']] = np.nan

      # Build a list of spectra to be fitted
      specs = [data[:,i,row] * good_spaxels[i] for i in range(header1['NAXIS2'])]
      errs = [error[:,i,row] * good_spaxels[i] for i in range(header1['NAXIS2'])]
        
      if method == 'ppxf':
         raise Exception('Check me too!')
         if params['verbose']:
            sys.stdout.write('\r   Log-rebin spectra in row %2.i ...            '%row)
            sys.stdout.flush()
                
         # I need to do some special preparation for 'ppxf'. Namely, log-rebin the 
         # spectra, and crop them to an appropriate wavelength after de-redshifting. 
         # Let's get started. 
         # To save time, only log-rebin spectra that are not all nans
         specs = [brutifus_tools.log_rebin(lam_range0c, this_spec[mask],
                                     sampling=params['ppxf_sampling'])[1] 
                                     if not(np.all(np.isnan(this_spec))) else np.nan 
                                     for this_spec in specs]                      
         # Also take care of the error
         errs = [brutifus_tools.log_rebin(lam_range0c, this_err[mask], 
                                     sampling=params['ppxf_sampling'])[1]
                                     if not(np.all(np.isnan(this_err))) else np.nan  
                                     for this_err in errs]                                                                 
                                                                            
         # Combine both spectra and errors, to feed only 1 element to fit_func
         # Turn the error into a std to feed ppxf
         specs = [[specs[i],np.sqrt(errs[i])] for i in range(header1['NAXIS2'])]
        
      # Set up the multiprocessing pool of workers
      if params['multiprocessing']:
         # Did the user specify a number of processes to use ?
         if type(params['multiprocessing']) == np.int:
            nproc = params['multiprocessing']
                
         else: # Ok, just use them all ...
            nproc = multiprocessing.cpu_count()
                
         # Limiting the maximum number of cpus for ppxf, to avoid memory errors... 
         if nproc > 8 and method == 'ppxf':
            nproc = 8
                
         pool = multiprocessing.Pool(processes = nproc, 
                                     initializer = brutifus_tools.init_worker())    

         if params['verbose']:
            sys.stdout.write('\r   Fitting spectra in row %2.i, %i at a time ...' % 
                                 (row,nproc))
            sys.stdout.flush()

         # Launch the fitting ! Make sure to deal with KeyBoard Interrupt properly
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

      # Here, I need to save these results. Pickle could be fast and temporary,
      # Until I then re-build the entire cube later on ? Also allow for better 
      # row-by-row flexibility.
      fn = os.path.join(params['tmp_loc'],
                        suffix+'_'+params['target']+'_'+method+'_row_'+
                        str(np.int(row)).zfill(4)+'.pkl')
      file = open(fn,'wb')
      pickle.dump(conts,file)
      file.close()
        
   # And add the generic pickle filename to the dictionary of filenames
   fn_list[method+'_pickle'] = suffix+'_'+params['target']+'_'+method+'_row_'
        
   print(' done !')
    
   return fn_list
# ----------------------------------------------------------------------------------------

def run_make_continuum_cube(fn_list, params, suffix=None, method='lowess', 
                            do_sub_sky=True):   
   ''' 
   This function is designed to construct a "usable and decent" datacube out of the
   mess generated by the continuum fitting function, i.e. out of the many pickle 
   files generated.
   '''
    
   if params['verbose']:
      print('-> Constructing the datacube for the continuum fitting (%s).' % method)
       
   # First, load the original datacube. I need to know how much stuff was fitted.
   #if fn_list['galdered_cube'] is None:
   #     hdu = fits.open(os.path.join(params['data_loc'],params['data_fn']))
   # else:
   hdu = fits.open(os.path.join(params['data_loc'],params['data_fn']))
   
   header0 = hdu[0].header
   data = hdu[ffmt[params['inst']]['data']].data
   header1 = hdu[ffmt[params['inst']]['data']].header
   error = hdu[ffmt[params['inst']]['data']].data
   header2 = hdu[ffmt[params['inst']]['data']].header
   hdu.close()
    
   nrows = header1['NAXIS1']
    
   cont_cube = np.zeros_like(data) * np.nan

   # In the case of ppxf, also extract the stellar moments.
   # Not the most trustworthy, but interesting nonetheless.
   '''
   if params['ppxf_moments']>0:
      ppxf_sol_map = np.zeros((params['ppxf_moments'],header1['NAXIS2'],
                                                        header1['NAXIS1'])) * np.nan
   '''
   # For ppxf, also create a spectrum with the de-logged rebin original spectra. Useful 
   # to characterize the error associated with log-bin-delog-rebin process employed here.
   delog_raw_cube = np.zeros_like(data) * np.nan
                                                        
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
         elif method =='ppxf':
            # Alright, in this case, the pickle file contains the output from ppxf
            
            # Extract the various wavelength ranges, etc ...
            fn = os.path.join(params['tmp_loc'],fn_list['ppxf_lams'])
            f = open(fn,'r')
            [lams,lams0,lams0c,log_lams0c] = pickle.load(f)
            f.close()
                
            # The fit outputs are all on the log_lams0c grid, and I want them back on
            # lams.  
            # First, extract the bestfit spectra. 
            # For test purposes, also extract the galaxy spectra with pp.galaxy
            bestfits = [pp[0] if not(pp is None) else None for pp in conts]
            galaxys = [pp[1] if not(pp is None) else None for pp in conts]
                
            # And the stellar moments
            ppxf_sol_map[:,:,row] = np.array([pp[2] if not(pp is None) else 
                                                  np.zeros(params['ppxf_moments'])*np.nan 
                                                  for pp in conts]).T
                                                  
            # Now, I need storage that spans the original data range
            full_bestfits = [np.zeros_like(lams0) for this_spec in bestfits]
            full_galaxys = [np.zeros_like(lams0) for this_spec in galaxys]
                
            # Now, I want to delog-rebin the spectra back to the original grid
            for (k,cube) in enumerate([bestfits, galaxys]):
               cube = [brutifus_tools.delog_rebin(log_lams0c, this_spec, lams0c, 
                                                    sampling=params['ppxf_sampling'])[1] 
                       if not(this_spec is None) else None for this_spec in cube]
                            
               # Now, I want to un-crop this, by adding nans elsewhere. Make a loop, 
               # but there's probably a more Pythonic way of doing it. 
               # TODO: try again when my brain won't be that fried ...
               mask = (lams0 >=lams0c[0]) * (lams0<=lams0c[-1])
                
               for (f,fit) in enumerate(cube):
                  if not(fit is None):
                     [full_bestfits,full_galaxys][k][f][mask] = fit
                
               # And finally, I need to de-redshift that spectra. No need to touch 
               # the spectra because I did not touch it before.
                
               # Ready to save the continuum cube
               [cont_cube,delog_raw_cube][k][:,:,row] = \
                                         np.array([full_bestfits, full_galaxys][k]).T
                                     
   # Very well, now let's create a fits file to save this as required.
   hdu0 = fits.PrimaryHDU(None,header0)
   hdu1 = fits.ImageHDU(cont_cube)
   # Make sure the WCS coordinates are included as well
   hdu1 = brutifus_tools.hdu_add_wcs(hdu1,header1)
   hdu1 = brutifus_tools.hdu_add_lams(hdu1,header1)
   # Also include a brief mention about which version of brutifus is being used
   hdu1 = brutifus_tools.hdu_add_brutifus(hdu1,suffix)
   hdu = fits.HDUList(hdus=[hdu0,hdu1])
   fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_'+method+'.fits')
   hdu.writeto(fn_out, overwrite=True)
    
   # And add the filename to the dictionary of filenames
   fn_list[method+'_cube'] = suffix+'_'+params['target']+'_'+method+'.fits'
   
   # Very well, do Iw ant to do the sky subtraction ?
   if do_sub_sky:
      
      if params['verbose']:
         print(' ')
         print('-> Subtracting the sky continuum (%s).' % method)
      
      # Since this is a substraction, and I assume no error on the sky, the errors remain unchanged
      data -= cont_cube 
     
      # And save this to a new fits file
      hdu0 = fits.PrimaryHDU(None,header0)
      hdu1 = fits.ImageHDU(data)
      hdu2 = fits.ImageHDU(error)
      
      for hdu in [hdu1,hdu2]:
         # Make sure the WCS coordinates are included as well
         hdu = brutifus_tools.hdu_add_wcs(hdu,header1)
         hdu = brutifus_tools.hdu_add_lams(hdu,header1)
         # Also include a brief mention about which version of brutifus is being used
         hdu = brutifus_tools.hdu_add_brutifus(hdu,suffix)
      
      hdus = fits.HDUList(hdus=[hdu0,hdu1,hdu2])
      fn_out = os.path.join(params['prod_loc'],
                          suffix+'_'+params['target']+'_skysub-'+method+'.fits')
      hdus.writeto(fn_out, overwrite=True)
    
      # And add the filename to the dictionary of filenames
      fn_list[method+'_skysub'] = suffix+'_'+params['target']+'_skysub-'+method+'.fits'
   
   print(' ')
    
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
    
    raise Exception('Check and update this function!!!')
    
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