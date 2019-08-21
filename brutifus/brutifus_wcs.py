# -*- coding: utf-8 -*-
# pylint: disable=C0103,W0311

'''
 brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
 Copyright (C) 2018-2019,  F.P.A. Vogt

----------------------------------------------------------------------------------------------------

 This file contains WCS-related function for the higher-leve brutifus routines.

 Created August 2019, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# --------------------------------------------------------------------------------------------------

import warnings
import os
import numpy as np
from scipy import signal

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt

from astropy.io import fits
import astropy.units as u
from astropy.coordinates.sky_coordinate import SkyCoord
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from astropy.wcs import WCS

from astroquery.gaia import Gaia

from photutils import DAOStarFinder

from . import brutifus_metadata as bifus_m
from . import brutifus_plots as bifus_p

# --------------------------------------------------------------------------------------------------
def get_linear_WCS_corr(fn, obstime=None, verbose=True, suffix=None, target=None):
   ''' Identify the linear offset between Gaia and the WCS solution of a given 2D image.

   :param fn: 2D FITS image filename (and location)
   :type fn: string
   :param obstime: Observing time and date
   :type obstime: astropy.Time()
   :param suffix: the tag of this step (for plot filename)
   :type suffix: str
   :param target: target name (for plot filename)
   :type target: str
   :param verbose: print some info to screen
   :type verbose: bool

   :return: (dx,dy), the linear offsets in pixels
   :rtype: tuple

   .. note:: Assumes the image is in the primary extension.

   '''

   # Get the data
   hdu = fits.open(fn)
   header = hdu[0].header
   wl_im = hdu[0].data
   hdu.close()

   # First, let's get some of the statistics on this image
   wl_mean, wl_median, wl_std = sigma_clipped_stats(wl_im, sigma=3.0)

   if verbose:
      print('   Looking for stars in the field ...')

   # Source: https://photutils.readthedocs.io/en/stable/detection.html
   # Get ready to look for stars in there
   daofind = DAOStarFinder(fwhm=3.0, threshold=5. * wl_std)

   # Start the search
   sources = daofind(wl_im - wl_median)

   # Create a list of sources in pixels
   sources_list = np.array([sources['xcentroid'], sources['ycentroid']])

   # Reshape this properly as Nx2
   sources_list = sources_list.T

   # Store the 2D WCS info for later
   wcs2d = WCS(header)

   # Find the coordinate at the center of the image.
   # To be robust (E.g. OCam), should be CRVAl 1/2 + CDELT1/2*(0.5*NAXIS 1/2 - CRPIX1/2)
   # Fixed v. 2019.08.1
   # Credits to J. Suherli for stumbling upon and reporting the bug.
   core_coord = SkyCoord(ra=header['CRVAL1'] +
                            header['CD1_1'] * (0.5*header['NAXIS1'] - header['CRPIX1']),
                         dec=header['CRVAL2'] +
                             header['CD2_2'] * (0.5*header['NAXIS2'] - header['CRPIX2']),
                         unit=(u.deg, u.deg), frame='icrs')
   search_radius = 0.75 * max([np.abs(header['NAXIS1'] * header['CD1_1']*u.deg),
                               np.abs(header['NAXIS2'] * header['CD2_2']*u.deg)])

   if verbose:
      print('   Querying GAIA ...')

   # Make it a sync search, because I don't think I need more than 2000 targets ...
   j = Gaia.cone_search_async(core_coord, search_radius, verbose=False)
   r = j.get_results()

   if len(r) >= 2000:
      warnings.warn('   Maximum number of GAIA entries returned. Some stars might be missing.')

   # Some of these proper motions are NaN's ... need to deal with this before feeding them
   # to SkyCoord
   r['pmra'].fill_value = 0
   r['pmdec'].fill_value = 0

   # Turn this into a "catalogue" of SkyCoord
   gaia_cat = SkyCoord(ra=r['ra'], dec=r['dec'], unit=(u.deg, u.deg), frame='icrs',
                       pm_ra_cosdec=r['pmra'].filled(),
                       pm_dec=r['pmdec'].filled(),
                       obstime=[Time(item, format='decimalyear') for item in r['ref_epoch']],
                       equinox='J2000',
                       distance=1.e3 * u.pc, # Assumes the default distance
                       )

   # Propagate Gaia to the observation date!
   if not obstime is None:
      if verbose:
         print('   Propagating proper motions to %s ...' % (obstime))
      gaia_cat = gaia_cat.apply_space_motion(new_obstime=obstime)

   # Clean this up, to only show the points within the image
   gaia_cat = gaia_cat[wcs2d.footprint_contains(gaia_cat)]

   # Transform this in pixel coords
   gaia_cat_pix = wcs2d.wcs_world2pix(gaia_cat.ra, gaia_cat.dec, 0)

   # Reshape it properly
   gaia_cat_pix = np.array(gaia_cat_pix).T

   # Ok, let's make two fake image from the source list, that I then cross-correlate
   nx = header['NAXIS1']
   ny = header['NAXIS2']

   im_ref = np.zeros((ny, nx))
   im_obs = np.zeros((ny, nx))

   # Fill the two images with fake stars
   # v2019.08.2: this gets slow with large numbers of stars
   # (thanks Janette Suherli for reporting this)
   # let's be smarter about this ... only create stars on small sub-arrays,
   # and add these to the final image (rather than creating full images each time).
   # Note: the FWHM of the stars in the fake images is fixed in pixels ... should I change this ?

   fwhm_pix = 2

   # Let's get ready to fill the images.

   for (im, entries) in [(im_obs, sources_list), (im_ref, gaia_cat_pix)]:

      # Let's create the fake observed image ...
      for i in range(len(entries)):

         # Pixel ID closest from the star centroid
         xc = np.floor(entries[i][0])
         yc = np.floor(entries[i][1])

         # If we are too close from the edge, skip the star
         # (otherwise I have array size issues)
         # That's the easy way out of course ... !
         if (xc < 2 * fwhm_pix + 3) or (xc > nx - (2 * fwhm_pix + 3)) or \
            (yc < 2 * fwhm_pix + 3) or (yc > ny - (2 * fwhm_pix + 3)):
            continue

         # Find the x and y indices of the sub-array pixels
         x_ind, y_ind = np.meshgrid(np.arange(xc-(2 * fwhm_pix + 1), xc+(2 * fwhm_pix + 1), 1),
                                    np.arange(yc-(2 * fwhm_pix + 1), yc+(2 * fwhm_pix + 1), 1))

      # For each of these, compute the distance to the star center
      d = np.sqrt((x_ind - entries[i][0])**2 + (y_ind - entries[i][1])**2)

      # And now add the star to the full image
      im[y_ind.astype(int), x_ind.astype(int)] += np.exp(-((d)**2 / (2.0 * fwhm_pix**2)))


   # Having done that, let's convolve the two images
   # Credit to Nick Whyborn @ ALMA for the idea
   # Don't forget to flip one of the image to get the correlation
   im_corr = signal.convolve(im_ref, im_obs[::-1, ::-1], mode='same', method='auto')

   # Now, let's find the brightest pixel
   peak_loc = np.unravel_index(im_corr.argmax(), im_corr.shape)

   # Let's also find the "fine" position of the peak with photutils
   mean, median, std = sigma_clipped_stats(im_corr, sigma=10.0)
   daofind = DAOStarFinder(fwhm=5.0, threshold=10.*std)
   peaks = daofind(im_corr)

   # Add a clean error, in case I find no peak.
   if not peaks: # Neat trick ... empty sequences == False !
      raise Exception('Ouch! No WCS correlation peak identified.')

   # Ok, which one is the brightest peak ?
   peak_id = np.argmax(peaks['peak'])
   peak_loc_fine = (peaks['xcentroid'][peak_id], peaks['ycentroid'][peak_id])

   # Issue a warning if the brightest peak is too far from the brightest pixel
   if np.abs(peak_loc_fine[0] - peak_loc[1]) > 2 or \
      np.abs(peak_loc_fine[1] - peak_loc[0]) > 2:

      warnings.warn('WCS correlation peak does not coincide with the brightest pixel!')

   # Ok, so what are the offsets I should apply ?
   dx = peak_loc_fine[0]-np.floor(nx/2) # Added the floor... it is required by 'same' ?
   dy = peak_loc_fine[1]-np.floor(ny/2) # v2019.07.1

   if verbose:
      print('   Best offset: %+.2f  %+.2f pixels' % (dx, dy))

   # Very well, now apply the corrections to the image
   fits.setval(fn, 'CRPIX1', value=header['CRPIX1'] - dx, ext=0)
   fits.setval(fn, 'CRPIX2', value=header['CRPIX2'] - dy, ext=0)

   # Let's make a plot to show how good I'm doing with getting a reasonable offset
   plt.close(99)
   fig99 = plt.figure(99, figsize=(6.94, 7))

   gs = gridspec.GridSpec(1, 1, height_ratios=[1], width_ratios=[1],
                          left=0.12, right=0.98, bottom=0.12, top=0.9,
                          wspace=0.1, hspace=0.1)

   ax99 = plt.subplot(gs[0, 0])

   ax99.imshow(im_corr, cmap='gray', origin='lower', interpolation='nearest')

   ax99.plot(peak_loc_fine[0], peak_loc_fine[1], marker=bifus_p.crosshair(pa=0, inner_r=0.5),
             c='darkorange', markersize=50)

   ax99.set_aspect('equal')
   ax99.set_xlabel('x [pix]')
   ax99.set_ylabel('y [pix]')
   ax99.set_title('Correlation peak: (%+.2f,%+.2f) spaxels' % (dx, dy))

   fig99.savefig(os.path.join(bifus_m.plot_loc, suffix + '_' + target +
                              '_get-WCS-offsets.pdf'))
   plt.close(99)

   # Now plot the white-light image to show how good I did ...
   (fig, ax, ofn) = bifus_p.make_2Dplot(fn, ext=0,
                                        ofn=os.path.join(bifus_m.plot_loc, suffix + '_' + target +
                                                         '_check-WCS.pdf'),
                                        plims=[10, 99],
                                        cmap=None, cblabel=None)

   # And add them to the plot
   ax.scatter(gaia_cat.ra, gaia_cat.dec, marker='o', color='darkorange',
              s=20, facecolor='none', transform=ax.get_transform('world'))

   # Show the stars I found in the data ... these are in pix coordinates ...
   # -> always "right" irrespective of the WCS!
   ax.scatter(sources['xcentroid'], sources['ycentroid'],
              color='crimson', s=50, marker='o', facecolor='none')

   # Save the updated plot
   fig.savefig(ofn)

   return (dx, dy)

# --------------------------------------------------------------------------------------------------
