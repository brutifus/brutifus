# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018-2019,  F.P.A. Vogt
 
-----------------------------------------------------------------------------------------
 
This file contains tools for the brutifus routines to create pretty plots seamlessly.

Created November 2018, F.P.A. Vogt - - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

import os
import sys
import pickle

import numpy as np
import scipy.interpolate as interp
from scipy.ndimage.filters import gaussian_filter

import aplpy

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec 
from mpl_toolkits.axes_grid1.anchored_artists import AnchoredSizeBar

from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import astropy.units as u
from astropy.wcs import WCS
import astropy.visualization as astrovis

from . import brutifus_metadata as bifus_m

# Set the proper plotting style
plt.style.use(bifus_m.plotstyle)

# ----------------------------------------------------------------------------------------
def crosshair(inner_r=1, pa = 0):
   ''' The path of an emtpy cross, useful for indicating targets without crowding the field.
   
   :params inner_r: empty inner radius. Default =1 (gap = tick length = 1/3 marker width
   :type inner_r: float
   :params pa: position angle (degrees)   
   :type pa: float
   
   :return: path
   
   '''
   
   verts = [(-1.5, 0),
            (-0.5*inner_r, 0),
            (0, 0.5*inner_r),
            (0, 1.5),
            (0.5*inner_r, 0),
            (1.5, 0),
            (0, -0.5*inner_r),
            (0, -1.5),
            (-1.5, 0),
            (-1.5, 0),
           ]
   
   pa = np.radians(pa)
   rot_mat = np.matrix([[np.cos(pa),-np.sin(pa)],[np.sin(pa),np.cos(pa)]])
   
   for (v, vert) in enumerate(verts):
      verts[v] = (vert*rot_mat).A[0] 
      
   codes = [mpl.path.Path.MOVETO,
            mpl.path.Path.LINETO,
            mpl.path.Path.MOVETO,
            mpl.path.Path.LINETO,
            mpl.path.Path.MOVETO,
            mpl.path.Path.LINETO,
            mpl.path.Path.MOVETO,
            mpl.path.Path.LINETO,
            mpl.path.Path.MOVETO,
            mpl.path.Path.CLOSEPOLY
            ]

   path = mpl.path.Path(verts, codes)
    
   return path


# ----------------------------------------------------------------------------------------
def get_fig_dims(nx, ny):
   ''' Returns all the necessary figure dimensions for gridspec, given the size of the 
       image to plot.
   
   :param nx: number of pixels along x axis
   :type nx: int
   :param ny: number of pixels along y axis
   :type ny: int
   
   :return: list containing the figure parameters [width, height, height_ratios, 
            width_ratios, left, right, botom, top, wspace, hspace]
   :rtype: list
   
   '''
   
    
   fig_height = bifus_m.margin_top + bifus_m.cb_thick + bifus_m.margin_bottom + \
                ny/nx * (bifus_m.fig_width - bifus_m.margin_right - bifus_m.margin_left)
   
   height_ratios = [bifus_m.cb_thick/fig_height, 1]
   
   left = bifus_m.margin_left/bifus_m.fig_width
   right = 1-bifus_m.margin_right/bifus_m.fig_width
   bottom = bifus_m.margin_bottom/fig_height
   top = 1-bifus_m.margin_top/fig_height
   
   wspace = 0.02
   hspace = 0.02
   
   return [np.round(bifus_m.fig_width,3), np.round(fig_height,3), 
           np.round(height_ratios,3), np.round([1],3), 
           np.round(left,3), np.round(right,3), 
           np.round(bottom,3), np.round(top,3), np.round(wspace,3), np.round(hspace,3)]

# ----------------------------------------------------------------------------------------
def get_im_stretch(stretch):
   ''' Returns a stretch to feed the ImageNormalize routine from Astropy.
   
   :param stretch: short, catchy name for the stretch I want. Possibilities are
                   'arcsinh' or 'linear'.
   :type stretch: string
   
   :return: A :class:`astropy.visualization.stretch` thingy ...
   
   '''
   
   if stretch == 'arcsinh':
      return astrovis.AsinhStretch()
   elif stretch == 'linear':
      return astrovis.LinearStretch()
   else:
      raise Exception('Ouch! Stretch %s unknown.' % (stretch))
      
# ----------------------------------------------------------------------------------------    
def get_im_interval(pmin = 10, pmax = 99.9, vmin = None, vmax = None):
   ''' Returns an interval, to feed the ImageNormalize routine from Astropy.
   
   :param pmin: lower-limit percentile
   :type pmin: float
   :param pmax: upper-limit percentile
   :type pmax: float
   :param vmin: absolute lower limit
   :type vmin: float
   :param vmax: absolute upper limit
   :type vmax: float
   
   :return: an :class:`astropy.visualization.interval` thingy ...
   
   .. note:: Specifying *both* vmin and vmax will override pmin and pmax.
   
   '''
   
   if not(vmin is None) and not (vmax is None):
      
      return astrovis.ManualInterval(vmin,vmax)
   
   else:
      
      return astrovis.AsymmetricPercentileInterval(pmin,pmax)
      
# ----------------------------------------------------------------------------------------    
def finetune_WCSAxes(im_ax):
   ''' Fine-tune the look of WCSAxes plots to my liking.
   
   Currently, this changes the tick label format to show integer seconds, sets the 
   separators to hms/d ' ", enables minor ticks, and tweaks the axis labels to be R.A. and 
   Dec. (with dots!). 
   
   :param im_ax: the axes to be improved
   :type im_ax: :class:`matplotlib.axes`
   
   :return: (ra,dec), the axes coordinates.
   
   ''' 
   
   im_ax.set_facecolor((0.5,0.5,0.5))
   
   # Deal with the axes
   ra = im_ax.coords[0]
   dec = im_ax.coords[1]
   
   ra.set_major_formatter('hh:mm:ss')
   dec.set_major_formatter('dd:mm:ss')
   
   # This only works when I am using proper LaTeX
   ra.set_separator((r'$^\mathrm{h}$', r'$^\mathrm{m}$', r'$^\mathrm{s}$ '))
   dec.set_separator((r'$^{\circ}$', r'$^{\prime}$', r'$^{\prime\prime}$'))
   
   ra.display_minor_ticks(True)
   dec.display_minor_ticks(True)
   
   im_ax.set_xlabel('R.A. [J2000]')
   im_ax.set_ylabel('Dec. [J2000]', labelpad=0)
   
   return (ra,dec)

# ----------------------------------------------------------------------------------------  
# Here, I define the "alligator" colormap to give brutifus a unique feel, just because I 
# can. More importantly, I can decide to set the lower and upper bounds to B&W, so that 
# spaxels outside the colorbar range are directly identifiable. Also, I can make sure that 
# nan's show up in 50% grey. 

# Alligator: 6 nodes with linear evolution, suitable for B&W impression
cdict_alligator = {
'red'  :  ( (0.00, 0./255, 0./255),
            (0.00, 0./255, 0./255),      (0.2, 20./255., 20./255.),
            (0.40, 46./255., 46./255.),  (0.6, 108./255., 108./255.),
            (0.8, 207./255., 207./255.), (1.00, 255./255.,255./255.),
            (1.00, 255./255, 255./255),),

'green' : ( (0.00, 0./255, 0./255),
            (0.00, 25./255., 25./255.),  (0.2, 85./255., 85./255.),
            (0.4, 139./255., 139./255.), (0.6, 177./255., 177./255.),  
            (0.8, 234./255, 234./255),   (1.00, 248./255.,248./255.),
            (1.00, 255./255, 255./255),),
            
'blue':   ( (0.00, 0./255, 0./255),
            (0.00, 25./255., 25./255.),  (0.2, 81./255., 81./255.),
            (0.4, 87./255., 87./255.),   (0.6, 86./255., 86./255.),
            (0.8, 45./255, 45./255),     (1.00, 215./255.,215./255.),
            (1.00, 255./255, 255./255),),   
}

# Initiate the colorbar
alligator = plt.matplotlib.colors.LinearSegmentedColormap('alligator',cdict_alligator, 1024)
# Set the bad colors
alligator.set_bad(color=(0.5,0.5,0.5), alpha=1) 

# ----------------------------------------------------------------------------------------
def reverse_colourmap(cdict):
   ''' Returns the inverted colorbar dictionary.

   I most definitely got this from the web somewhere ... Stackoverflow?
      
   :params cdict: a colorbar dictionary
   :type cdict: dict
      
   :return: the flipped dictionary
   :rtype: dict
      
   '''

   new_cdict = {}
   for channel in cdict:
      data = []
      for t in cdict[channel]:
         data.append((1. - t[0], t[1], t[2]))
        
      new_cdict[channel] = data[::-1]
   return new_cdict

# ----------------------------------------------------------------------------------------
# Define the inverse alligator cmap as well
cdict_alligator_r = reverse_colourmap(cdict_alligator)
alligator_r = plt.matplotlib.colors.LinearSegmentedColormap('alligator_r', 
                                                            cdict_alligator_r, 1024) 
# Set the bad colors
alligator_r.set_bad(color=(0.5,0.5,0.5), alpha=1) 

# ----------------------------------------------------------------------------------------      
def make_galred_plot(lams, alams, etau, Ab, Av, ofn):
   ''' Plots the extinction Alambda and flux correction factor for galactic extinction.
    
   :param lams: The wavelength nodes.
   :type lams: 1D numpy array
   :param alams: The corresponding values of Alambda.
   :type lams: 1D numpy array
   :param etau: The flux correction factor, i.e. F_(unextinct)/F_(observed).
   :type lams: 1D numpy array
   :param Ab: The value of Ab.
   :type lams: float
   :param Av: The value of Av.
   :type lams: float
   :param ofn: path+name of where to save the plot
   :type ofn: string 
   
   '''
   
   # Start the plotting
   plt.close(1)
   plt.figure(1,figsize=(6.93,6))
   gs = gridspec.GridSpec(1,1, height_ratios=[1], width_ratios=[1])
   gs.update(left=0.14, right=0.85, bottom=0.12, top=0.93, wspace=0.02, hspace=0.03)
  
   ax = plt.subplot(gs[0,0])
   
   ax.plot(lams,alams,'k-',drawstyle='steps-mid')
   
   ax.set_xlabel(r'Wavelength [\AA]')
   ax.set_ylabel(r'A$_\lambda$ [mag]')
   
   ax2 = ax.twinx()
   ax2.plot(lams,etau,'--',c='firebrick',drawstyle='steps-mid')
   ax2.set_ylabel(r'F$_{\lambda}$/F$_{\lambda,obs}= e^{\tau}$', color='firebrick',
                  labelpad=10)
   
   # make the other y-axis crimson for clarity
   ax2.tick_params(axis='y',colors='firebrick', which='both')
   ax2.spines['right'].set_color('firebrick')
   
   # Add the legend
   ax.text(0.8,0.9,'A$_{B}$=%.3f\n A$_{V}$=%.3f' % (Ab,Av), transform=ax.transAxes,
           verticalalignment='center',
           horizontalalignment='center',
         )
   
   # Tweak the limits
   ax.set_xlim((np.min(lams),np.max(lams)))
   
   plt.savefig(ofn)
   plt.close() 
   
# ----------------------------------------------------------------------------------------         
def show_scale(ax, scale_length = 10.*u.arcsec, scale_text = None, scale_loc = 'lower left'):
   ''' Adds a scalebar to a given 2D plot. 
   
   :param ax: The axes to which to add the sale.
   :param scale_length: The scale length in astropy.units
   :param scale_text: the string to print as the scalebar, or None
   :param scale_loc: string [default: 'top left']
                 The location of the scale bar. e.g 'top right', 'bottom left', etc ...
                 
   '''
   
   raise Exception('Ouch! The scalebar is not working yet ...!')
   
   # Make a default scale text in case none is specified
   if scale_text is None:
      
      the_unit = scale_length.unit.to_string()
      
      if the_unit == 'arcsec':
         the_unit = r'$^{\prime\prime}$'
      elif the_unit == 'arcmin':
         the_unit = r'$^{\prime}$'
      else:
         the_unit = ' ' + the_unit
                
      scale_text = '%.1f%s' % (scale_length.value, the_unit)
   
   # TODO:
   # Looks like I still need to divide the length by the "pixel scale" ? 
   # Unless I need to divide by the Dec ?
   # But how do I extract this from the ax ???
      
   # Create the scalebar
   bar = AnchoredSizeBar(ax.get_transform('world'), scale_length.to(u.degree).value, 
                         scale_text, scale_loc, 
                         sep = 5, pad = 0.5, borderpad = 0.4)
    
   # Adjust the bar thickness
   bar.size_bar.get_children()[0].set_linewidth(2)
   # Add it to the plot
   ax.add_artist(bar)
   
# ----------------------------------------------------------------------------------------      

def make_2Dplot(fn, # path to the data (complete!)
                ext = 0, # Which extension am I looking for ?
                ofn = 'plot.pdf', # Savefig filename
                stretch = 'arcsinh',
                vmin = None, 
                vmax = None,
                pmin = 10.0,
                pmax = 99.99,
                gauss_blur = None,
                cmap = None,
                cblabel = None,
                scalebar = None,
                ):
   ''' Creates an image from a 2-D fits image with WCS info.

   :param fn: The filename (+path!) fo the fits file to display.
   :type fn: string
   :param ext: The extension of the image to display. The data in this extension MUST be 2-D.
   :type ext: int
   :param ofn: The output filename (+path!)
   :type ofn: string
   :param stretch: The stretch of the image, fed to aplpy.FITSFigure. 'linear', 'arcsinh', ...
   :type stretch: string
   :param vmin: If set, the lower bound of the colorbar
   :type vmin: float
   :param vmax: If set, the upper bound of the colorbar
   :type vmax: float
   :param pmin: If set, the lower percentile bound of the colorbar
   :type pmin: float
   :param pmax: If set, the upper percentile bound of the colorbar
   :type pmax: float
   :param cmap: If set, the colormap to use. Use 'alligator' for the special brutifus cmap.
   :type cmap: string
   :param cblabel: If set, the label of the colorbar.
   :type cblabel: string
   :param scalebar: If set, adds a scale bar to the plot. Format: [lenght arcsec, length kpc, loc, unit]
   :type scalebar: list of len(3)
            
   :return: a list containing the figure, the ax1 and the plot filename.
   :rtype: list
        
   .. note:: This function absolutely requires WCS coordinates, and a 2-D array.
      
   '''
   
   # Alright, let's start by computing the proper figure height, so that my plots all
   # have the same ratio, irrespective of whether I have a colorbar, or not.
   
   # First get the data
   
   hdu = fits.open(fn)
   data = hdu[ext].data
   header = hdu[ext].header
   hdu.close()
   
   # What is the size of my image ?
   (ny,nx) = np.shape(data)
   
   # What is the associated fig dimensions I need ?
   fig_dims = get_fig_dims(nx, ny)
   
   # If requested, smooth the array
   if not(gauss_blur is None):
      data = gaussian_filter(data, sigma = gauss_blur)
      
   # Let's create a plot, and forget about using aplpy
   plt.close(1)
   fig = plt.figure(1, figsize = (fig_dims[0], fig_dims[1]))
   
   # Set the scene with gridspec
   gs = gridspec.GridSpec(2,1, 
                          height_ratios = fig_dims[2],
                          width_ratios = fig_dims[3], 
                          left = fig_dims[4], 
                          right = fig_dims[5], 
                          bottom = fig_dims[6], 
                          top = fig_dims[7], 
                          wspace = fig_dims[8], 
                          hspace = fig_dims[9])
   
   ax1 = plt.subplot(gs[1,0], projection=WCS(header))
   
   # What is the colorscheme selected ?
   if cmap == 'alligator':
      mycmap = alligator
   elif cmap == 'alligator_r':
      mycmap = alligator_r
   elif cmap is None:
      mycmap = 'Greys'
   else:
      mycmap = cmap
   
   # What is the normalization ?
   norm = astrovis.ImageNormalize(data, 
                                  interval = get_im_interval(pmin = pmin, pmax = pmax,
                                                             vmin = vmin, vmax = vmax),
                                  stretch = get_im_stretch(stretch),
                                  clip = False)
   
   # Plot the image
   im = ax1.imshow(data, origin = 'lower', interpolation = 'nearest', cmap = mycmap,
                   norm = norm)
   
   #ax1.set_nan_color((0.5,0.5,0.5))
   # This does not work !!!!
   # TODO:
   ax1.set_facecolor('salmon')
   
   # Add the scalebar
   #if scalebar is not None:
   #   show_scale(ax1, scale_length = scalebar[0], scale_text = scalebar[1], 
   #              scale_loc = scalebar[2])
   
   # Make the axes look pretty
   (ra,dec) = finetune_WCSAxes(ax1)
   
   # Deal with the colorbar
   if not(cmap is None):
      
      ax0 = plt.subplot(gs[0,0])
   
      cb = plt.colorbar(cax = ax0, mappable=im, orientation='horizontal', 
                        ticklocation = 'top')
      cb.set_label(cblabel, labelpad = 10)
      #cb.ax.xaxis.set_ticks_position('top')
      ax0.tick_params(axis='x', which='major', pad = 5)
      
      # BUG: can't get the minor ticks working properly - turn them off for now
      cb.ax.minorticks_off()
   
   fig.savefig(ofn)
   
   return (fig, ax1, ofn)
   
# ----------------------------------------------------------------------------------------  
def make_RGBplot(fns, ofn, 
                 ext = [0,0,0],
                 stretch = ['linear','linear','linear'], 
                 plims = [0.25, 99.75, 0.25,99.75, 0.25, 99.75], 
                 vlims = [None, None, None, None, None, None],
                 gauss_blur = [None,None,None],
                 title = None,
                 scalebar = None,
                ):
   ''' Creates an RGB image from three fits files.
   
   :param fns: The filename (+path!) fo the  3 fits file to display (in R, G and B orders).
   :type fns: list
   :param ofn: The filneame (+path) of the output file.
   :type ofn: str
   :param ext: What HDU extension to read ?
   :type ext: list
   :param stretch: The stretch to apply to the data, e.g. 'linear', 'log', 'arcsinh'.
   :type stretch: list
   :param plims: The limiting percentiles for the plot, as 
                 [pmin_r, pmax_r, pmin_g, pmax_g, pmin_b, pmax_b]
   :type plims: list
   :param vlims: The limting values for the plot (superseeds plims), as
                 [vmin_r, vmax_r, vmin_g, vmax_g, vmin_b, vmax_b]
   :type vlims: list
   :param gauss_blur: list of radius (in pixels) for gaussian_blur. None = no blur
   :type gauss_blur: list
   :param title: Title to display above the image
   :type title: str
   :param scalebar: If set, adds a scale bar to the plot. 
                    Format: [lenght arcsec, length kpc, loc, unit]   
   :type scalebar: list 

   :return: a list containing the figure, the ax1 and the plot filename.
   :rtype: list
   
   '''  
    
   # If I was given a single stretch for all images, enlarge it
   if type(stretch) == np.str:
      stretch = [stretch,stretch,stretch]
   
   # If I was given a single extension, assume it is the same everywhere
   if type(ext) == np.int:
      ext = [ext,ext,ext]
   
   # Open the data and headers
   data = [] 
   header = []
   
   for (f,fn) in enumerate(fns):
      hdu = fits.open(fn)
      
      # Get the data
      this_data = hdu[ext[f]].data
      
      # If requested, smooth the array
      if not(gauss_blur[f] is None):

         this_data = gaussian_filter(this_data, sigma = gauss_blur[f])
      
      # Normalize it
      this_data = get_im_interval(pmin = plims[2*f], pmax = plims[2*f+1], 
                                  vmin = vlims[2*f], vmax = vlims[2*f+1])(this_data)
      
      # Stretch it
      this_data = get_im_stretch(stretch[f])(this_data)
      
      data += [this_data]
      header += [hdu[ext[f]].header]
   
    # What is the size of my images ?
   (ny,nx) = np.shape(data[0])
   
   # What is the associated fig dimensions I need ?
   fig_dims = get_fig_dims(nx, ny)
   
   # Let's create a plot, and forget about using aplpy
   plt.close(1)
   fig = plt.figure(1, figsize = (fig_dims[0], fig_dims[1]))
   
   # Set the scene with gridspec
   gs = gridspec.GridSpec(2,1, 
                          height_ratios = fig_dims[2],
                          width_ratios = fig_dims[3], 
                          left = fig_dims[4], 
                          right = fig_dims[5], 
                          bottom = fig_dims[6], 
                          top = fig_dims[7], 
                          wspace = fig_dims[8], 
                          hspace = fig_dims[9])
   
   ax1 = plt.subplot(gs[1,0], projection = WCS(header[0]))
   
   im = ax1.imshow(np.transpose(np.array(data), (1,2,0)), 
                   origin = 'lower', 
                   interpolation = 'nearest', zorder = 0)
   
   
   # Assume an astro image, and set the ticks to white
   ax1.tick_params(axis='both', which='both', color='w', direction='in')
   
   # Make the axes look pretty
   (ra,dec) = finetune_WCSAxes(ax1)
   
   
   if not(title is None):
      ax1.set_title(title)

   fig.savefig(ofn)
        
   return (fig, ax1, ofn)
