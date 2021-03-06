# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018-2020,  F.P.A. Vogt
Copyright (C) 2021, F.P.A. Vogt & J. Suherli
All the contributors are listed in AUTHORS.

Distributed under the terms of the GNU General Public License v3.0 or later.

SPDX-License-Identifier: GPL-3.0-or-later

This file contains tools for the brutifus routines to create pretty plots seamlessly.

Created November 2018, F.P.A. Vogt - - frederic.vogt@alumni.anu.edu.au
'''
# --------------------------------------------------------------------------------------------------

import numpy as np
from scipy.ndimage.filters import gaussian_filter

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
import astropy.visualization as astrovis

from . import brutifus_metadata as bifus_m

# Set the proper plotting style
plt.style.use(bifus_m.plotstyle)

# --------------------------------------------------------------------------------------------------
def reverse_colourmap(cdict):
    ''' Returns the inverted colorbar dictionary.

    I most definitely got this from the web somewhere ... Stackoverflow?

    Args:
        cdict (dict): a colorbar dictionary

    Returns:
        dict: the flipped dictionnary

    '''

    new_cdict = {}
    for channel in cdict:
        data = []
        for t in cdict[channel]:
            data.append((1. - t[0], t[1], t[2]))

        new_cdict[channel] = data[::-1]

    return new_cdict

# --------------------------------------------------------------------------------------------------
# Here, I define the "alligator" colormap to give brutifus a unique feel, just because I
# can. More importantly, I can decide to set the lower and upper bounds to B&W, so that
# spaxels outside the colorbar range are directly identifiable. Also, I can make sure that
# nan's show up in 50% grey.

# Alligator: 6 nodes with linear evolution, suitable for B&W impression
cdict_alligator = {
    'red': ((0.00, 0/255, 0/255),
            (0.00, 0/255, 0/255),
            (0.20, 20/255, 20/255),
            (0.40, 46/255, 46/255),
            (0.60, 108/255, 108/255),
            (0.80, 207/255, 207/255),
            (1.00, 255/255, 255/255),
            (1.00, 255/255, 255/255)),

    'green': ((0.00, 0/255, 0/255),
              (0.00, 25/255, 25/255),
              (0.20, 85/255, 85/255),
              (0.40, 139/255, 139/255),
              (0.60, 177/255, 177/255),
              (0.80, 234/255, 234/255),
              (1.00, 248/255, 248/255),
              (1.00, 255/255, 255/255)),

    'blue': ((0.00, 0/255, 0/255),
             (0.00, 25/255, 25/255),
             (0.20, 81/255, 81/255),
             (0.40, 87/255, 87/255),
             (0.60, 86/255, 86/255),
             (0.80, 45/255, 45/255),
             (1.00, 215/255, 215/255),
             (1.00, 255/255, 255/255))
    }

# Initiate the colorbar
alligator = plt.matplotlib.colors.LinearSegmentedColormap('alligator', cdict_alligator, 1024)
# Set the bad colors
alligator.set_bad(color=(0.5, 0.5, 0.5), alpha=1)

# --------------------------------------------------------------------------------------------------
# Define the inverse alligator cmap as well
cdict_alligator_r = reverse_colourmap(cdict_alligator)
alligator_r = plt.matplotlib.colors.LinearSegmentedColormap('alligator_r', cdict_alligator_r, 1024)
# Set the bad colors
alligator_r.set_bad(color=(0.5, 0.5, 0.5), alpha=1)

# --------------------------------------------------------------------------------------------------
def crosshair(inner_r=1, pa=0):
    ''' The path of an emtpy cross, useful for indicating targets without crowding the field.

    Args:
        inner_r (float, optional): empty inner radius. Defaults to 1.
            Note: full gap = tick length = 1/3 marker width.
        pa (float): position angle (degrees). Defaults to 0.

    Returns:
        matplotlib.path.Path: the marker path.
    '''

    verts = [(-1.5, 0),
             (-0.5 * inner_r, 0),
             (0, 0.5 * inner_r),
             (0, 1.5),
             (0.5 * inner_r, 0),
             (1.5, 0),
             (0, -0.5 * inner_r),
             (0, -1.5),
             (-1.5, 0),
             (-1.5, 0),
            ]

    pa = np.radians(pa)
    rot_mat = np.matrix([[np.cos(pa), -np.sin(pa)], [np.sin(pa), np.cos(pa)]])

    for (v, vert) in enumerate(verts):
        verts[v] = (vert * rot_mat).A[0]

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

# --------------------------------------------------------------------------------------------------
def get_fig_dims(nx, ny):
    ''' Returns all the necessary figure dimensions for gridspec, given the size of the
        image to plot.

    Args:
        nx (int): number of pixels along x axis.
        ny (int): number of pixels along y axis.

    Returns:
        list: list containing the figure parameters [width, height, height_ratios,
            width_ratios, left, right, botom, top, wspace, hspace]
    '''

    fig_height = bifus_m.margin_top + bifus_m.cb_thick + bifus_m.margin_bottom + \
                 ny/nx * (bifus_m.fig_width - bifus_m.margin_right - bifus_m.margin_left)

    height_ratios = [bifus_m.cb_thick/fig_height, 1]

    left = bifus_m.margin_left/bifus_m.fig_width
    right = 1 - bifus_m.margin_right/bifus_m.fig_width
    bottom = bifus_m.margin_bottom/fig_height
    top = 1 - bifus_m.margin_top/fig_height

    wspace = 0.02
    hspace = 0.02

    return [np.round(bifus_m.fig_width, 3), np.round(fig_height, 3),
            np.round(height_ratios, 3), np.round([1], 3),
            np.round(left, 3), np.round(right, 3),
            np.round(bottom, 3), np.round(top, 3), np.round(wspace, 3), np.round(hspace, 3)]

# --------------------------------------------------------------------------------------------------
def get_im_stretch(stretch):
    ''' Returns a stretch to feed the ImageNormalize routine from Astropy.

    Args:
        stretch (str): short name for the stretch I want. Possibilities are 'arcsinh' or 'linear'.

    Returns:
        astropy.visualization.stretch: A :class:`astropy.visualization.stretch` thingy ...
    '''

    if stretch == 'arcsinh':
        return astrovis.AsinhStretch()
    if stretch == 'linear':
        return astrovis.LinearStretch()

    raise Exception('Ouch! Stretch %s unknown.' % (stretch))

# --------------------------------------------------------------------------------------------------
def get_im_interval(pmin=10, pmax=99.9, vmin=None, vmax=None):
    ''' Returns an interval, to feed the ImageNormalize routine from Astropy.

    Args:
        pmin (float, optional): lower-limit percentile. Defaults to 10.
        pmax (float, optional): upper-limit percentile. Defaults to 99.9
        vmin (float, optional): absolute lower limit. Defaults to None.
        vmax (float, optional): absolute upper limit. Defaults to None.

    Returns:
        astropy.visualization.interval: a :class:`astropy.visualization.interval` thingy ...

   .. note:: Specifying *both* vmin and vmax will override pmin and pmax.

   '''

    if vmin is not None and vmax is not None:
        return astrovis.ManualInterval(vmin, vmax)

    return astrovis.AsymmetricPercentileInterval(pmin, pmax)

# --------------------------------------------------------------------------------------------------
def finetune_WCSAxes(im_ax):
    ''' Fine-tune the look of WCSAxes plots to my liking.

    Currently, this changes the tick label format to show integer seconds, sets the
    separators to hms/d ' ", enables minor ticks, and tweaks the axis labels to be R.A. and
    Dec. (with dots!).

    Args:
        im_ax (matplotlib.axes): ther axes to be improved.

    Returns:
        matplotlib.ax.coords: (ra, dec), the ax.coords.
    '''

    im_ax.set_facecolor((0.5, 0.5, 0.5))

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

    return (ra, dec)

# --------------------------------------------------------------------------------------------------
def make_galred_plot(lams, alams, etau, Ab, Av, ofn):
    ''' Plots the extinction Alambda and flux correction factor for galactic extinction.

    Args:
        lams (ndarray): The wavelength nodes.
        alams (ndarray): The corresponding values of Alambda.
        etau (ndarray): The flux correction factor, i.e. F_(unextinct)/F_(observed).
        Ab (float): The value of Ab.
        Av (float): The value of Av.
        ofn (str): path+name of where to save the plot
    '''

    # Start the plotting
    plt.close(1)
    plt.figure(1, figsize=(6.93, 6))
    gs = gridspec.GridSpec(1, 1, height_ratios=[1], width_ratios=[1])
    gs.update(left=0.14, right=0.85, bottom=0.12, top=0.93, wspace=0.02, hspace=0.03)

    ax = plt.subplot(gs[0, 0])

    ax.plot(lams, alams, 'k-', drawstyle='steps-mid')

    ax.set_xlabel(r'Wavelength [\AA]')
    ax.set_ylabel(r'A$_\lambda$ [mag]')

    ax2 = ax.twinx()
    ax2.plot(lams, etau, '--', c='firebrick', drawstyle='steps-mid')
    ax2.set_ylabel(r'F$_{\lambda}$/F$_{\lambda,obs}= e^{\tau}$', color='firebrick',
                   labelpad=10)

    # Make the other y-axis crimson for clarity
    ax2.tick_params(axis='y', colors='firebrick', which='both')
    ax2.spines['right'].set_color('firebrick')

    # Add the legend
    ax.text(0.8, 0.9, 'A$_{B}$=%.3f\n A$_{V}$=%.3f' % (Ab, Av), transform=ax.transAxes,
            verticalalignment='center',
            horizontalalignment='center')

    # Tweak the limits
    ax.set_xlim((np.min(lams), np.max(lams)))

    plt.savefig(ofn)
    plt.close()

# --------------------------------------------------------------------------------------------------
def show_scale(ax, scale_length=10.*u.arcsec, scale_text=None, scale_loc='lower left'):
    ''' Adds a scalebar to a given 2D plot.

    Args:
        ax (matplotlib.axes): The figure axes to which to add the scale.
        scale_length (astropy.units, optional): The scale length. Defaults to 10*u.arcsec
        scale_text (str, optional): the string to print as the scalebar. Defaults to None.
        scale_loc (str, optional): The location of the scale bar. e.g 'top right', 'bottom left',
            etc ... Defaults to 'lower left'.
    '''

    # TODO: this function is broken after leaving aplpy behind ...
    ## Make a default scale text in case none is specified
    #if scale_text is None:
    #
    #    the_unit = scale_length.unit.to_string()
    #
    #    if the_unit == 'arcsec':
    #        the_unit = r'$^{\prime\prime}$'
    #    elif the_unit == 'arcmin':
    #        the_unit = r'$^{\prime}$'
    #    else:
    #        the_unit = ' ' + the_unit
    #
    #    scale_text = '%.1f%s' % (scale_length.value, the_unit)
    #
    ## TODO:
    ## Looks like I still need to divide the length by the "pixel scale" ?
    ## Unless I need to divide by the Dec ?
    ## But how do I extract this from the ax ???
    #
    ## Create the scalebar
    #bar = AnchoredSizeBar(ax.get_transform('world'), scale_length.to(u.degree).value,
    #                      scale_text, scale_loc,
    #                      sep = 5, pad = 0.5, borderpad = 0.4)
    #
    ## Adjust the bar thickness
    #bar.size_bar.get_children()[0].set_linewidth(2)
    ## Add it to the plot
    #ax.add_artist(bar)

    raise Exception('Ouch! The scalebar is not working yet ...!')

# --------------------------------------------------------------------------------------------------
def make_2Dplot(fn, # path to the data (complete!)
                ext=0, # Which extension am I looking for ?
                ofn='plot.pdf', # Savefig filename
                stretch='arcsinh',
                vlims=[None, None],
                plims=[10.0, 99.99],
                gauss_blur=None,
                cmap=None,
                cblabel=None,
                scalebar=None,
                ):
    ''' Creates an image from a 2-D fits image with WCS info.

    Args:
        fn (str): The filename (+path!) fo the fits file to display.
        ext (int, optional): The extension of the image to display. The data in this extension MUST
            be 2-D. Defaults to 0.
        ofn (str, optional): The output filename (+path!). Defaults to 'plot.pdf'.
        stretch (str, optional): The stretch of the image, fed to aplpy.FITSFigure. 'linear',
            'arcsinh', ... Defaults to 'arcsinh'.
        vlims (list, optional): the lower and upper absolute bound of the colorbar.
            Defaults to [None, None].
        plims (list, optional): the lower and upper percentile bound of the colorbar.
            Defaults to [10.0, 99.9].
        cmap (str, optional): If set, the colormap to use. Use 'alligator' for the special brutifus
            cmap. Defaults to None.
        cblabel (str, optional): If set, the label of the colorbar. Defaults to None.
        scalebar (list, optional): If set, adds a scale bar to the plot.
            Format: [lenght arcsec, length kpc, loc, unit]. Defaults to None.

    Returns:
        list: a list containing the figure, the ax1 and the plot filename.

    .. note:: This function absolutely requires WCS coordinates, and a 2-D array.

    '''

    # First get the data
    hdu = fits.open(fn)
    data = hdu[ext].data
    header = hdu[ext].header
    hdu.close()

    # What is the size of my image ?
    (ny, nx) = np.shape(data)

    # What is the associated fig dimensions I need ?
    fig_dims = get_fig_dims(nx, ny)

    # If requested, smooth the array
    if gauss_blur is not None:
        data = gaussian_filter(data, sigma=gauss_blur)

    # Let's create a plot, and forget about using aplpy
    plt.close(1)
    fig = plt.figure(1, figsize=(fig_dims[0], fig_dims[1]))

    # Set the scene with gridspec
    gs = gridspec.GridSpec(2, 1, height_ratios=fig_dims[2], width_ratios=fig_dims[3],
                           left=fig_dims[4], right=fig_dims[5], bottom=fig_dims[6], top=fig_dims[7],
                           wspace=fig_dims[8], hspace=fig_dims[9])

    ax1 = plt.subplot(gs[1, 0], projection=WCS(header))

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
                                   interval=get_im_interval(pmin=plims[0], pmax=plims[1],
                                                            vmin=vlims[0], vmax=vlims[1]),
                                   stretch=get_im_stretch(stretch),
                                   clip=False)

    # Plot the image
    im = ax1.imshow(data, origin='lower', interpolation='nearest', cmap=mycmap, norm=norm)

    # TODO: setting the background color of RGB plots does not work
    ax1.set_facecolor('salmon')

    # Add the scalebar
    #if scalebar is not None:
    #   show_scale(ax1, scale_length = scalebar[0], scale_text = scalebar[1],
    #              scale_loc = scalebar[2])

    # Make the axes look pretty
    (_, _) = finetune_WCSAxes(ax1)

    # Deal with the colorbar
    if cmap is not None:

        ax0 = plt.subplot(gs[0, 0])

        cb = plt.colorbar(cax=ax0, mappable=im, orientation='horizontal', ticklocation='top')
        cb.set_label(cblabel, labelpad=10)
        #cb.ax.xaxis.set_ticks_position('top')
        ax0.tick_params(axis='x', which='major', pad=5)

        #TODO: can't get the minor ticks working properly - turn them off for now
        cb.ax.minorticks_off()

    fig.savefig(ofn)

    return (fig, ax1, ofn)

# --------------------------------------------------------------------------------------------------
def make_RGBplot(fns, ofn, ext=[0, 0, 0],
                 stretch=['linear', 'linear', 'linear'],
                 plims=[0.25, 99.75, 0.25, 99.75, 0.25, 99.75],
                 vlims=[None, None, None, None, None, None],
                 gauss_blur=[None, None, None],
                 title=None,
                 scalebar=None,
                ):
    ''' Creates an RGB image from three fits files.

    Args:
        fns (list): The filename (+path!) fo the  3 fits file to display (in R, G and B orders).
        ofn (str): The filneame (+path) of the output file.
        ext (list of int, optional): What HDU extension to read ? Defaults to [0, 0 ,0]
        stretch (list of str, optional): The stretch to apply to the data for each image. E.g.
            'linear', 'log', 'arcsinh'. Defaults to ['linear', 'linear', 'linear'].
        plims (list of float, optional): The limiting percentiles for the plot, as
            [pmin_r, pmax_r, pmin_g, pmax_g, pmin_b, pmax_b]. Defaults to
            [0.25, 99.75, 0.25, 99.75, 0.25, 99.75].
        vlims (list of floats), optional: The limting values for the plot (superseeds plims), as
            [vmin_r, vmax_r, vmin_g, vmax_g, vmin_b, vmax_b]. Defaults to
            [None, None, None, None, None, None].
        gauss_blur (list of float, optional): list of radius (in pixels) for gaussian_blur.
            None = no blur. Defaults to [None, None, None].
        title (str, optional): Title to display above the image. Defaults to None.
        scalebar (list, optional): If set, adds a scale bar to the plot.
            Format: [lenght arcsec, length kpc, loc, unit]. Defaults to None.

    Returns:
        list: a list containing the figure, the ax1 and the plot filename.
    '''

    # If I was given a single stretch for all images, enlarge it
    if isinstance(stretch, str):
        stretch = [stretch, stretch, stretch]

    # If I was given a single extension, assume it is the same everywhere
    if isinstance(ext, int):
        ext = [ext, ext, ext]

    # Open the data and headers
    data = []
    header = []

    for (f, fn) in enumerate(fns):
        hdu = fits.open(fn)

        # Get the data
        this_data = hdu[ext[f]].data

        # If requested, smooth the array
        if gauss_blur[f] is not None:

            this_data = gaussian_filter(this_data, sigma=gauss_blur[f])

        # If the image if full of NaN's ... issue a proper error, 'cause I can't normalize it!
        if np.size(this_data) == np.size(this_data[np.isnan(this_data)]):
            raise Exception('Ouch ... the %s image if full of NaNs!' % (['R', 'G ', 'B'][f]))

        # Normalize it
        this_data = get_im_interval(pmin=plims[2 * f], pmax=plims[2 * f + 1],
                                    vmin=vlims[2 * f], vmax=vlims[2 * f + 1])(this_data)

        # Stretch it
        this_data = get_im_stretch(stretch[f])(this_data)

        data += [this_data]
        header += [hdu[ext[f]].header]

    # What is the size of my images ?
    (ny, nx) = np.shape(data[0])

    # What is the associated fig dimensions I need ?
    fig_dims = get_fig_dims(nx, ny)

    # Let's create a plot, and forget about using aplpy
    plt.close(1)
    fig = plt.figure(1, figsize=(fig_dims[0], fig_dims[1]))

    # Set the scene with gridspec
    gs = gridspec.GridSpec(2, 1, height_ratios=fig_dims[2], width_ratios=fig_dims[3],
                           left=fig_dims[4], right=fig_dims[5], bottom=fig_dims[6], top=fig_dims[7],
                           wspace=fig_dims[8], hspace=fig_dims[9])

    ax1 = plt.subplot(gs[1, 0], projection=WCS(header[0]))

    ax1.imshow(np.transpose(np.array(data), (1, 2, 0)), origin='lower',
                    interpolation='nearest', zorder=0)

    # Assume an astro image, and set the ticks to white
    ax1.tick_params(axis='both', which='both', color='w', direction='in')

    # Make the axes look pretty
    (_, _) = finetune_WCSAxes(ax1)

    if title is not None:
        ax1.set_title(title)

    # Save the figure
    fig.savefig(ofn)

    return (fig, ax1, ofn)
