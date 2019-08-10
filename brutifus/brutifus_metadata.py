# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018-2019,  F.P.A. Vogt

-----------------------------------------------------------------------------------------

This file contains some global metadata used throughout the brutifus code.

Created November 2018-2019, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

import numpy as np
import os

import matplotlib.pyplot as plt

# ---| Some basic parameters |------------------------------------------------------------

# Where are we located ?
bifus_dir = os.path.dirname(__file__) # Get the project directory from the file location!

# Name of the parameters files
bifus_params = 'params_brutifus.yaml'
bifus_procsteps = 'procsteps_brutifus.yaml'

# Name of default storage spaces
plot_loc = 'brutifus_plots'
prod_loc = 'brutifus_products'
tmp_loc =  'brutifus_tmp'


def get_fn_list_fn(target):
   ''' Returns the filename of the storage pickle for all the filenames used by the code
      
      :param target: name of the target/object to be processed
      :type target: string 
   
   
      :return: the pickled dictionary filename
      :rtype: string
   '''
      
   return 'bifus_fn-list_%s.pkl' % (target) # Name of the dictionary for filenames

# ---| Plotting parameters |--------------------------------------------------------------

usetex = True
plotstyle = os.path.join(bifus_dir,'mpl_styles','brutifus_plots.mplstyle')

fig_width = 6.94 # In inches
margin_left = 1.3 # in inches
margin_right = 0.2 # in inches
margin_bottom = 0.75 # in inches
margin_top = 0.7 # in inches
cb_thick =  0.3 # in inches

# ---| Information about the file format for the different instruments |------------------

ffmt = {'MUSE': {'data':1, 'var':2, 'badpix':None, 
                 'funit':r'F$_\lambda$ [$10^{-20}$ erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]'},
        'WiFeS':{'data':0, 'error':1, 'badpix':2, 
                 'funit':r'F$_\lambda$ [erg s$^{-1}$ cm$^{-2}$ \AA$^{-1}$]'},
       }


      