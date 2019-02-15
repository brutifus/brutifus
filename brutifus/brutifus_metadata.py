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

__version__ = '2019.02.2'

# Where are we located ?
bifus_dir = os.path.dirname(__file__) # Get the project directory from the file location!

# And get all the different useful locations associated
#refdata_dir = os.path.join(brutifus_dir,'..','reference_data')

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


      