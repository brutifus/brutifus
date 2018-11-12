# -*- coding: utf-8 -*-
'''
brutus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2016-2018,  F.P.A. Vogt

-----------------------------------------------------------------------------------------

This file contains some global metadata used throughout the brutus code to fit IFU data 
cubes. Includes reference wavelengths, the speed of light, etc ...

Created April 2016, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
Updated November 2018, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

import numpy as np
import os

import matplotlib.pyplot as plt

# ---| Some basic parameters |------------------------------------------------------------

__version__ = '2018.11.0'

# Where are we located ?
bifus_dir = os.path.dirname(__file__) # Get the project directory from the file location!

# And get all the different useful locations associated
#refdata_dir = os.path.join(brutifus_dir,'..','reference_data')

# Specify the plotting style
usetex = True
plotstyle = os.path.join(bifus_dir,'mpl_styles','brutifus_plots.mplstyle')

# ---| Information about the file format for the different instruments |------------------

ffmt = {'MUSE': {'data':1, 'var':2, 'badpix':None},
        #'WiFeS':{'data':0, 'error':1, 'badpix':2},
       }

# ---| Some constants |-------------------------------------------------------------------

# Some emission line rest-wavelengths (in air and Ansgtroem)
'''
o2l = 3726.09
o2h = 3728.83
hd = 4101.73
hc = 4340.47
hb = 4861.32
o3l = 4958.91
o3h = 5006.73
he1 = 5875.66
nal = 5889.95
nah = 5895.92
o1 = 6300.304
n2l = 6548.04
ha = 6562.81
n2h = 6583.46
s2l = 6716.44
s2h = 6730.81
'''

# ---| Information about the different stellar models for ppxf|---------------------------
'''
sl_models ={'MILES_ppxf_default':{'sl_loc':os.path.join(refdata_dir,
                                                        'stellar_templates',
                                                        'miles_models'),
                                  'nAges':26, 'nMetal':6, 
                                  'metal':[-1.71, -1.31, -0.71, -0.40, 0.00, 0.22],
                                  'logAge':np.linspace(np.log10(1),np.log10(17.7828),26),
                                  'metal_str':['m1.71','m1.31','m0.71','m0.40','p0.00',
                                               'p0.22'],
                                  'fwhm':2.51,
                                }, 
            'MILES_Padova00_kb_1.30':{'sl_loc':os.path.join(refdata_dir,
                                                            'stellar_templates',
                                                            'MILES_Padova00_kb_1.30'),
                                  'nAges':50, 'nMetal':7, 
                                  'metal':[-2.32, -1.71, -1.31, -0.71, -0.40, 0.00, 0.22],
                                  'logAge':np.linspace(np.log10(0.0631),
                                                       np.log10(17.7828),
                                                       50),
                                  'metal_str':['m2.32', 'm1.71', 'm1.31', 'm0.71', 
                                               'm0.40', 'p0.00', 'p0.22'],
                                  'fwhm':2.51,
                                },
            'MILES_Padova00_kb_1.30_sub':{'sl_loc':os.path.join(refdata_dir,
                                                            'stellar_templates',
                                                            'MILES_Padova00_kb_1.30_sub'),
                                  'nAges':25, 'nMetal':6, 
                                  'metal':[-1.71, -1.31, -0.71, -0.40, 0.00, 0.22],
                                  'logAge':np.linspace(np.log10(0.0631),
                                                       np.log10(15.8489),
                                                       25),
                                  'metal_str':['m1.71', 'm1.31', 'm0.71', 
                                               'm0.40', 'p0.00', 'p0.22'],
                                  'fwhm':2.51,
                                },  
            'MILES_Padova00_kb_1.30_sub_sub':{'sl_loc':os.path.join(refdata_dir,
                                                            'stellar_templates',
                                                            'MILES_Padova00_kb_1.30_sub_sub'),
                                  'nAges':13, 'nMetal':5, 
                                  'metal':[-1.31, -0.71, -0.40, 0.00, 0.22],
                                  'logAge':np.linspace(np.log10(0.0631),
                                                       np.log10(15.8489),
                                                       13),
                                  'metal_str':['m1.31', 'm0.71', 'm0.40', 'p0.00', 'p0.22'],
                                  'fwhm':2.51,
                                },                   
           }

'''         




      