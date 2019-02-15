# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018-2019  F.P.A. Vogt

-----------------------------------------------------------------------------------------
 
'''

# ---| Import the packages |--------------------------------------------------------------

import sys
import os
import pickle
import datetime
import argparse

from matplotlib import pyplot as plt

import brutifus as bifus
import brutifus.brutifus_metadata as bifus_m


# ---| Toggle on-off all the different processing steps to be run |-----------------------
# This list effectively defines the "recipes" to be run by brutifus. In practice, it should 
# be flexible, and could/should be adaptable to anyone's specific needs.

proc_steps = [
   
   # Correct the WCS values to check with Gaia.
   {'step':'adjust_WCS', 'run':True, 'suffix':'00',
    'args':{'name_in':'raw_cube',
            'name_out':'wcs-corr_cube', 
            'pmin': 10.0, # pmin for the plot
            'pmax': 99.0, # pmax for the plot
           },
   },
   # First, compute some S/N maps for the continuum and emission lines of interest ...
   {'step':'crude_snr_maps', 'run':False, 'suffix':'01', 
    'args':{'name_in':'raw_cube', # What cube do I want to process ?
            'do_plot':True,  # Do I want to save some pretty pictures ?
            'zcorr_lams': False, # True = correct wavelength ranges for redshift
           },
   },
   # Make some RGB plots from the datacube ...
   {'step':'plot_RGB', 'run':False, 'suffix':'02', 
    'args':{'name_in':'raw_cube', # What cube do I want to plot ?
            'stretches':['arcsinh','arcsinh','arcsinh'],
            'stretch_plims': [[1.,98.,1.,98.,1.,98.],],
            'bands':[[[6560,6570],[7500,8750],[4750,5500]],],
            },
   },
   # Perform some sky subtraction, if needed
   {'step':'sky_sub', 'run':True, 'suffix':'03', 
    'args':{'name_in':'raw_cube',
            'name_out':'skysub_cube'},
   },
   # Correct for galactic reddening via the values given on NED
   {'step':'gal_dered', 'run':False, 'suffix':'04',
    'args':{'name_in':'skysub_cube',
            'name_out':'galdered_cube',
            'do_plot':True},
   },
   
   # ----- Continuum fitting ------------------------------------------------------------
   # First, fit the continuum using a LOWESS filter. Much more elegant than a polynomial!
   {'step':'fit_continuum', 'run':False, 'suffix':'05',
    'args':{'name_in':'galdered_cube', 
            'start_row': 0, # Where to start the fitting ? None = 0
            'end_row': None, # Where to end the fitting ? None = max
            'method': 'lowess', # What kind of fitting is this ?
           },
   },
   # Gather the output in a single datacube
   {'step':'make_continuum_cube', 'run':False, 'suffix':'06',
    'args':{'method':'lowess',
           },
   },
   # Make some RGB plots ...
   {'step':'plot_RGB', 'run':False, 'suffix':'07', 
    'args':{'name_in':'lowess_cube', # What cube do I want to plot ?
            'stretches':['linear','linear','linear'],
            'stretch_plims': [[10.,90,10.,90.,10.,90.],],
            'bands':[[[6532,6562],[6532,6562],[6532,6562]],],
            'conts':[[[6470,6500],[6470,6500],[6470,6500]],],
            
            },
   },
   # Subtract the sky continuum
   {'step':'subtract_continuum', 'run':False, 'suffix':'08',
    'args':{'name_in': 'galdered_cube',
            'name_out':'consub_cube',
            'method':'lowess',
           },
   },
   # Make some RGB plots ...
   {'step':'plot_RGB', 'run':False, 'suffix':'09', 
    'args':{'name_in':'consub_cube', # What cube do I want to plot ?
            'stretches':['linear','linear','linear'],
            'stretch_plims': [[10.,99.,10.,99.,10.,99.],],
            'bands':[[[6532,6562],[7500,8750],[5284,5324]],],
            'conts':[[[6470,6500],None,[5350,5390]],],
            },
   },
   # Make some BW plots ...
   {'step':'plot_BW', 'run':True, 'suffix':'11', 
    'args':{'name_in':'galdered_cube', # What cube do I want to plot ?
            'stretches':['linear','linear', 'linear', 'linear'],
            'stretch_plims': [[0.5,99.9],
                             ],
            'stretch_vlims': [[None]*2,
                             ],
            'gauss_blurs':[None,
                           ],
            'bands':[[6556,6579],
                    ],
            'conts':[[6356,6379],
                     ],
            },
   },
   
   # ----- Emission line fitting --------------------------------------------------------
   # TBD ...
   
   # ----- Fit inspection ---------------------------------------------------------------
   # Interactive window displaying the fit, and comparison between the different 
   # continuum removal.
   #{'step':'inspect_fit','run':False, 'suffix':'10',
   # 'args':{'irange':[1,1e4], 
   #         'vrange':[7200,7400],
   #        },
   #},
   
   # ----- Structures and apertures business --------------------------------------------
   # TBD ...
   
   # ----- Emission line analysis -------------------------------------------------------
   # Compute the electron density from the [SII] line ratio
   #{'step':'get_ne', 'run':False, 'suffix':'15',
   # 'args':{'do_plot':True,
   #         'ratio_range':[0.4,1.4],
   #        },
   #},
   ]
    

# ---| Set up the scene, get the params, etc ... |----------------------------------------
start_time = datetime.datetime.now()

# --- | Use argparse to make brutifus user friendly |-------------------------------------
parser = argparse.ArgumentParser(description='''High-level recipe handling the processing
                                                of IFU datacubes. ''',
                                 epilog =' Full documentation: %s \n \n \
                                           Feedback, questions, comments: \
                                           frederic.vogt@alumni.anu.edu.au \n' % ('http://fpavogt.github.io/brutifus'))

parser.add_argument('--version', action='version', version=('brutifus %s'% bifus.__version__))

parser.add_argument('--params', action='store', metavar='string', 
                    default = None,
                    help='Pickled parameter filename (generated with brutifus_params.py)')

parser.add_argument('--no-systemtex', action='store_true',
                    help='disable the use of the system-wide LaTeX')

# Get all the supplied arguments
args = parser.parse_args()

# Check that I was fed some parameters
if args.params is None:
   raise Exception('Failed to load the pickled parameter file. Have you supplied any?')

# Disable the use of system-Latex if required ... (Why would anyone ask such a thing?!)
if args.no_systemtex:
   bifus_m.usetex = False
   bifus_m.plotstyle = os.path.join(bifus_m.bifus_dir,'mpl_styles',
                                         'brutifus_plots_nolatex.mplstyle')
else:
   bifus_m.usetex = True
   bifus_m.plotstyle = os.path.join(bifus_m.bifus_dir,'mpl_styles',
                                         'brutifus_plots.mplstyle')
                                         
# Set the chosen plot style
plt.style.use(bifus_m.plotstyle)

# Load the parameter file
f = open(args.params, 'rb')
bp = pickle.load(f)
f.close()

# Is there a dictionary of filenames already in place ? If not, create one
fn = os.path.join(bp['prod_loc'],bp['fn_list_fn'])

if not(os.path.isfile(fn)):
   if bp['verbose']:
      print('-> This looks like a fresh run. Creating an emtpy dictionary of filenames.')
   
   # For now, all I know is where the raw data is located
   fn_list = {'raw_cube': os.path.join(bp['data_loc'],bp['data_fn']),
             }
      
   f = open(fn,'wb')
   pickle.dump(fn_list,f) # For now, just an empty dictionary!
   f.close()

# ---| Execute the recipe, by calling all the individual master step functions | ---------
prev_suffix = None
for step in proc_steps:
   step_name   = step['step']
   step_run    = step['run']
   step_suffix = step['suffix']
   step_args   = step['args']
   func_name = 'run_'+step_name
   # Call a function based on a string ... pretty sweet !
   func = getattr(bifus,func_name) 
   
   if step_run:
      # Here, I want to maintain a dictionary of filenames, to be used accross functions
      # For each step, load the dictionary, feed it to the function, and save it update
      # it when it's all done. Each function returns the updated dictionary !
      fn = os.path.join(bp['prod_loc'],bp['fn_list_fn'])
      f = open(fn,'rb')
      fn_list = pickle.load(f)
      f.close()
      
      # Launch the function
      fn_list = func(fn_list, bp, suffix = step_suffix, **step_args)
      
      # Save the updated dictionary of filenames
      f = open(fn,'wb')
      fn_list = pickle.dump(fn_list,f)
      f.close()

# ---| All done ! |-----------------------------------------------------------------------
duration = datetime.datetime.now() - start_time
print('All done in %.01f seconds.' % duration.total_seconds())
 
# ---| End of the World as we know it |---------------------------------------------------
