# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018,  F.P.A. Vogt
 
-----------------------------------------------------------------------------------------
 
This file contains some global variables and other metadata used throughout the
brutifus set of Python modules to process IFU datacubes. 

This is where (if all goes according to plan) users can modify the different
fit parameters to suit their needs. 

Created November 2018, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au

'''
# ----------------------------------------------------------------------------------------

import os
import pickle

from astropy import units as u

# What is the target name (used in filenames, etc...)
target_name = 'G266'

# What is the name of the resulting pickle file going to be ?
pickle_fn = 'bifus-params_%s.pkl' % (target_name)

# Where is everything ?
# Get the project directory from the file location !
proj_dir = os.path.dirname(__file__)

brutifus_params = {

   # ---| General variables |------------------------------------------------------------
   'target': target_name, # The generic name for files created by brutifus
   'z_target': 0., # Best-guess redshift of the target
   'inst': 'MUSE', # Which instrument took the data ?
   'multiprocessing': True, # Specify the number of cores, or set to True for using all those available.
   'verbose': True,

   # ---| Generic plotting parameters |---------------------------------------------------
   # The location of the scalebar. Set to None if you don't want any. 
   # Format: [length*astropy.units, r'text', loc]
   # 'scalebar': [10.*u.arcsec, r'10$^{\prime\prime}$', 'bottom left'],
                                     
   # ---| Location and name of the data file to process |---------------------------------
   'data_loc': os.path.join(proj_dir,'reduc'), # relative path from this file loc !
   'data_fn': '0509-67_muse.fits',

   # ---| Generic locations for plots and products |-------------------------------------
   'plot_loc': os.path.join(proj_dir,'brutifus_plots') ,
   'prod_loc': os.path.join(proj_dir,'brutifus_products'),
   'tmp_loc': os.path.join(proj_dir,'brutifus_products','tmp'),
    
   'fn_list_fn': 'bifus_fn-list_%s.pkl' % (target_name), # Name of the dictionary for filenames
   
   # ---| Constructing SNR maps |--------------------------------------------------------
   # lambda is in *REST* frame !
   'snr_ranges':[[7400., 8500., 'c'], # lam_min, lam_max, 'c'ontinuum, or 'e'mission
                 [6560., 6570., 'e'], # lam_min, lam_max, 'c'ontinuum, or 'e'mission
                ], 
   
   # ---| Manual sky subtraction |--------------------------------------------------------
   # Here we provide the sky positions in pixel values in the PYTHON scheme (starts at 0).
   'sky_regions':[[290, 327, 3], # [x0,y0,radius] (=circle), or [x0,y0,dx,dy] (=rectangle)
                  [182, 319, 3],
                  [216,  61, 3],
                  [306, 205, 2],
                  [ 72, 118, 2],
                  [171,  82, 3],
                  [170,  65, 3],
                  [ 60, 192, 3],
                  [312, 292, 2],
                  [ 66, 183, 3],
                  [ 93, 226, 3],
                  [168, 297, 3],
                  [265, 134, 3],
                 ],
   
   # ---| Correcting for the Galactic extinction |---------------------------------------
   'gal_curve': 'f99', # Set to 'f99' to follow NED.
   'gal_rv': 3.1, # Set to 3.1 to follow NED.
   'Ab': 0.272, # Ab (Landolt) from NED.
   'Av': 0.206, # Av (Landolt) from NED.
   
   # ---| Continuum fitting |------------------------------------------------------------
   # Lowess specific variables
   'lowess_it': 10, # Number of iteration for sigma-clipping to get rid of outliers
   'lowess_frac': 0.15, # % of the array used for deriving each point. 0.05 = sweet spot?
   
   # ---| Emission line fitting | -------------------------------------------------------
   # TBD ...
            
   # ---| Structure detection and aperture extraction |----------------------------------
   # TBD ...
   
}
  
# Export this as a pickle file    
f = open(pickle_fn, 'wb')
pickle.dump(brutifus_params, f)
f.close()  




          
