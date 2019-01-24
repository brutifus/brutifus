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

# What is the name of the resulting pickle file going to be ?
pickle_fn = 'bifus-params_0509-67.pkl'

# Where is everything ?
# Get the project directory from the file location !
proj_dir = os.path.dirname(__file__)

brutifus_params = {

   # ---| General variables |------------------------------------------------------------
   'target': '0509-67', # The generic name for files created by brutifus
   'z_target': 0., # Best-guess redshift of the target
   'inst': 'MUSE', # Which instrument took the data ?
   'multiprocessing': True, # Specify the number of cores, or set to True for using all those available.
   'verbose':True,
   #'warnings':'default', # 'ignore' = silence known (?) warnings about all-nan's spaxels
                          # 'default' = default system behavior.

   # ---| Generic plotting parameters |---------------------------------------------------
   # The location of the scalebar. Set to None if you don't want any. 
   # Format: [length*astropy.units, r'text', loc]
   'scalebar': [10.*u.arcsec, r'10$^{\prime\prime}$', 'bottom left'],
                                     

   # ---| Location and name of the data file to process |---------------------------------
   'data_loc': os.path.join(proj_dir,'reduc'), # relative path from this file loc !
   'data_fn': '0509-67_muse.fits',

   # ---| Generic locations for plots and products |-------------------------------------
   'plot_loc': os.path.join(proj_dir,'brutifus_plots') ,
   'prod_loc': os.path.join(proj_dir,'brutifus_products'),
   'tmp_loc': os.path.join(proj_dir,'brutifus_products','tmp'),
    
   'fn_list_fn': 'brutifus_fn-list_0509-67.pkl', # Name of the dictionary for filenames
   
   # ---| Constructing SNR maps |--------------------------------------------------------
   # lambda is in *REST* frame !
   'snr_ranges':[[7400., 8500., 'c'], # lam_min, lam_max, 'c'ontinuum, or 'e'mission
                 [6560., 6570., 'e'], # lam_min, lam_max, 'c'ontinuum, or 'e'mission
                ], 
   
   # ---| Manual sky subtraction |--------------------------------------------------------
   'sky_regions':[[290, 327,3], # x0,y0, radius, or x0,y0,dx,dy
                  [182,319,3],
                  [216,61,3],
                  [306,205,2],
                  [72,118,2],
                  [171,82,3],
                  [170,65,3],
                  [60,192,3],
                  [312,292,2],
                  [66,183,3],
                  [93,226,3],
                  [168,297,3],
                  [265,134,3],
                 ],
   
   # ---| Correcting for the Galactic extinction |---------------------------------------
   'gal_curve': 'f99', # Set to 'f99' to follow NED.
   'gal_rv': 3.1, # Set to 3.1 to follow NED.
   'Ab': 0.272, # Ab (Landolt) from NED.
   'Av': 0.206, # Av (Landolt) from NED.
   
   # ---| Continuum fitting |------------------------------------------------------------
   # Lowess specific variables
   #'lowess_snr_id': 0, # 0 = only check if spaxel exists, 1,2,3 = items in snr_ranges above
   #'lowess_snr_min': 0.5,  # What's the lowest SNR to run the LOWESS fit ?
   #'lowess_snr_max': None, # What's the highest SNR to run the LOWESS fit ? None = max
   'lowess_it': 10, # Number of iteration for sigma-clipping to get rid of outliers
   'lowess_frac': 0.15, # % of the array used for deriving each point. 0.05 = sweet spot?
   
   # ---| Emission line fitting | -------------------------------------------------------
   # How do I want to subtract the continuum: give the range of snr (in the continuum!)
   # and the corresponding technique.
   #'which_cont_sub':{'0->3':'lowess',
   #                  '3->max':'ppxf',
   #                 },
   
   # For the inital delta-v guess, the code looks at one single line. Which one is the
   # strongest in your case ? Give here the un-redshifted wavelength.
   #'ref_dv_line':o3h,
   
   # What kind of profile for emission lines ?
   # 'gauss_hist' = gaussian, accounting for bin size (the correct way).
   #'line_profile':'gauss', # 'gauss', 'gauss-herm'
   
   #'elines_snr_min': 3,  # What's the lowest SNR to run the line fit ?
   #'elines_snr_max': None, # What's the highest SNR to run the line fit ? None = max
   
   # What lines do I want to fit, and how ?
   # This is very close to defining the mpfit "parinfo" structure, but not quite.
   # For each line, first provide the reference wavelength, the offset from the "main"
   # component (NOT the total offset - used when fitting multiple gaussians to 1 line),
   # and the guess velocity dispersion (in km/s).
   # Then, individual dictionaries for the intensity, velocity and sigma of each lines
   # allow to define "fixed" parameters (that will NOT befitted), set lower/upper bounds,
   # and tie a parameters with another one (e.g. line intensity, velocity, etc ...)
   
   #'elines':{'a':[[ha,0,20,5], # Line reference wavelength, guess dispersion (km/s), min SNR
   #               {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, # Intensity >0
   #               {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':''}, # Velocity
   #               {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':''}, # Dispersion >0
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4
   #              ],
   #          'b':[[n2h,0,20,3], 
   #               {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(a)'}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(a)'},
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4 
   #              ],
   #          'c':[[n2l,0,20,3],
   #               {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':'1./3.*p(b)'}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(a)'}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(a)'},
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4 
   #              ],
   #          'd':[[s2h,0,20,3], 
   #               {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(a)'}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(a)'}, 
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4
   #              ],
   #          'e':[[s2l,0,20,3], 
   #               {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(a)'}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(a)'}, 
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4
   #              ],
   #          'f':[[hb,0,20,3],
   #               {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(a)'},
   #               {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':''},
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4 
   #              ],
   #          'g':[[o3h,0,20,3], 
   #               {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':''}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(a)'}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(f)'},
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4 
   #              ],
   #          'h':[[o3l,0,20,3], 
   #               {'fixed':0, 'limited':[1,0], 'limits':[0.,0.], 'tied':'1./2.98*p(g)'},
   #               {'fixed':0, 'limited':[1,1], 'limits':[6900.,7600.], 'tied':'p(a)'}, 
   #               {'fixed':0, 'limited':[1,1], 'limits':[1.,100.], 'tied':'p(f)'},
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h3
   #               {'fixed':1, 'limited':[0,0], 'limits':[0.,0.], 'tied':''}, # h4 
   #              ],
   #         },
            
   # ---| Structure detection and aperture extraction |----------------------------------
   # What line flux do I want to use to extract my aperture map ? Use a list of the 
   # elines keys for identification.
   #'ap_map_lines':['a','g'],
   
   # The default radius for the apertures in the automated peak detection routine, in pixel.
   # Choosing a value corresponding to the seeing is a good first guess.
   #'ap_radius':3.0,             
   
   # ---| Extragalactic reddening/attenuation correction |------------------------------- 
   #'hahb_0': 2.85, # theoretical reference flux ratio of Halpha/Hbeta. 2.85 for Case B 
                   # recombination, 3.1 for AGNs.
   #'egal_curve': 'fd05', # Reddening/attenuation law to use. Valid values are:
                         # 'fd05': Fischera & Dopita (2005) use fv=3.08 and rva=4.5 for a 
                         #         curve similar to Calzetti (2000).
                         # 'cal00': Calzetti (2000), use with rv=4.05 for default.
                         # 'cal00*': Calzetti (2000), with 0.44 correction factor.
                         # 'f99': Fitzpatrick (1999)
                         # 'ccm89': Cardelli, Clayton & Mathis (1989)
   #'egal_rv': 3.08, # The value of Rv. Use 3.08 for a Calzetti-like law, when using 'fd05'.
   #'egal_rva':4.3,  # The value of Rva, if using 'fd05'.   
   
   # ---| fit_kinematic_pa |-------------------------------------------------------------
   
   #'fit_kin': ['ppxf','a'], # Which velocity map to fit ? ppxf = stars, 'a' = elines codes
   #'fit_kin_nsteps': 361, # The angular resolution (361= default=0.5 degrees)
   #'fit_kin_center': [146, 123], # In pixels, the location of the kinematic center
   
}
  
# Export this as a pickle file    
f = open(pickle_fn, 'wb')
pickle.dump(brutifus_params, f)
f.close()  




          
