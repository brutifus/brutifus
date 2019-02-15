# -*- coding: utf-8 -*-
'''
 brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
 Copyright (C) 2018-2019,  F.P.A. Vogt
 
 -----------------------------------------------------------------------------------------
 
 This file contains functions related to the continuum fitting inside brutifus. 

 Created November 2018, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# ----------------------------------------------------------------------------------------

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

# ----------------------------------------------------------------------------------------

def lowess_fit(spec, lams, frac=0.05, it=5):
   ''' Fit a spectrum using a Locally Weighted Scatterplot Smoothing approach.
    
   Wraps around statsmodels.nonparametric.smoothers_lowess.lowess().
    
   :param spec: 1-D numpy array
            The input spectrum.
   :param lams: 1-D numpy array
            The corresponding wavelength array.
   :param frac: float [default:0.05]
            Between 0 and 1. The fraction of the data used when estimating each y-value.
            [From the statsmodel lowess function]
   :param it: int [default:5]
         The number of residual-based reweightings to perform.
         [From the statsmodel lowess function]
     
   :return: The fitted array, with size equal to spec.  
   :rtype: numpy array 
    
   .. note:: This function fits a spectrum using a LOWESS (Locally Weighted Scatterplot 
             Smoothing) technique, described in: 
             Cleveland, W.S. (1979) Robust Locally Weighted Regression and Smoothing 
             Scatterplots. Journal of the American Statistical Association 74 (368): 829-836.
    
             This is robust to outliers (hot pixels, cosmics), and is also efficient to 
             ignore emission lines. frac=0.05 and it=5 seem to work very fine for spectra 
             of any SNR, both lousy with no continuum, and good ones in the center of 
             galaxies - modulo the stellar absorption features which are of course 
             "ignored" by the LOWESS routine.
   
   ..warning:: If you have broad emission line in the spectrum, you want to be *very*
               careful with LOWESS !!!
               
   '''
    
   # Only do the fit if there is some signal. Avoid an ugly warning in the prompt.
   if np.all(np.isnan(spec)):
      fit = np.zeros_like(spec) * np.nan
   else:
      fit = lowess(spec,lams,frac=frac, it=it, is_sorted=True, missing = 'drop', 
                   return_sorted=False)
                                            
   return fit
   
# ----------------------------------------------------------------------------------------      
