# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018-2020,  F.P.A. Vogt
Copyright (C) 2021, F.P.A. Vogt & J. Suherli
All the contributors are listed in AUTHORS.

Distributed under the terms of the GNU General Public License v3.0 or later.

SPDX-License-Identifier: GPL-3.0-or-later

This file contains functions related to the continuum fitting inside brutifus.

Created November 2018, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# --------------------------------------------------------------------------------------------------

import numpy as np
from statsmodels.nonparametric.smoothers_lowess import lowess

# --------------------------------------------------------------------------------------------------

def lowess_fit(spec, lams, frac=0.05, it=5):
    ''' Fit a spectrum using a Locally Weighted Scatterplot Smoothing approach.

    Wraps around statsmodels.nonparametric.smoothers_lowess.lowess().

    Args:
        spec (ndarray): The input spectrum.
        lams (ndarray): The corresponding wavelength array.
        frac (float, optional): Between 0 and 1. The fraction of the data used when estimating each
            y-value. See the statsmodel lowess function for details. Defaults to 0.05
        it (int, optional): The number of residual-based reweightings to perform.
            See the statsmodel lowess function for details. Defaults to 5.

    Returns:
        ndarray: The fitted array, with size equal to spec.

    .. note:: This function fits a spectrum using a LOWESS (Locally Weighted Scatterplot
              Smoothing) technique, described in:
              Cleveland, W.S. (1979) Robust Locally Weighted Regression and Smoothing
              Scatterplots. Journal of the American Statistical Association 74 (368): 829-836.

              This is robust to outliers (hot pixels, cosmics), and is also efficient to
              ignore emission lines. frac=0.05 and it=5 seem to work very fine for spectra
              of supernova remnants, including for those with strong foreground/background stellar
              continuum.

    .. warning:: Users should not forget that this method will NOT remove stellar absorption lines !

    .. warning:: If you have broad emission lines in your spectra, you will want to be *very*
               careful with LOWESS !!!

    '''

    # Only do the fit if there is some signal. Avoid an ugly warning in the prompt.
    if all(np.isnan(spec)):
        fit = np.zeros_like(spec) * np.nan
    else:
        fit = lowess(spec, lams, frac=frac, it=it, is_sorted=True, missing='drop',
                     return_sorted=False)

    return fit

# --------------------------------------------------------------------------------------------------
