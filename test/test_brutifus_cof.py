# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018-2020,  F.P.A. Vogt
Copyright (C) 2021, F.P.A. Vogt & J. Suherli
All the contributors are listed in AUTHORS.

Distributed under the terms of the GNU General Public License v3.0 or later.

SPDX-License-Identifier: GPL-3.0-or-later

This file contains test functions related to brutifus_cof.py

Created November 2018, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''

# Import from python
import numpy as np

# Import from brutifus
from brutifus.brutifus_cof import lowess_fit

def test_lowess_fit():
    """ Tests the lowess fit function """

    out = lowess_fit(np.arange(10), np.ones(10)*np.nan)

    assert all(np.isnan(out))
    assert len(out) == 10
