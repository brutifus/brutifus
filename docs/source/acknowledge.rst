.. |DOI_latest| image:: https://zenodo.org/badge/157203434.svg
   :target: https://zenodo.org/badge/latestdoi/157203434
.. |ASCL| image:: https://img.shields.io/badge/ascl-1903.004-blue.svg?colorB=262255
   :target: http://ascl.net/1903.004 
   
Acknowledging brutifus
======================

1. Only use lower case letters when mentioning brutifus, and always include the release number.
Ideally, you should also include a) the DOI associated with any of the Github releases, 
and b) the code's ASCL entry number, e.g.:

    brutifus |release|: |DOI_latest| |ASCL|

2. If you use brutifus for your data analysis (and remember that you did so by the time you
   reach the publication stage!), please cite:
   
   Vogt, *brutifus: Python module to post-process datacubes from integral field spectrographs*,
   ASCL 1903.004 (2019). `ADS entry <http://adsabs.harvard.edu/abs/2019ascl.soft03004V>`_
 

   brutifus also uses several packages that **should also be acknowledged in their own right.** 
   The following Tex-formatted acknowledgment is one way to do so::

    This research has made use of \textsc{brutifus}, a Python module to process data cubes 
    from integral field spectrographs (Vogt, 2019). \textsc{brutifus} relies on 
    \textsc{statsmodel} (Seabold & Perktold 2010),
    \textsc{matplotlib} (Hunter 2007), \textsc{astropy}, a community-developed core Python 
    package for Astronomy (Astropy Collaboration et al., 2013, 2018), and \textsc{photutils}, 
    an affiliated package of \textsc{astropy} for photometry (DOI:10.5281/zenodo.2533376).

   Finally, you also ought to cite the following works, depending on your use of brutifus:

    a) Cleveland (1979): 
        the reference for the Locally Weighted Scatterplot Smoothing (LOWESS) algorithm used 
        by brutus (via statsmodels) to fit the continuum.
            
    b) The reddening laws:
        Either the Cardelli, Clayton & Mathis (1989) law, the Calzetti et al. (2000) law or 
        the theoretical model of a turbulent dust screen of **Fischera & Dopita (2005) 
        [default]** for the extragalactic attenuation corrections, and
        the **Fitzpatrick (1999) law [default]** for the galactic extinction.
        
        If you use the extinction values :math:`A_V` and :math:`A_B` from the NASA 
        Extragalactic Database (NED) to correct for the Galactic extinction [default], then
        to be thorough, you should mention that::
        
            The Galactic extinction is derived using NED from the Schlafly & Finkbeiner 
            (2011) recalibration of the Schlegel, Finkbeiner & Davis (1998) infrared-based 
            dust map. The map is based on dust emission from COBE/DIRBE and IRAS/ISSA; 
            the recalibration assumes a Fitzpatrick (1999) reddening law with Rv = 3.1 and 
            different source spectrum than Schlegel, Finkbeiner & Davis (1998).
        
        and you also ought to acknowledge NED itself::
        
            This research has made use of the NASA/IPAC Extragalactic Database (NED) 
            which is operated by the Jet Propulsion Laboratory, California Institute of 
            Technology, under contract with the National Aeronautics and Space Administration. 
        
References:
 - `Astropy Collaboration et al. (2013) <http://cdsads.u-strasbg.fr/abs/2013A%26A...558A..33A>`_
 - `Astropy Collaboration et al. (2018) <http://adsabs.harvard.edu/abs/2018arXiv180102634T>`_
 - `Cardelli, Clayton & Mathis (1989) <http://adsabs.harvard.edu/abs/1989ApJ...345..245C>`_
 - `Calzetti et al. (2000) <http://adsabs.harvard.edu/abs/2000ApJ...533..682C>`_
 - Cleveland (1979)::
    
    @article{doi:10.1080/01621459.1979.10481038,
        author = { William S.   Cleveland },
        title = {Robust Locally Weighted Regression and Smoothing Scatterplots},
        journal = {Journal of the American Statistical Association},
        volume = {74},
        number = {368},
        pages = {829-836},
        year = {1979},
        doi = {10.1080/01621459.1979.10481038},
    }
 
 - `Fischera & Dopita (2005) <http://adsabs.harvard.edu/abs/2005ApJ...619..340F>`_
 - `Hunter (2007) <http://cdsads.u-strasbg.fr/abs/2007CSE.....9...90H>`_    
 - Seabold & Perktold (2010)::
 
    @inproceedings{seabold2010,
        title={Statsmodels: Econometric and statistical modeling with python},
        author={Seabold, Skipper and Perktold, Josef},
        booktitle={9th Python in Science Conference},
        year={2010},
    }
    
 - `Schlafly & Finkbeiner (2011) <http://adsabs.harvard.edu/abs/2011ApJ...737..103S>`_  
 - `Schlegel, Finkbeiner & Davis (1998) <http://adsabs.harvard.edu/abs/1998ApJ...500..525S>`_ 
 - `Vogt (2019) <http://adsabs.harvard.edu/abs/2019ascl.soft03004V>`_
