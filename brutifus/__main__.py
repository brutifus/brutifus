# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018-2019,  F.P.A. Vogt
 
-----------------------------------------------------------------------------------------
 
This file contains some convenience functions for brutifus. 

Created July 2019, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
 
-----------------------------------------------------------------------------------------
  
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
'''

from . import brutifus_metadata as bifus_m
from .brutifus_version import __version__

import os
import shutil
import argparse

# Use argparse to make brutifus user friendly ---------------------------------------------
parser = argparse.ArgumentParser(description=''' Aids in the post-processing of 3D datacubes. ''',
                                 epilog =' Full documentation: %s \n \n \
                                           Feedback, questions, comments: \
                                           frederic.vogt@alumni.anu.edu.au \n' % ('http://fpavogt.github.io/brutifus'))

parser.add_argument('--version', action='version', version=('%s'%__version__))                   
                    
parser.add_argument('--setup', action = 'store_true', 
                    help = 'Create a local copy of the parameter and execution files needed to run brutifus.',)


def main():
   # What did the user type in ?
   args = parser.parse_args()
   
   # Setup the scene ?
   if args.setup:

      # Very well, let's copy the files at the current location
      for f in ['brutifus_params.py', 'brutifus_execute.py']:
         if not os.path.isfile(os.path.join('.',f)):
          
            shutil.copyfile(os.path.join(bifus_m.bifus_dir,'exec_scripts',f), f)
         
      # Let's also create the working directories if needed
      for d in ['brutifus_plots','brutifus_products']:
         if not os.path.isdir(os.path.join('.',d)):
            
            os.mkdir(os.path.join('.',d))
            
      # Don't forget the tmp one
      tmp_dir = os.path.join('.','brutifus_products','tmp')
      if not os.path.isdir(tmp_dir):
         os.mkdir(tmp_dir)
             


# Start of the interactive part ----------------------------------------------------------
if __name__ == "__main__":

   main()
            
            
      
      
   
   