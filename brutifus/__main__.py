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
from . import brutifus as bifus

import os
import shutil
import argparse

# Use argparse to make brutifus user friendly ---------------------------------------------
parser = argparse.ArgumentParser(description=''' Aids in the post-processing of 3D datacubes. ''',
                                 epilog =' Full documentation: %s' % ('http://fpavogt.github.io/brutifus'),
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-v', '--version', action='version', version=('%s'%__version__))                   
                    
parser.add_argument('-s', '--setup', action = 'store_true', 
                    help = 'Create a local copy of the parameter and execution files '+
                            'needed to run brutifus.')

parser.add_argument('-e', '--execute', action='store', nargs = 2, 
                    metavar = ('procsteps_file', 'params_file'),
                    #type=open, 
                    help='run the steps specified in "procsteps_file" '+
                         'with the high-level parameters specified in "params_file".') 

def main():
   # What did the user type in ?
   args = parser.parse_args()
   
   # Setup the scene ?
   if args.setup:
      
      print('')
      
      # Very well, let's copy the files at the current location
      for f in [bifus_m.bifus_params, bifus_m.bifus_procsteps]:
         if not os.path.isfile(os.path.join('.',f)):
          
            print ('   - creating the local copy of %s ...' % (f))
            shutil.copyfile(os.path.join(bifus_m.bifus_dir,'exec_scripts',f), f)
            
         
      # Let's also create the working directories if needed
      for d in [bifus_m.plot_loc, bifus_m.prod_loc]:
         if not os.path.isdir(os.path.join('.',d)):
            
            print ('   - creating the local directory "%s" ...' % (d))
            os.mkdir(os.path.join('.',d))
            
      # Don't forget the tmp one
      tmp_dir = os.path.join('.', bifus_m.prod_loc, bifus_m.tmp_loc)
      if not os.path.isdir(tmp_dir):
      
         print ('   - creating the local directory "%s" ...' % (tmp_dir))
         os.mkdir(tmp_dir)
       
             
      print('')
      print('brutifus setup complete.')
   
   if args.execute:
      # Start the processing
      bifus.run(args.execute[0], args.execute[1])

   
# Start of the interactive part ----------------------------------------------------------
if __name__ == "__main__":

   main()
            
            
      
      
   
   