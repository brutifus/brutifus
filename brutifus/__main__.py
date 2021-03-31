# -*- coding: utf-8 -*-
'''
brutifus: a set of Python modules to process datacubes from integral field spectrographs.\n
Copyright (C) 2018-2020,  F.P.A. Vogt
Copyright (C) 2021, F.P.A. Vogt & J. Suherli
All the contributors are listed in AUTHORS.

Distributed under the terms of the GNU General Public License v3.0 or later.

SPDX-License-Identifier: GPL-3.0-or-later

This file contains some convenience functions for brutifus.

Created July 2019, F.P.A. Vogt - frederic.vogt@alumni.anu.edu.au
'''
# --------------------------------------------------------------------------------------------------

import os
import shutil
import argparse

from . import brutifus_metadata as bifus_m
from .brutifus_version import __version__
from . import brutifus as bifus

# Use argparse to make brutifus user friendly ------------------------------------------------------
parser = argparse.ArgumentParser(description=''' Aids in the post-processing of 3D datacubes. ''',
                                 epilog=' Full documentation: %s' %
                                 ('http://brutifus.github.io/'),
                                 formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-v', '--version', action='version', version=('%s' % (__version__)))

parser.add_argument('-s', '--setup', action='store_true',
                    help='Create a local copy of the parameter and execution files '+
                    'needed to run brutifus.')

parser.add_argument('-e', '--execute', action='store', nargs=2,
                    metavar=('procsteps_file', 'params_file'),
                    help='run the steps specified in "procsteps_file" ' +
                    'with the high-level parameters specified in "params_file".')

def main():
    ''' The main function, run when brutifus is started with the high-level entry point. '''

    # What did the user type in ?
    args = parser.parse_args()

    # Setup the scene ?
    if args.setup:

        print('')

        # Very well, let's copy the files at the current location
        for f in [bifus_m.bifus_params, bifus_m.bifus_procsteps]:
            if not os.path.isfile(os.path.join('.', f)):

                print('   - creating the local copy of %s ...' % (f))
                shutil.copyfile(os.path.join(bifus_m.bifus_dir, 'exec_scripts', f), f)

        # Let's also create the working directories if needed
        for d in [bifus_m.plot_loc, bifus_m.prod_loc]:
            if not os.path.isdir(os.path.join('.', d)):

                print('   - creating the local directory "%s" ...' % (d))
                os.mkdir(os.path.join('.', d))

        # Don't forget the tmp one
        tmp_dir = os.path.join('.', bifus_m.prod_loc, bifus_m.tmp_loc)
        if not os.path.isdir(tmp_dir):

            print('   - creating the local directory "%s" ...' % (tmp_dir))
            os.mkdir(tmp_dir)

        print('')
        print('brutifus setup complete.')

    if args.execute:
        # Start the processing
        bifus.run(args.execute[0], args.execute[1])

# Start of the interactive part --------------------------------------------------------------------
if __name__ == "__main__":

    main()
