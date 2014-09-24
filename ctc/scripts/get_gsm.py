#!/usr/bin/env python2

"""
Author: Isaac Domagalski
Email: idomagalski@berkeley.edu
Get GSM data in FITS format.
"""

import os
import sys
import aipy
import glob
import argparse
import numpy as np

if __name__ == '__main__':
    # Read options from the command line
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--lower-bound', required=True,
                        help='Starting frequency of the GSM.')
    parser.add_argument('-s', '--step-size', required=True,
                        help='Step in frequencies.')
    parser.add_argument('-n', '--num-steps', required=True,
                        help='Number of frequency steps.')
    parser.add_argument('-p', '--prefix', required=True,
                        help='Filename prefix of GSM files.')
    parser.add_argument('-o', '--outdir', required=True,
                        help='Output directory to store the files to.')
    parser.add_argument('-d', '--save-dat', action='store_true',
                        help='Save intermediate .dat files.')
    args = parser.parse_args()

    # Store the options
    start_freq = args.lower_bound
    step_size  = args.step_size
    num_steps  = args.num_steps
    prefix     = args.prefix
    outdir     = os.path.abspath(args.outdir)

    # Change to the directory containing the fortran binaries.
    os.chdir(os.path.dirname(os.path.abspath(sys.argv[0])))

    # Create the arguments file for the fortran routine
    with open('args.dat', 'w') as f:
        f.write(' '.join([prefix, start_freq, step_size, num_steps]) + '\n')

    # Generate the GSM
    try:
        if os.system('./gsmmf.sh'):
            print 'ERROR: Cannot generate GSM files.'
            sys.exit(1)

        # Move the data files to the output directory.
        os.system('mkdir -p ' + outdir)
        os.system('mv ' + ' '.join(glob.glob(prefix + '*.dat')) + ' ' + outdir)

        # Convert the data files to FITS format.
        for name in glob.glob(os.path.join(outdir, prefix + '*.dat')):
            fitsname = name.replace('.dat', '.fits')
            h = aipy.healpix.HealpixMap(nside=512)
            h.map = np.loadtxt(name)
            h.to_fits(fitsname)
            if not args.save_dat:
                os.system('rm -f ' + name)
            print os.path.basename(fitsname) + ' completed'

        # Move the args.dat file to the output directory
        os.system('mv args.dat ' + outdir)

    except KeyboardInterrupt:
        print 'GSM generation terminated by the user.'
        sys.exit(1)
