#!/usr/bin/env python

"""

NAME: 
      update_comp_list.py
PURPOSE: 
      Updates a list of .uv files to be compressed based on whether they have already been successfully compressed.
EXAMPLE CALL:
      ./update_comp_list.py --filenames psa6620xx01
AUTHOR:
      Carina Cheng

"""

import aipy
import numpy
import pylab
import pyfits
import matplotlib.pyplot as plt
import optparse
import os, sys
import glob

o = optparse.OptionParser()
o.set_usage('update_comp_list.py [options]')
o.set_description(__doc__)
o.add_option('--filenames', dest='filenames', default='psa6620xx01',
             help='Name of directory containing filenames.')
opts, args = o.parse_args(sys.argv[1:])

c_files = []

for cfile in glob.glob(os.path.join(opts.filenames[:7]+'/*.uvcRRE')):
	c_files.append(cfile[8:])

files = open(opts.filenames)


num = int(opts.filenames[-1])
newfile = open(opts.filenames[:-1]+str(num+1),'w') 

for f in files:
	
	fname = f.split('/')[-1][:-1]+'cRRE'

	#fprefix = '/'.join(f.split('/')[:-1])

	if fname in c_files:
		continue
	else:
		newfile.write(f)
		
	
