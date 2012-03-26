#!/usr/bin/env python
#
#  find_pointing.py
#  
#
#  Created by Danny Jacobs on 1/11/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse, pyfits as pf,ephem

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
o.add_option('-v',dest="verbose",action='store_true',
    help="Print more")
o.add_option('-J',dest="J",action='store_true',
    help="Use truncated J names.")
o.add_option('--sep',dest='sep',
    help="Seperation character [CR]")
    
opts, args = o.parse_args(sys.argv[1:])

def hname(ra):
	one = str((int(str(ephem.hours(ra*n.pi/180)).split(':')[0])+100000))[-2::]
	two =  str((int(str(ephem.hours(ra*n.pi/180)).split(':')[1])+100000))[-2::]
	ret = one + two
	return ret
def dname(dec):
	one = str((int(str(ephem.degrees(dec*n.pi/180)).split(':')[0])+100000))[-2::]
	two = str((int(str(ephem.degrees(dec*n.pi/180)).split(':')[1])+100000))[-2::]
	if dec < 0: add = '-'
	else: add = '+'
	ret = add + one + two
	return ret

for file in args:
    data,kwds = a.img.from_fits(file)
    if opts.verbose: print file,
    if opts.J: print 'J'+hname(kwds['ra'])+dname(kwds['dec']),
    else: 
        if kwds['dec']>0: print str(ephem.hours(kwds['ra']*n.pi/180))+\
                '+'+str(ephem.degrees(kwds['dec']*n.pi/180)),
        else:  print str(ephem.hours(kwds['ra']*n.pi/180)) + \
                str(ephem.degrees(kwds['dec']*n.pi/180)),
    if not opts.sep is None: print opts.sep,
    else: print
    