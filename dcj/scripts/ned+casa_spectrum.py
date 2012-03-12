#!/usr/bin/env python
#
#  ned+casa_spectrum.py
#  
#
#  Created by Danny Jacobs on 1/17/10.
#  PAPER Project
#

import aipy as a, numpy as n,math as m
import sys, optparse
from pylab import *

o = optparse.OptionParser()
o.set_usage("ned+casa_spectrum.py ned_spectrum.txt casa_spectrum.txt")
#a.scripting.add_standard_options(o)
o.add_option('--freq_range',dest='freq_range',default=None,
    help="Specify a frequency range as fmin_fmax in Hz. [None]")
o.add_option('--ylmin',dest='ylim',default=None,
    help='Speficy a ylimit')
opts, args = o.parse_args(sys.argv[1:])
if not opts.freq_range is None: freq_range = map(float,opts.freq_range.split('_'))

nedfile = args[0]
casafile = args[1]
nedlines = open(nedfile).readlines()

#get the column labels in the ned spectrum file
for l in nedlines:
    if l[0:2]=='No':
        nedcols = dict(zip(l.split('\t'),range(len(l.split('\t')))))
        break
for l in nedlines:
    if l.startswith('Photometric'):
        srcname = l[21:]

nedlinesplit = [l.split('\t') for l in nedlines]
#print freq_range
#print nedlinesplit[30][nedcols['Frequency']]>freq_range[0] and nedlinesplit[30][nedcols['Frequency']]<freq_range[1]
#print  nedlinesplit[30][nedcols['Frequency']]
if not opts.freq_range is None:
    nedspec = [(l[nedcols['Frequency']],l[nedcols['NED Photometry Measurement']],l[nedcols['NED Uncertainty']][2:]) for \
    l in nedlinesplit if len(nedcols)==len(l) and l[nedcols['NED Units']].startswith('Jy')\
    and l[nedcols['NED Photometry Measurement']]!='' and  l[nedcols['Uncertainty']]!=''\
     and float(l[nedcols['Frequency']])>freq_range[0] and float(l[nedcols['Frequency']])<freq_range[1]]
else: 
    nedspec = [(l[nedcols['Frequency']],l[nedcols['NED Photometry Measurement']],l[nedcols['NED Uncertainty']][2:]) for \
        l in nedlinesplit if len(nedcols)==len(l) and l[nedcols['NED Units']].startswith('Jy')\
        and l[nedcols['NED Photometry Measurement']]!='' and  l[nedcols['Uncertainty']]!='']
for i in range(len(nedspec)):
    if nedspec[i][2]=='': nedspec[i]= (nedspec[i][0],nedspec[1][1],0)
print nedspec
nedspec = n.array(nedspec).astype(n.float)
nedspec[:,2]= n.abs(n.array(nedspec[:,2]))

casalines = open(casafile).readlines()
casaspec = n.array([map(float,l.split()) for l in casalines if l[0]!='#']).astype(n.float)
casaspec[:,0] *= 1e9
figure()
axes(xscale='log',yscale='log')
errorbar(nedspec[:,0],nedspec[:,1],yerr=nedspec[:,2],fmt='+')
#plot(casaspec[:,0],casaspec[:,1],'.',alpha=.5,markersize=2)
errorbar(n.average(casaspec[:,0]),n.average(casaspec[:,1]),yerr=n.sqrt(n.var(casaspec[:,1])),fmt='x',lw=2)
#xlabel(input("xlabel:"))
#ylabel(input("ylabel:"))
xlabel('Frequency [Hz]')
ylabel('flux [Jy]')
if not opts.ylim is None:
    ylim(float(opts.ylim.split('_')[0]),float(opts.ylim.split('_')[0]))
else:
    ylim([n.min(nedspec[:,1])/10,n.max(nedspec[:,1])*10])
title(srcname)
ioff()
show()