#!/usr/bin/env python
#
#  pspec_all_sky.py
#  
#
#  Created by Danny Jacobs on 1/26/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
import healpy as hp
from pylab import *
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
o.add_option('--no_cache',dest='no_cache',
   help = """Dont try and load cached cls from a file for listed files.  
   Force recalculation. eg  --no_cache=summap.fits,vlss""")
o.add_option('--conv',
   help="binary toggle for conversion from Jys to T. Each bit corresponds with a file.")

opts, args = o.parse_args(sys.argv[1:])
if not opts.conv is None:
    opts.conv = n.array(map(int,opts.conv.split(',')))
else:
    opts.conv = n.zeros(len(args))
print opts.conv
cls ={}
if not opts.no_cache is None: no_cache = opts.no_cache.split(',')
else:no_cache = []

for i,file in enumerate(args):
    print file
    name = file.split('/')[-1]
    try:
        if not name in no_cache and not ('all' in no_cache): 
            s_cl = n.loadtxt(name+'.cl')
            print "loading",name+'.cl'
    except(IOError):
        no_cache.append(name)
    if name in no_cache or 'all' in no_cache:
        s = hp.read_map(file)
        print "computing alms"
        s_alm = hp.map2alm(s)
        s_cl = hp.alm2cl(s_alm)
        n.savetxt(name+'.cl',s_cl)
    cls[file] = s_cl * 0.002**(opts.conv[i])
    
    
figure(88)    
for cl in cls:
    s_cl = cls[cl]
    loglog(s_cl*n.linspace(0,len(s_cl),len(s_cl))**2,label=cl)
legend(loc='lower right')
xlabel('$\ell$')
ylabel('$\ell^2 C_\ell$ $[K^2]$')

figure(89)    
for cl in cls:
    s_cl = cls[cl]
    loglog(s_cl,label=cl)
legend(loc='lower right')
xlabel('$\ell$')
ylabel('$C_\ell$ $[K^2]$')


print  "done"


show()