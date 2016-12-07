#! /usr/bin/env python

"""
Plot first-cal phase solutions as a function of
pysical antenna position
"""

import os, optparse, numpy as np, aipy, sys
from matplotlib import pyplot as plt
from matplotlib import colors
import matplotlib.cm as cmx

o = optparse.OptionParser()
o.set_description(__doc__)
o.set_usage('map_fcPhase.py [options] *fc.npz')
aipy.scripting.add_standard_options(o, pol=True, cal=True)
o.add_option('-b','--badants',dest='b',default=None,help='Comma separated list of antennae to exlude')
#TODO: badants file as input
opts,args = o.parse_args(sys.argv[1:])

#Array data
print 'reading, %s'%opts.cal
exec("import {calfile} as cal".format(calfile=opts.cal))
antpos = cal.prms['antpos']
antpos = cal.prms['antpos']
nants = len(antpos.keys())
if not opts.b is None: ba = opts.b.split(',')
else: ba = []

x,y = [],[]
for ant in antpos.keys():
    if str(ant) in ba: continue
    if ant>111: continue #avoid outriggers
    pos = antpos[ant]
    x.append(pos['top_x'])
    y.append(pos['top_y'])

XMIN,XMAX = 50+5.408e5,320+5.408e5
YMIN,YMAX = 6.60113e6,30+6.60113e6
VMIN,VMAX = -30,30
for npz in args:
    nm = os.path.basename(npz)
    d=[]
    for ant in antpos.keys():
        if str(ant) in ba: continue
        if ant>111: continue
        data = np.load(npz)
        d.append(data[str(ant)+'d'][0])
    f,axarr = plt.subplots(1,2,sharey=True)
    
    axarr[0].plot(x,d,'ko')
    axarr[0].set_xlim(XMIN,XMAX)
    axarr[0].set_ylim(VMIN,VMAX)
    axarr[0].set_xlabel('E-W position',size=15)
    axarr[0].set_ylabel('Delay (ns)',size=15)
    
    axarr[1].plot(y,d,'ko')
    axarr[1].set_xlim(YMIN,YMAX)
    axarr[1].set_xlabel('N-S position',size=15)
    
    #sc=plt.scatter(x,y,c=d,s=40,vmin=VMIN,vmax=VMAX)
    #plt.colorbar(sc)
    #plt.suptitle(nm)
    plt.savefig('/home/saulkohn/testimg_xy/%s.png'%nm)#TODO separate directory
    #plt.show()
    plt.close()

#import IPython;IPython.embed()
