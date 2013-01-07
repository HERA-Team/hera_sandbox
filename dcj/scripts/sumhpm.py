#!/usr/bin/env python
#
#  sumhpm.py
#  
#
#  Created by Danny Jacobs on 2/15/11.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
import healpy as hpy

#o = optparse.OptionParser()
#o.add_option('--smooth', dest='smooth', default=0.25,
#    help='Smooth the weights to this scale. [0.25 deg,fwhm]')
#o.add_option('--wcut', default=99.0,
#    help='Smooth the weights to this scale. [0.25 deg,fwhm]')
#
#opts, args = o.parse_args(sys.argv[1:])



cals = {'psa455_smap_v1.fits':{'px':2506187,'flux':72.,'name':'0521-365'},
#        'psa336_v11.fits':{'px':2298059,'flux':21.,'name':'1419-272'}}
        'psa336_aipy_dev_smap.fits':{'px':2359502,'flux':19,'name':'1422-297'}}
args = cals.keys()
smooth = 0.75
wcut = 99.0
beam_sigma = 2.5 #this parameter is really really important and can screw everything up if its wrong. FYI
outfile = 'psa455_336_cal_v2.fits'
nside = 512
Weight,Flux = None,None
for file in args:
    px = cals[file]['px']
    calflux = cals[file]['flux']
    F = hpy.read_map(file,field=0)
    W = hpy.read_map(file,field=1)
    if hpy.get_nside(F) != nside:
        print "changing nside of %s from %d to %d"%(file,hpy.get_nside(F),nside)
        F = hpy.ud_grade(F,nside)
        W = hpy.ud_grade(W,nside)
    #smooth and mask the weights to get rid of insane pixels
#    W = hpy.smoothing(W,fwhm=smooth,degree=True)
#    WPDF,B = n.histogram(W,bins=n.sqrt(len(W)))
#    WCDF = n.cumsum(WPDF)/n.sum(WPDF)
#    wcutvali = n.argwhere(n.abs(WCDF-wcut/100.)==n.min(n.abs(WCDF-wcut/100.))).squeeze()
#    wcutval = B[wcutvali]
#    print "masking %3.1f %% of the data where weights are below %f"%(len(WCDF[WCDF<wcut/100.])/len(WCDF)*100,wcutval)
#    wcutval = 1/exp(beam_sigma) * W.max()
#    print "masking %3.1f %% of the data where weights are below %f"%(len(W<wcutval)/len(W)*100,wcutval)
#    W = n.ma.masked_where(W<B[wcutval],W)
    Sky = F/W
    
    if Flux is None:
        Flux = n.zeros_like(F)
        Weight = n.zeros_like(F)
    print "calibrating %s at %4.1f Jys in %s to be %4.1f Jys"%(cals[file]['name'],
        Sky[px],file,calflux)
#    Flux += F * n.sqrt(calflux/Sky[px])
#    Weight += W / n.sqrt(calflux/Sky[px])
    Flux  += F * calflux/Sky[px]
    Weight += W

print "WARNING: there are %d NaNs in the sky"%(n.sum(n.isnan(Flux/Weight)),)    
hpy.write_map(outfile,n.vstack([Flux,Weight]),colnames=['TEMPERATURE','WEIGHT'])
