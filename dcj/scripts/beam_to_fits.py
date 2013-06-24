#!/usr/bin/env python
#
#  beam_to_fits.py
#  
#
#  Created by Danny Jacobs on 11/3/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse
from pylab import *

"""
Generate a fits spectral cube of the beam.
"""

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--freq_range',default='110_190',
    help='A range of frequencies in MHz [default=110_190]')
o.add_option('--nchan',default=10,type='int',
    help='Number of channels [default = 10]')    
o.add_option('--pol',default='x',
    help='Polarization default=x')
opts, args = o.parse_args(sys.argv[1:])

#    a.img.to_fits(filename, i, clobber=True,
#        object=src.src_name, obs_date=str(aa.date),
#        ra=cen.ra*a.img.rad2deg, dec=cen.dec*a.img.rad2deg, epoch=2000.,
#        d_ra=L[-1,-1]*a.img.rad2deg, d_dec=M[1,1]*a.img.rad2deg,
#        freq=n.average(aa[0].beam.afreqs)*1e9,history=history,
#        axes=('ra---sin','dec--sin'))


flim = n.array(map(float,opts.freq_range.split('_')))/1000
freqs = n.linspace(flim[0],flim[1],num=opts.nchan)
print "beam_to_fits starting"
print "using freqs"
print freqs
print "loading beam model"
aa = a.cal.get_aa(opts.cal, freqs)
uvsize = 100
uvres = 0.5
print "generating image"
im = a.img.Img(size=uvsize, res=uvres)
res = 1/float(uvsize) * a.img.rad2deg


tx,ty,tz = im.get_top()
print tx.shape #200,200
invalid = tx.mask.copy()
tx = tx.filled(0).flatten()
ty = ty.filled(0).flatten()
tz = tz.filled(0).flatten()
print "evaluating beam at image coords"
resp = aa[0].bm_response((tx,ty,tz), pol=opts.pol)
print "found coords of shape:",resp.shape 
invalid = n.array([invalid]*len(freqs))
resp.shape = invalid.shape
resp = n.where(invalid, 0, resp/resp[opts.nchan/2,0,0])
print "flagging NaNs"
print resp.shape
history = ' '.join(sys.argv)
resp = n.array([a.img.recenter(R,[resp.shape[1]/2,resp.shape[2]/2]) for R in resp])
print "saving file %s"%args[0]
a.img.to_fits(args[0], resp, clobber=True,
    object='beam', obs_date='2455504',
    ra=0, dec=0, epoch=2000.,
    d_ra=res, d_dec=res,d_freq=(flim[1] - flim[0])/opts.nchan*1e9,
    freq=n.average(aa[0].beam.afreqs)*1e9,history=history,
    axes=('ra---sin','dec--sin','freq'))
