#! /usr/bin/env python
"""
Generate a Data cube of Q+iU, Perform Rotation Measure Synthesis on it.
"""

import numpy as np
import aipy as a
import os,sys,pyfits
import RotMeasTools as RMT
from glob import glob
from pylab import *
import time; t0 = time.time()

f0 = 0.129354
df = 0.00189
Nchan=32

fq = np.arange(f0,f0+(Nchan-1.)*df,df)

#Read in, sort fits files
DIR='synth_maps'
Ifiles = glob(DIR+'/RotSyn_FB*I0000_1:00_-30.fits')
Qfiles = glob(DIR+'/RotSyn_FB*Qrot0000_1:00_-30.fits')
Ufiles = glob(DIR+'/RotSyn_FB*Urot0000_1:00_-30.fits')
Vfiles = glob(DIR+'/RotSyn_FB*V0000_1:00_-30.fits')

I,Q,U,V = {},{},{},{}
for i,q,u,v in zip(Ifiles,Qfiles,Ufiles,Vfiles):
    I[int((i.split('FB')[1]).split('_')[0])] = i
    Q[int((q.split('FB')[1]).split('_')[0])] = q
    U[int((u.split('FB')[1]).split('_')[0])] = u
    V[int((v.split('FB')[1]).split('_')[0])] = v

bands = Q.keys()
bands.sort()

#Generate the data cube (fq)
print "Generating Data Cube."
Ic,Qc,Uc,Vc = None,None,None,None
for b in bands:
    i,kwds = a.img.from_fits(I[b])
    q,kwds = a.img.from_fits(Q[b])
    u,kwds = a.img.from_fits(U[b])
    v,kwds = a.img.from_fits(V[b])
    if Ic is None: 
        Ic = [i]
        Qc = [q]
        Uc = [u]
        Vc = [v]
    else: 
        Ic.append(i)
        Qc.append(q)
        Uc.append(u)
        Vc.append(v)
Ic = np.array(Ic)
Qc = np.array(Qc)
Uc = np.array(Uc)
Vc = np.array(Vc)

(Nfreq,Nra,Ndec) = Ic.shape
t1 = time.time()
print "\t --- Done in %3.3f s" %(t1-t0)

print "Writing to FITS"
a.img.to_fits('Fq_cube_I_1:00_-30.fits',Ic, clobber=True,
                axes=('ra--sin','dec--sin','freq'), 
                ra=kwds['ra'], dec=kwds['dec'], d_ra=kwds['d_ra'],d_dec=kwds['d_dec'],
                freq=fq[Nchan/2], d_freq=fq[1]-fq[0] )
a.img.to_fits('Fq_cube_Q_1:00_-30.fits',Qc, clobber=True,
                axes=('ra--sin','dec--sin','freq'), 
                ra=kwds['ra'], dec=kwds['dec'], d_ra=kwds['d_ra'],d_dec=kwds['d_dec'],
                freq=fq[Nchan/2], d_freq=fq[1]-fq[0] )
a.img.to_fits('Fq_cube_P_1:00_-30.fits',np.sqrt(Qc**2 + Uc**2), clobber=True,
                axes=('ra--sin','dec--sin','freq'), 
                ra=kwds['ra'], dec=kwds['dec'], d_ra=kwds['d_ra'],d_dec=kwds['d_dec'],
                freq=fq[Nchan/2], d_freq=fq[1]-fq[0] )
a.img.to_fits('Fq_cube_U_1:00_-30.fits',Uc, clobber=True,
                axes=('ra--sin','dec--sin','freq'), 
                ra=kwds['ra'], dec=kwds['dec'], d_ra=kwds['d_ra'],d_dec=kwds['d_dec'],
                freq=fq[Nchan/2], d_freq=fq[1]-fq[0] )
a.img.to_fits('Fq_cube_V_1:00_-30.fits',Vc, clobber=True,
                axes=('ra--sin','dec--sin','freq'), 
                ra=kwds['ra'], dec=kwds['dec'], d_ra=kwds['d_ra'],d_dec=kwds['d_dec'],
                freq=fq[Nchan/2], d_freq=fq[1]-fq[0] )
