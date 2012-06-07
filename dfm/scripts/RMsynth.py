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
Qfiles = glob(DIR+'/RotSyn_FB*Qrot0000_1:00_-30.fits')
Ufiles = glob(DIR+'/RotSyn_FB*Urot0000_1:00_-30.fits')

Q,U = {},{}
for q,u in zip(Qfiles,Ufiles):
    Q[int((q.split('FB')[1]).split('_')[0])] = q
    U[int((u.split('FB')[1]).split('_')[0])] = u
bands = Q.keys()
bands.sort()

#Generate the data cube (fq)
print "Generating Data Cube."
Pcube = None
for i in bands:
    q,kwds = a.img.from_fits(Q[i])
    u,kwds = a.img.from_fits(U[i])
    P = q+1.j*u
    if Pcube is None: Pcube = [P]
    else: Pcube.append(P)
Pcube = np.array(Pcube)
(Nfreq,Nra,Ndec) = Pcube.shape
t1 = time.time()
print "\t --- Done in %3.3f s" %(t1-t0)

print "Writing to FITS"
a.img.to_fits('Fq_cube_Q_1:00_-30.fits',Pcube.real, clobber=True,
                axes=('ra--sin','dec--sin','freq'), 
                ra=kwds['ra'], dec=kwds['dec'], d_ra=kwds['d_ra'],d_dec=kwds['d_dec'],
                freq=fq[Nchan/2], d_freq=fq[1]-fq[0] )
a.img.to_fits('Fq_cube_U_1:00_-30.fits',Pcube.imag, clobber=True,
                axes=('ra--sin','dec--sin','freq'), 
                ra=kwds['ra'], dec=kwds['dec'], d_ra=kwds['d_ra'],d_dec=kwds['d_dec'],
                freq=fq[Nchan/2], d_freq=fq[1]-fq[0] )

#Do the RM transform
print "Beginning Rotation Measure Synthesis"
RMs,W = RMT.RMTmat(fq)

RMcube = np.zeros_like(Pcube)
for i in range(Nra):
    for j in range(Ndec):
        RMcube[:,i,j] = RMT.RMT(Pcube[:,i,j],W)
t2 = time.time()
print "\t --- Done in %3.3f s" %(t2-t1)

plot(RMs,RMcube[:,Nra/2,Ndec/2].real)
show()

print "Writing to FITS file"
a.img.to_fits('RM_cube_real_1:00_-30.fits',RMcube.real, clobber=True,
                axes=('ra--sin','dec--sin','freq'), 
                ra=kwds['ra'], dec=kwds['dec'], d_ra=kwds['d_ra'],d_dec=kwds['d_dec'],
                freq=0, d_freq=RMs[1]-RMs[0] )
a.img.to_fits('RM_cube_imag_1:00_-30.fits',RMcube.imag, clobber=True,
                axes=('ra--sin','dec--sin','freq'), 
                ra=kwds['ra'], dec=kwds['dec'], d_ra=kwds['d_ra'],d_dec=kwds['d_dec'],
                freq=0, d_freq=RMs[1]-RMs[0] )
















