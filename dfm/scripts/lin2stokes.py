#!/usr/bin/env python

"""
This module adds two fits images together. It should go into some place in the structure of aipy or another.
"""

import aipy as a
import numpy as np
import optparse,sys

#DEBUG
from pylab import *

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
o.add_option('-f','--freq',dest='freq',default=0.15,type='float',help='Frequency at which to set the beam')
opts,args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal,0.001,opts.freq,1)

#Organize all the files
images = {}
pol_list = ['xx','xy','yx','yy']
stokes = ['I','Q','U','V']
for arg in args:
    fstr = arg.split('.')
    if fstr[-1] != 'fits': continue
    type = fstr[-2]
    if not type in images: images[type] = {}
    for pol in pol_list:
        if pol in arg:
            names = arg.split(pol)
            names = (names[0],names[1])
            if not names in images[type]: images[type][names] = []
            images[type][names].append(pol)

A = None
for type in images:
    for names in  images[type]:
        lin_images = {}
        for pol in pol_list:
            imfile = names[0]+pol+names[1]
            print 'Reading',imfile
            im,kwds = a.img.from_fits(imfile)
            lin_images[pol] = im

        #Calculate the beam
        if A is None:
            print 'Calculating beam weighting.'
            A = {}
            DIM = lin_images['xx'].shape[0]
            RES = 1. / np.abs(kwds['d_ra'] * a.img.deg2rad * DIM)
            im = a.img.Img(DIM*RES,RES,mf_order=0)
            
            tx,ty,tz = im.get_top(center=(DIM/2,DIM/2))
            tx = tx.flatten()
            ty = ty.flatten()
            tz = tz.flatten()
            
            for p in ('x','y'):
                A[p] = aa.ants[0].bm_response((tx,ty,tz),pol=p)
                A[p] = np.where( tz <= 0., 0., np.abs(A[p]))
                A[p] = A[p].reshape(im.shape)
                A[p] = A[p].data

        for S in stokes:
            Sname = names[0]+S+names[1]
            print 'Writing',Sname
            if S == 'I':
                if type == 'dbm': Sim = lin_images['xx']+lin_images['yy']
                else: 
                    Sim = A['y']*A['y']*lin_images['xx']
                    Sim += A['x']*A['x']*lin_images['yy']
                    Sim = np.where(A['x'] == 0, 0., Sim/(A['x']*A['y']))
            elif S == 'Q':
                if type == 'dbm': Sim = lin_images['xx']+lin_images['yy']
                else: 
                    Sim = A['y']*A['y']*lin_images['xx']
                    Sim -= A['x']*A['x']*lin_images['yy']
                    Sim = np.where(A['x'] == 0, 0., Sim/(A['x']*A['y']))
            elif S == 'U': 
                Sim = lin_images['xy'] + lin_images['yx']
            else: 
                if type == 'dbm': Sim = lin_images['xy']+lin_images['yx']
                else: Sim = lin_images['xy'] - lin_images['yx']

            a.img.to_fits(Sname, Sim, clobber=True, object=kwds['object'], obs_date=kwds['obs_date'],
                        ra=kwds['ra'], dec=kwds['dec'], epoch=kwds['epoch'],
                        d_ra=kwds['d_ra'], d_dec=kwds['d_dec'])
