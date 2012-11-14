#! /usr/bin/env python 
"""
Do the correct rotaiton between meaasured Q and U and intrinsic Q and U.
"""

import numpy as np
import aipy as a
import optparse,ephem,sys,os

o = optparse.OptionParser()
o.set_usage('QUrot.py zen.<JD>.[Q,U]0000.bim.fits')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])


Qfits = []
Ufits = []
for arg in args:
    if 'U_00' in arg: Ufits.append(arg)
    if 'Q_00' in arg: Qfits.append(arg)
Qfits.sort()
Ufits.sort()

for q_fit,u_fit in zip(Qfits,Ufits):
    qimg,kwds = a.img.from_fits(q_fit)
    uimg,kwds = a.img.from_fits(u_fit)
    s = ephem.Equatorial(kwds['ra']*a.img.deg2rad,kwds['dec']*a.img.deg2rad,epoch=ephem.J2000)
    s = ephem.Equatorial(s,epoch=kwds['obs_date'])
    LST,lat = s.get()
    print '='*30
    print 'Reading files %s %s' %(q_fit,u_fit)
    print 'LST,lat = ',LST,lat
    #Coordinates
    DIM = qimg.shape[-1]
    RES = 1. / (kwds['d_ra'] * a.img.deg2rad * DIM)
    im = a.img.Img(DIM*RES,RES,mf_order=0)
    ra,dec = a.coord.eq2radec(im.get_eq(LST,lat,center=(DIM/2,DIM/2)))

    #Calculating 
    X = a.pol.ParAng(LST-ra,dec,lat)
    sX,cX = np.sin(2.*X),np.cos(2.*X)

    Qint = cX*qimg - sX*uimg
    Uint = cX*uimg + sX*qimg
    Pint = np.ma.sqrt(Qint**2 + Uint**2)
    Xint = 0.5*np.arctan2(Uint,Qint)
    #Housekeeping
    names = q_fit.split('Q')
    newQname = names[0]+'Qrot'+names[1]
    newUname = names[0]+'Urot'+names[1]
    newPname = names[0]+'P'+names[1]
    newXname = names[0]+'X'+names[1]
    #write files
    a.img.to_fits(newQname, Qint.data, clobber=True, object=kwds['object'], obs_date=kwds['obs_date'],
                ra=kwds['ra'], dec=kwds['dec'], epoch=kwds['epoch'],
                d_ra=kwds['d_ra'], d_dec=kwds['d_dec'])
    a.img.to_fits(newUname, Uint.data, clobber=True, object=kwds['object'], obs_date=kwds['obs_date'],
                ra=kwds['ra'], dec=kwds['dec'], epoch=kwds['epoch'],
                d_ra=kwds['d_ra'], d_dec=kwds['d_dec'])
    a.img.to_fits(newPname, Pint.data, clobber=True, object=kwds['object'], obs_date=kwds['obs_date'],
                ra=kwds['ra'], dec=kwds['dec'], epoch=kwds['epoch'],
                d_ra=kwds['d_ra'], d_dec=kwds['d_dec'])
    a.img.to_fits(newXname, Xint.data, clobber=True, object=kwds['object'], obs_date=kwds['obs_date'],
                ra=kwds['ra'], dec=kwds['dec'], epoch=kwds['epoch'],
                d_ra=kwds['d_ra'], d_dec=kwds['d_dec'])
