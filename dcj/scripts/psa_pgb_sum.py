#!/usr/bin/env python
#
#  psa_pgb_sum.py
#  
#
#  Created by Danny Jacobs on 9/7/10.
#  PAPER Project
#
"""
Add skymaps from northern and southern hemisphere. Weight each in declination by
the width of the primary beam

Produces a file that can be interpreted in the nprmal way (divides by weights).
"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,healpy as hp

o = optparse.OptionParser()
#a.scripting.add_standard_options(o, cal=True)
opts, args = o.parse_args(sys.argv[1:])


northlat = 38.5*a.img.deg2rad
southlat = -30.*a.img.deg2rad
beamwidth = 60*a.img.deg2rad #fwhm

print "loading: ",args[0]
southsky = hp.read_map(args[0])
southweights = hp.read_map(args[0],field=1)
NSIDE_s = hp.get_nside(southsky)
print "NSIDE = ",NSIDE_s

print "loading: ",args[1]
northsky = hp.read_map(args[1])
northweights = hp.read_map(args[1],field=1)
NSIDE_n = hp.get_nside(northsky)
print "NSIDE = ",NSIDE_n


if NSIDE_n > NSIDE_s:
    print "degrading northern map to match southern resolution:", NSIDE_s
    northsky = hp.ud_grade(northsky, NSIDE_s, pess=False, order_in='RING')
    northweights = hp.ud_grade(northweights, NSIDE_s, pess=False, order_in='RING')
elif NSIDE_n < NSIDE_s:
    print "degrading southern map to match northern resolution:", NSIDE_n
    southsky = hp.ud_grade(southsky, NSIDE_n, pess=False, order_in='RING')
    southweights = hp.ud_grade(southweights, NSIDE_n, pess=False, order_in='RING')
NSIDE = n.min([NSIDE_n,NSIDE_s])
inds = n.arange(hp.nside2npix(NSIDE))
theta,phi = hp.pix2ang(NSIDE, inds)
dec,ra = n.pi/2-theta,phi
#print inds,(northlat-dec)**2/beamwidth
#p.plot(dec)
#p.show()
southweights /= n.exp(-(southlat-dec)**2/beamwidth)
northweights /= n.exp(-(northlat-dec)**2/beamwidth)*25

sky = a.map.Map(nside=NSIDE,scheme='RING')
sky.map.map = northsky + southsky
sky.wgt.map = southweights + northweights
sfile = args[0].split('/')[-1]
nfile = args[1].split('/')[-1]
filename = sfile[:-len('.fits')] + '_'+nfile
print "writing: ",filename
sky.to_fits(filename,clobber=True)