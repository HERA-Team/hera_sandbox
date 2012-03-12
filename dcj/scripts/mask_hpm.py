#!/usr/bin/env python
#
#  mask_hpm.py
#  
#
#  Created by Danny Jacobs on 1/28/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse, healpy as hp,ephem
"""
Apply a mask to a healpix map (hpm) by zeroing values.
Mask can be:
declination cut min,max in degrees
an existing hpm mask with 1s for mask
second weight column, expressed as percentage (weight value level that will mask % pixels)

Output either a masked version or a hpm of the mask rule.

"""
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
o.add_option('--mask',dest='type',
    help="type:value. [dec:0,-60],hpm:<filename.fits> or wgt:50 ")
o.add_option('--gen',dest='gen',
    help="Output mask file of this filename.  [optional]")

opts, args = o.parse_args(sys.argv[1:])
if opts.type is None:
    opts.type = "dec:0,-60"
for file in args:
    print file
    h = hp.read_map(file)
    npix = len(h)
    nside = hp.npix2nside(npix)
    mask = n.zeros_like(h)
    type = opts.type.split(':')[0]
    value = opts.type.split(':')[1]
    outfile = file.split('.')
    #print type,value
    print "calc mask"
    if type=="dec":
        lim = map(lambda x: float(x)*a.img.deg2rad,value.split(','))
        in_range = lambda x: x<lim[0] and x>lim[1]
        decs = n.array([ephem.Galactic(hp.pix2ang(nside,i)[1],n.pi/2 - hp.pix2ang(nside,i)[0]).to_radec()[1] for i in range(npix)])
        mask[n.where([in_range(dec) for dec in decs])[0]] = 1
    elif type=="hpm":
        if not opts.gen is None: sys.exit()
        else: mask = hp.read_map(value)
        if len(mask)!=ma: print "mask file dimension mismatch"; sys.exit()
    elif type=="wgt":
        wgts = hp.read_map(file,field=1)
        i = int(float(value)/100*npix)
        wgts = n.sort(wgts)
        wgt_min = wgts[i] 
        mask[n.argwhere(wgts>wgt_min)] = 1
    if not opts.gen is None:
        print opts.gen
        hp.write_map(opts.gen,mask)  
    else:
        if len(outfile)>2:
             outfile = '.'.join(outfile[:-1])+'m'+'.fits'
        else: outfile = '.'.join(outfile[:-1])+'.m.fits'
        print outfile
        hp.write_map(outfile,h*mask)

    
    