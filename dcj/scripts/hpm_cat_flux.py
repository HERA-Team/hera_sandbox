#!/usr/bin/env python
#
#  hpm_cat_flux.py
#  
#
#  Created by Danny Jacobs on 11/23/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,healpy as hpy, ephem
import pickle
"""
Find brightest flux near each input source.

"""
def pix2radec(nside,pix):
    """input a pixel index into an nside healpix vector
        return ra,dec (in radians)
    """
    try:
        hp_loc = [hpy.pix2ang(nside,px) for px in pix]
    except(TypeError):
        hp_loc = [hpy.pix2ang(nside,pix)]
    ra = [loc[1] for loc in hp_loc]
    dec = [n.pi/2 - loc[0] for loc in hp_loc]
    return n.array(ra),n.array(dec)
def radec2pix(nside,ra,dec):
    if n.abs(dec)>n.pi or ra>2*n.pi: raise ValueError('Input coordinates must be in radians')
    theta = n.pi/2 - dec
    phi = ra
    return hpy.ang2px(nside,theta,phi)
def pixsep(nside,pix1,pix2):
    """return the seperation between two healpix pixels in degrees"""
    src1_loc = hpy.pix2ang(nside,pix1)
    src2_loc = hpy.pix2ang(nside, pix2)
    aa = a.scripting.get_null_aa()
    src1 = a.phs.RadioFixedBody(src1_loc[1],n.pi/2 - src1_loc[0])
    src2 = a.phs.RadioFixedBody(src2_loc[1],n.pi/2 - src2_loc[0])
    src1.compute(aa)
    src2.compute(aa)
    return ephem.separation(src1,src2)    
def pos2name(pos):
    raname=''.join(str(ephem.hours(pos['RA'])).split(':')).split('.')[0]
    decname=''.join(str(ephem.degrees(pos['DEC'])).split(':')).split('.')[0]
    if n.sign(pos['DEC'])>0.: decname = '+'+decname
    return raname+decname
def R(x):
    return n.round(x,2)
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.add_option('-r',dest='radius',default=1.,type='float',
    help="Analysis region around each input source, in degrees. [1]")
o.add_option('-o',dest='outfile',
    help="Output the catalog in a 'standard format' [None]")
o.add_option('--hist',action='store_true',
    help="Output the histograms of pixels in the rms annulus (pickled dict)")
opts, args = o.parse_args(sys.argv[1:])
outer_radius = 4
inner_radius = 1
wcut = 90
clist,ccoff,ccats = a.scripting.parse_srcs(opts.src,opts.cat)
print clist,ccoff,ccats
if not opts.cal is None:
    print "cal"
    ccat = a.cal.get_catalog(opts.cal,srcs=clist,cutoff=ccoff,catalogs=ccats)
else:
    ccat = a.src.get_catalog(srcs=clist,cutoff=ccoff,catalogs=ccats)
aa = a.scripting.get_null_aa()
for src in ccat:
    ccat[src].compute(aa)
n.set_printoptions(precision=2,suppress=True)
for file in args:
    sky = hpy.read_map(file)
    nside = hpy.get_nside(sky)
    try: 
        Weight = hpy.read_map(file,field=1)
        #mask the lowest wcut% of the weights
#        WPDF,B = n.histogram(Weight,bins=n.sqrt(len(Weight)))
#        WCDF = n.cumsum(WPDF)/n.sum(WPDF)
#        wcutvali = n.argwhere(n.abs(WCDF-wcut/100.)==n.min(n.abs(WCDF-wcut/100.))).squeeze()
#        wcutval = B[wcutvali]
#        print "masking %3.1f of the data where weights are below %f"%(len(WCDF[WCDF<wcut/100.]),wcutval)
#        Weight = n.ma.masked_where(Weight<B[wcutval],Weight)
        sky /= Weight 
<<<<<<< HEAD
    except(IndexError): print "No weights found"
=======
    except(IndexError): 
        print "No weights found"
        Weight = n.ones_like(sky)
>>>>>>> b0ddfc5bcc04943dcd0263788e1031c580bc93f2
#    #get the analysis subregion
#    sky_an = n.zeros_like(sky)
#    print "getting the analysis subregion"
    #mask all the NaNs
    sky = n.ma.masked_where(n.isnan(sky),sky)
    okpx = []
    catpx = []
    print "extracting %d sources from %d catalog sources"%(len(okpx),len(ccat))
    head='\t'.join(('srcpx','_Ra','_Dec',
        'Name','S_nu_','nu','rms','DR',
<<<<<<< HEAD
        'catpx','dist','_Ra_cat','_Dec_cat','S_nu_cat_'))
=======
        'catpx','dist','_Ra_cat','_Dec_cat','S_nu_cat_','weight'))
>>>>>>> b0ddfc5bcc04943dcd0263788e1031c580bc93f2
    print head
    if not opts.outfile is None:
            outfile = open(opts.outfile,'w')
            outfile.write('#'+head+'\n')
    if opts.hist:
        histfile = open(opts.outfile[:-3]+'_hist.pkl','w')
    
    ccat = sorted(ccat.iteritems(),key=lambda src: src[1].jys,reverse=True)
#    for srcname,src in ccat.iteritems():
    foundsrcs = []
    hists = {}
    for srcname,src in ccat:
    
        catpx = hpy.ang2pix(nside,n.pi/2-src.dec,src.ra)
        v = hpy.pix2vec(nside,catpx)
        px = hpy.query_disc(nside,v,opts.radius)
#        px = n.array(list(decpx.intersection(px)))
        #print "pixels with nans",n.argwhere(n.logical_not(n.isnan(sky[px])) )
        px = px[n.argwhere(n.logical_not(n.isnan(sky[px])))]
        if len(px)==0: "no available points for %s"%(srcname,); continue
        srcpx = n.argwhere(sky==n.max(sky[px]))[0,0]
        if srcpx in foundsrcs: continue
        else: foundsrcs.append(srcpx)
        v = hpy.pix2vec(nside,srcpx)
        hole_pix = hpy.query_disc(nside,
                    v,inner_radius)
        outer_pix = set(hpy.query_disc(nside,
                    v,outer_radius))
        annulus_pix = n.array(list(outer_pix.difference(set(hole_pix))))
        if opts.hist:
            rmspix = sky[annulus_pix]
            rmspix = rmspix[n.logical_not(rmspix.mask)]
            srchist = n.histogram(rmspix,bins=int(n.sqrt(len(rmspix))))
            hists[srcname] = srchist
        pos = n.array(pix2radec(nside,srcpx)).squeeze()*a.img.rad2deg
#        pos *= a.img.rad2deg
#        pos.dtype = [('RA',n.float),('DEC',n.float)]
        line=(srcpx,R(pos[0]),R(pos[1]),
         srcname,R(sky[srcpx]),0.15,
        R(n.ma.std(sky[annulus_pix])),R(sky[srcpx]/n.ma.std(sky[annulus_pix])),
<<<<<<< HEAD
        catpx,R(pixsep(nside,srcpx,catpx)*a.img.rad2deg),R(src.ra*a.img.rad2deg),R(src.dec*a.img.rad2deg),R(src.jys.squeeze()))
=======
        catpx,R(pixsep(nside,srcpx,catpx)*a.img.rad2deg),R(src.ra*a.img.rad2deg),R(src.dec*a.img.rad2deg),R(src.jys.squeeze()),n.round(Weight[srcpx],4))
>>>>>>> b0ddfc5bcc04943dcd0263788e1031c580bc93f2
        print '\t'.join(map(str,line))
        if not opts.outfile is None:
            outfile.write('\t'.join(map(str,line))+'\n')
    if opts.hist:
        pickle.dump(hists,histfile)
        histfile.close()
if not opts.outfile is None: outfile.close()

        
        
