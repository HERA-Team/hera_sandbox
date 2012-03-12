#!/usr/bin/env python
#
#  extract_spectrum.py
#  
#
#  Created by Danny Jacobs on 1/19/10.
#  PAPER Project
#
"""
given an input source and fits data cube, extract the spectrum of that source
"""
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,pyfits as pf,ephem
from kapteyn import wcs


o = optparse.OptionParser()
o.set_usage('extract_spectrum.py -s srcs facets*.fits')
o.set_description(__doc__)

a.scripting.add_standard_options(o, src=1)
o.add_option('--sep',dest='sep',default=5,type='float',
   help='Selection criteria for source distance from img center. [5 degrees]')
o.add_option('-v',dest='verb',action='store_true',
    help='Print more')
o.add_option('--prefix',dest='prefix',
   help="Thumbnail filename prefix [optional]")
o.add_option('--width',dest='width',default=10,type='float',
 help="Search within this square region for the brightest continuum source. [10 deg]")
o.add_option('--plot',dest='plot',
  help="Plot in either: normal, glyph. [optional]")
o.set_description(__doc__)
opts, args = o.parse_args(sys.argv[1:])
opts.sep *= a.img.deg2rad
def update_pos(c):
    date=ephem.J2000
    for s in c.keys():
        try: ephem.FixedBody.compute(c[s], date)
        except(TypeError):
            if opts.juldate is None: del(c[s])
            else: ephem.Body.compute(c[s], date)


srcs,coff,catalogs = a.scripting.parse_srcs(opts.src,opts.cat)
cat = a.src.get_catalog(srcs=srcs,catalogs=catalogs)
update_pos(cat)



for file in args:
    print file
    hdulist = pf.open(file)
    hdr = hdulist[0].header
    width = n.abs(opts.width/hdulist[0].header.get('CDELT1'))
    proj = wcs.Projection(hdulist[0].header)
    img_ra = hdulist[0].header.get('CRVAL1')*a.img.deg2rad
    img_dec = hdulist[0].header.get('CRVAL2')*a.img.deg2rad
    center = a.phs.RadioFixedBody(img_ra,
            img_dec)
    ephem.FixedBody.compute(center,ephem.J2000)
    for name,src in cat.iteritems():
        src_sep = ephem.separation(src,center)
        if src_sep<opts.sep:
            if opts.verb: print "getting \n",src
            ra = src.ra * a.img.rad2deg
            dec = src.dec * a.img.rad2deg
            if opts.verb:print src_sep,ra,dec
            px = n.round(proj.topixel((ra,dec,1,1)))
            if opts.verb:print "at",px
#            if opts.verb:print "subimg"
#            if opts.verb: print BL, TR

#            im = n.flipud(hdulist[0].data.squeeze()[sub[1]+1,sub[0]+1])
##            imshow(im,aspect='equal')
            BL = n.round((px[0]-width/2,px[1]-width/2))
            TR = n.round((px[0]+width/2,px[1]+width/2))
            sub = n.meshgrid(range(BL[0],TR[0]), range(BL[1],TR[1]))
            maxpow = 0     
            bpix = px       
            for crd in n.vstack((n.concatenate(sub[0]),n.concatenate(sub[1]))).transpose():
                if opts.verb: print crd,
                spec = n.array(hdulist[0].data.squeeze()[:,crd[1],crd[0]])
                spec = n.ma.array(spec,mask=n.isnan(spec))
                pow = n.ma.average(spec)
                if opts.verb: print pow,
                if pow>maxpow: 
                    bpix = crd
                    maxpow = pow
                if opts.verb: print maxpow
            spec = hdulist[0].data.squeeze()[:,bpix[1],bpix[0]]
            spec = n.ma.array(spec,mask=n.isnan(spec))
            bpix = tuple(bpix)
            bpix += (1,)
            bpix += (1,)
            print "given pixel coords:", px
            print "found better pixel coords:", bpix
            loc = proj.toworld(bpix)
            #save file
            freqs = n.linspace(hdr.get('CRVAL3'),
                    hdr.get('NAXIS3')*hdr.get('CDELT3')+hdr.get('CRVAL3'),
                    num=hdr.get('NAXIS3'))
            spectrum = n.vstack((freqs,spec)).transpose()
            outfile = name+"_spectrum.txt"
            if not opts.prefix is None: outfile = opts.prefix+outfile
            print ' > '+ outfile
            #file = open(outfile,'w')
            if loc[1]>0: sgn = '+'
            else: sgn = ''
            src = str(ephem.hours(loc[0]*a.img.deg2rad)) + sgn + \
                str(ephem.degrees(loc[1]*a.img.deg2rad)) + "\n"
            #file.write(src) 
            print src
            n.savetxt(outfile,spectrum)
            outfile = '.'.join(outfile.split('.')[:-1])+'.'+'png'
            print outfile
            if not opts.plot is None:
                fig = p.figure()
                p.loglog(freqs,spec,'.')
                p.xlabel('Hz')
                p.ylabel('flux [Jy]')
                p.title(src)
                fig.savefig(outfile)
                p.close()
    