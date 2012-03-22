#!/usr/bin/env python
#
#  im_src_find.py
#  
#
#  Created by Danny Jacobs on 2/3/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse, healpy as hpy,ephem


def find_src_in_cat(src,cat,sep):
    """
    Given a src (src = aipy.amp.RadioFixedObject or lower) 
    and a catalog (cat = aipy.amp.RadioCatalog or lower), find the closest
    source in cat that falls within sep [deg].
    Return: closest source object and seperation (ephem.Angle).
    If no source found within sep, return None,closest_sep
    """
    last_sep=n.pi
    closest_src=None
    for nk,test_src in cat.iteritems():
        cursep = ephem.separation(src, test_src)
        if cursep <= last_sep:
            last_sep = cursep
            closest_src = test_src
    if last_sep<sep * a.img.deg2rad: return closest_src,last_sep
    else: return None,last_sep
def update_pos(c):
    date=ephem.J2000
    for s in c.keys():
        try: ephem.FixedBody.compute(c[s], date)
        except(TypeError):
            if opts.juldate is None: del(c[s])
            else: ephem.Body.compute(c[s], date)
def pos2name(pos):
    raname=''.join(str(ephem.hours(pos['RA'])).split(':')).split('.')[0]
    decname=''.join(str(ephem.degrees(pos['DEC'])).split(':')).split('.')[0]
    if n.sign(pos['DEC'])>0.: decname = '+'+decname
    return raname+decname
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
def im_px2world(proj,px):
    """
    Convert a linear pixel number into a world coordinate.
    """
    npix = proj.naxis[0]*proj.naxis[1]
    if px>npix or px<0: return None
    grid = n.zeros(proj.naxis)
    grid = n.reshape(grid,npix)
    return n.argwhere(grid).squeeze()
    
o = optparse.OptionParser()
a.scripting.add_standard_options(o)
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
o.add_option('-o',dest='outfile',
    help="Output the catalog in a 'standard format' [None]")
o.add_option('--flux',dest='lower_flux',default=10,type='float',
    help="Lower flux limit in Jys. [10]")
o.add_option('-x','--exclude', dest='exclude', 
    help='A string that is parsed similarly to "srcs", but excludes sources that may otherwise have been selected by the "srcs" parameter.')
o.add_option('--xcat',dest='xcat',default='helm',
    help='A catalog in which to look up the excluded source. default = helm')
o.add_option('--gal_cut',dest='gal_cut',default=10,type='float',
    help='Galactic latitude to exclude')
o.add_option('-c','--center',dest='center',default='0_0',
    help='Center of region to analyze. Use format as for -s [0_0]')
o.add_option('--ccat',dest='ccat',default='helm,misc',
    help='Catalog selection for analysis region center [helm,misc]')
o.add_option('-r',dest='radius',default=180,type='float',
    help="Analysis region in degrees. [180]")
opts, args = o.parse_args(sys.argv[1:])
print opts.exclude
if opts.exclude != None:
    xlist,xoff,cats = a.scripting.parse_srcs(opts.exclude,opts.xcat)
    xcat = a.src.get_catalog(srcs=xlist, cutoff=xoff,catalogs=cats)
    update_pos(xcat)
else: xcat = {}
clist,ccoff,ccats = a.scripting.parse_srcs(opts.center,opts.ccat)
ccat = a.src.get_catalog(srcs=clist,cutoff=ccoff,catalogs=ccats)
update_pos(ccat)
for file in args:
    #load the data
#    sky = hpy.read_map(file)
#    nside = hpy.get_nside(sky)
#    try: 
#        wgts = hpy.read_map(file,field=1)
#        sky /= wgts #parsons includes data weights
#    except(IndexError): print "No weights found"
#    #get the analysis subregion
#    sky_an = n.zeros_like(sky)
#    for src in clist:
#        v = hpy.pix2vec(nside,hpy.ang2pix(nside,n.pi/2-src.dec,src.ra))
#        px = hpy.query_disc(nside,v,opts.radius)
#        sky_an[px] = sky[px]
#    sky = sky_an
#    del(sky_an)
#        
    hdulist = pf.open(file)
    proj = wcs.Projection(hdulist[0].header)
    sky = hdulist[0].data
    sky = n.reshape(sky,proj.naxis[0]*proj.naxis[1])
    #start the source selection algorithm
        #pair the healpix number and flux value
    sky_npx = n.array(zip(range(len(sky)),sky),dtype=[('npix',n.int),('flux',n.float)])
    print "sorting pixels"
    sky_npx_srt = n.sort(sky_npx,order='flux')
    #pick the top 1%. Should be well above my fluxlimit
    cut_flux = n.flipud(sky_npx_srt[-1*len(sky)/100:])
    cut_flux_pos = []
    #the standard dtype for the source array
    src_type =[('npix',n.int),('flux',n.float),('int_flux',n.float),('RA',n.float),
        ('DEC',n.float)]
    #get the locations of the sources, make sure to convert between 
    #ange and declination dec = pi/2-phi
    #also compute 
    for i,pix in enumerate(cut_flux['npix']):
#        pos = hpy.pix2ang(nside,pix)
        pos = im_px2world(proj,pix)
        cut_flux_pos.append((pix,cut_flux['flux'][i],
#            n.sum(sky[hpy.get_all_neighbours(nside,pix)])+sky[pix],
            cut_flux['flux'],
            pos[0]*n.pi/180,pos[1]*n.pi/180))
    srcs = n.array(list(cut_flux_pos),dtype=src_type)
    #sort by flux (sort is increasing, so flip it)
    srcs_fluxsort = n.flipud(n.sort(srcs,order='flux'))
    #use the flux sorted list to get those above my flux limit
    print "all sources with flux > ",opts.lower_flux, " Jy"
    print "all sources must also have peak flux>m*int_flux"
    srcs_fluxlim = []
    for t in srcs_fluxsort:
        if t['flux']<opts.lower_flux: break
        srcs_fluxlim.append(t)
#        srcline = "%i \t %s \t %s \t %s \t  %2.2f \n" % (
#            t['npix'],
#            str(ephem.hours(t['RA'])),
#            str(ephem.degrees(t['DEC'])),
#            pos2name(t),
#            t['flux'])
#        print srcline,

    srcs_fluxlim = n.array(srcs_fluxlim,dtype=src_type)

    #remove all dimmer sources within 0.5 degrees of the brightest source
    #NB this method should be improved to make this distance flux dependent
    peak_select_distance = 0.5 #select peak flux within this 
    remaining = len(srcs_fluxlim)-1
    remaining_srcs = n.flipud(srcs_fluxlim.copy())
    found_srcs = []
    while remaining:
        nearby_pix = hpy.query_disc(nside,
            hpy.pix2vec(nside,
            remaining_srcs[remaining]['npix']),peak_select_distance)
        nearby_srcs = n.array([s for s in remaining_srcs if s['npix'] in nearby_pix],
            dtype=src_type)

        cur_src = nearby_srcs[-1]

        found_srcs.append(nearby_srcs[-1])
        nearby_srcs = n.delete(nearby_srcs,-1)
        if len(nearby_srcs): 
            print "cleaning around a %2.1f Jy source pix=%i "%(cur_src['flux'],cur_src['npix'])
            print "found ",len(nearby_srcs)," nearby source(s)"
            remaining-=1
            deleted_flux = []
            for s in nearby_srcs:
                #print "deleting source: ",s
                dist = pixsep(nside,cur_src['npix'],s['npix'])/a.const.arcsec #arcsecond distance
                #if s['flux']<(0.48/dist * cur_src['flux']):
                if True:
                    print "deleting a %2.1f Jy source %2.2f \' away" %(
                        s['flux'],dist/60)
                    remaining_srcs = n.delete(remaining_srcs,
                        n.where(s['npix']==remaining_srcs['npix'])[0],axis=0)
                    deleted_flux.append(s['flux'])
                remaining-=1
            print  "deleted ",len(deleted_flux),"sources and a total of ",n.sum(deleted_flux), "Jys"
        else:remaining-=1
    found_srcs = n.array(found_srcs,
            dtype=src_type)
    #remove "sources" in the galaxy
    prov_srcs = []
    if not opts.gal_cut is None:
        obs = ephem.Observer()
        for s in found_srcs:
            glat = ephem.Galactic(ephem.Equatorial(s['RA'],s['DEC'])).lat
            if n.abs(glat)>opts.gal_cut*a.img.deg2rad: prov_srcs.append(s)
    print "deleted",len(found_srcs)-len(prov_srcs)," in galactic cut"
    found_srcs = n.array(prov_srcs,dtype=src_type)

    #remove sources in the "exclude" list
    final_srcs = []
    if not opts.exclude is None:
        obs = ephem.Observer()
        for s in found_srcs:
            psrc = a.phs.RadioFixedBody(s['RA'],s['DEC'])
            psrc.compute(obs)
            xsrc,dist = find_src_in_cat(psrc,xcat,peak_select_distance)
            if not xsrc is None and dist<peak_select_distance*n.pi/180:
                 xsrc.update_jys(0.15)
                 s['flux'] = xsrc.jys
                 s['int_flux']=xsrc.jys
                 print "using catalog flux for:",xsrc.src_name
            #else: final_srcs.append(s)
            final_srcs.append(s)
        final_srcs = n.array(final_srcs,dtype=src_type)
#        print "excluding ", len(found_srcs)-len(final_srcs)," as requested"
    else:final_srcs=found_srcs
    
    #scan for sources with obvious cleaning 
    problem_srcs = []
    for i,s in enumerate(final_srcs):
        if s['int_flux']<0.9*s['flux']:
            problem_srcs.append(i)
    
    if not opts.outfile is None:
        file = open(opts.outfile,'w')
        file.write("#healpix pix \t RA \t DEC \t pk flux [Jy] \t int flux [Jy]\n")
        for i,t in enumerate(final_srcs):
            if not i in problem_srcs:
                srcline = "%10i \t %s \t %s \t %s \t  %2.2f \t %2.2f \n" % (
                    t['npix'],
                    str(ephem.hours(t['RA'])),
                    str(ephem.degrees(t['DEC'])),
                    pos2name(t),
                    t['flux'],
                    t['int_flux'])
                print srcline,
                file.write(srcline)
        file.close()
    print "output %i sources to catalog"%(len(final_srcs)-len(problem_srcs),)
    print "the following %i sources were flagged as being odd" %(len(problem_srcs),)
    for i in problem_srcs:
        t = final_srcs[i]
        srcline = "%10i \t %s \t %s \t %s \t  %02.2f \t %2.2f \n" % (
            t['npix'],
            str(ephem.hours(t['RA'])),
            str(ephem.degrees(t['DEC'])),
            pos2name(t),
            t['flux'],
            t['int_flux'])
        print srcline,
                  
            
        


   # srcs_decsort = n.sort(srcs,order='DEC')
#    for t in srcs_decsort:
#        print "%2.2f /t %s /t %s" % (t['flux'],str(ephem.hours(t['RA'])),
#        str(ephem.degrees(t['DEC'])))
#        print f,ephem.hours(pos[1]),ephem.degrees(pos[0])
        
