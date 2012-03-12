#!/usr/bin/env python
#
#  cat_twopoint.py
#  
#
#  Created by Danny Jacobs on 2/10/10.
#  PAPER Project
#
#from matplotlib import use; use('agg')
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,time,ephem,healpy as hp
from pylab import *
o = optparse.OptionParser()
a.scripting.add_standard_options(o, src=True)
o.add_option('--fig',dest='fig',default=65,type='int',
    help='figure number')
o.add_option('--mbuff',dest='mbuff',default=1e4,type='float',
    help='Amount to buffer before doing a histogram.')
o.add_option('--no_cache',dest='no_cache',
   help = """Dont try and load cached histogram from a file for listed catalogs.  
   Force recalculation. eg  --no_cache=paper,vlss""")
o.add_option('--area',dest='area',default='coverage',
    help="""specify the area.  options are center_rad using 
    <src spec>-<degrees> or coverage map. [coverage]  """)
opts, args = o.parse_args(sys.argv[1:])

src_spec = a.scripting.parse_srcs(opts.src,opts.cat)
cat_names = src_spec[2]
srcs = src_spec[0]
flux_cut = src_spec[1]
def shift_bins(bins):
    return bins[:-1] + n.diff(bins)[0]

def update_pos(c):
    date=ephem.J2000
    for s in c.keys():
        try: ephem.FixedBody.compute(c[s], date)
        except(TypeError):
            if opts.juldate is None: del(c[s])
            else: ephem.Body.compute(c[s], date)
if not opts.no_cache is None: no_cache = opts.no_cache.split(',')
else:no_cache = []
def get_coverage(name):
    "returns a healpixmap of catalog coverage"
    coveragecachename = name+'_coverage.fits'
    try: 
        h = hp.read_map(coveragecachename)
        print "loading", coveragecachename
    except(IOError): 
        print "calculating coverage"
        if not opts.nside is None:
            h = n.zeros(hp.nside2npix(opts.nside))
            nside = opts.nside
        else: 
            nside = hp.npix2nside(hpm_lengths[n.searchsorted(hpm_lengths,
                    len(flux))])
            h = n.zeros(hp.nside2npix(nside))
        for i,s in enumerate(pos):
            if flux[i]>med: h[hp.ang2pix(nside,n.pi/2-s[1],s[0])] += 1
        #h = hp.smoothing(h,10,degree=True)
        hp.write_map(coveragecachename,h)
    print "saving coverage map to ",coveragecachename
    return h
if opts.area=='coverage':
    def rand_src(coverage):
        "return a source (RadioFixedBody) randomly drawn from hpm coverage pdf"
        global obs
        nside = hp.npix2nside(len(coverage))
        good_pix = n.where(coverage>1)[0]
        px = good_pix[int(n.random.uniform(0,len(good_pix)))]
        ang = hp.pix2ang(nside,px)
        src = a.phs.RadioFixedBody(ang[1],n.pi/2 - ang[0])
        src.compute(obs)
        return src
else:
    def rand_src(area):
        "return a random point from within the area specified by <cen>_<rad>"
        global obs
        rad = float(area.split('-')[1])
        cen = map(float,area.split('-')[0].split('_'))
        r = n.random.uniform(0,rad)*a.img.deg2rad
        t = n.random.uniform(0,2*n.pi)
        ra = cen[0]*n.pi/12 + r * cos(t)
        dec = cen[1]*a.img.deg2rad + r * sin(t)
        src = a.phs.RadioFixedBody(ra,dec)
        src.compute(obs)
        return src
        
obs = ephem.Observer()    
fig = figure(opts.fig)
for name in cat_names:
    cachename = name+'_tps.txt'
    if opts.area=='coverage':coverage = get_coverage(name)
    try:
        if not name in no_cache and not ('all' in no_cache): 
            tps = n.loadtxt(cachename)
            dt,w = tps[:,0],tps[:,1]
            print "loading cache" + cachename
    except(IOError):
        no_cache.append(name)
    if name in no_cache or 'all' in no_cache:
        cat = a.src.get_catalog(srcs=srcs,cutoff=flux_cut,catalogs=[name])
        if len(cat)==0: cat=a.src.get_catalog(catalogs=[name])
        if len(cat)==0: continue
        update_pos(cat)
        print name,len(cat)
        DDseps = []
        DRseps = []

        nbins = len(cat)
        DD = n.zeros(nbins)
        DR = n.zeros(nbins)
        t2=0
        DDs = []
        nseps = len(cat)/2.*(len(cat)-1)
        print time.strftime("%H:%M:%S")
        print "finding %d correlations"%(nseps)
        t1 = time.time()
        bcount = 0
        for i in range(nseps):
            try:
                srcA = cat.keys()[int(n.random.uniform(0,len(cat)))]
                srcB = cat.keys()[int(n.random.uniform(0,len(cat)))]
                if opts.area=='coverage': srcR = rand_src(coverage)
                elif len(opts.area.split('-'))>1: srcR = rand_src(opts.area)
        #        for srcA in cat:
        #            for srcB in cat:
                if srcA!=srcB and not (srcA,srcB) in DDs: 
                    DDseps.append(ephem.separation(cat[srcA], cat[srcB]))
                    DDs.append((srcA,srcB))
                    DRseps.append(ephem.separation(cat[srcA],srcR))
    #                if opts.area=='coverage':
    #                    DRseps.append(ephem.separation(cat[srcA],rand_src(coverage)))
    #                elif len(opts.area.split('-'))>1: 
    #                    DRseps.append(ephem.separation(cat[srcA],rand_src(opts.area)))
                if len(DRseps)==opts.mbuff:
                    print n.average(DDseps),n.average(DRseps),
                    dDD,sepbins = n.histogram(DDseps,bins=nbins)
                    dDR,sepbins = n.histogram(DRseps,bins=nbins)
                    sys.stdout.flush()
                    DR = DR+dDR
                    DD =DD+dDD
                    cnt =0
    #                del(DDseps)
    #                del(DRseps)
                    DDseps = []
                    DRseps = []
                    if not mod(bcount,10):
                        DD = (DD + dDD)/2.
                        DR = (DR + dDR)/2.
                        dt = shift_bins(sepbins)
                        dt *= a.img.rad2deg
                        w = DD/DR - 1
                        clf()
                        plot(dt,w,label=name)
                        savefig(name+'_twopoint_'+str(bcount)+'.png')
                        print "saving figure",name+'_twopoint_'+str(bcount)+'.png'
                    if t2==0: 
                        t2 = time.time()
                        print "single buffer executed in %d seconds"%(t2-t1)
                        print "est execution time: %d hours or %d minutes"%\
                        ((t2-t1)*nseps/opts.mbuff/3600.,(t2-t1)*nseps/opts.mbuff/60.)
                        print '.',
                    else:
                        t2 = time.time()
                        print t2-t1,'.',            
                    t1 = time.time()
                    bcount +=1
            except(KeyboardInterrupt):
                    break
                    
                    
        t2 = time.time()
#        print "found %d seperations in %3f s" %(len(DDseps),(t2-t1))
        print "calculating histogram"
        print len(DDseps)
        if cnt: 
            DD = (DD + dDD)/2.
            DR = (DR + dDR)/2.
        
        dt = shift_bins(sepbins)
        dt *= a.img.rad2deg
        print "saving to "+cachename
        w = DD/DR - 1
        n.savetxt(cachename,n.vstack([dt,w]).transpose())
    plot(dt,w,label=name)
legend()    
show()
                