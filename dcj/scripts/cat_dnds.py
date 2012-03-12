#!/usr/bin/env python
#
#  cat_dnds.py
#  
#
#  Created by Danny Jacobs on 1/16/10.
#  PAPER Project
#

import aipy as a, numpy as n,math as m
import sys, optparse,healpy as hp,ephem
from pylab import *
"""
Computes and plots flux distribution dN/ds of a catalog.  All fluxes 
are scaled to --freq using a single global power law defined by:
F_scaled = F*(--freq/mfreq)^(--global_index)

"""


o = optparse.OptionParser()
a.scripting.add_standard_options(o, src=True)
o.add_option('--global_index', dest='global_index', default=-1,type='float',
    help='Global index used to scale flux distributions to --freq. [-1]')
o.add_option('--freq',dest='freq',default=0.15,type='float',
   help = 'Frequency to which to scale fluxes. [0.15] GHz')
o.add_option('--no_cache',dest='no_cache',
   help = """Dont try and load cached histogram from a file for listed catalogs.  
   Force recalculation. eg  --no_cache=paper,vlss""")
o.add_option('--nside',dest='nside',type='int',
    help="Estimate the survey area on a healpix grid of this nside.[auto]")
o.add_option('--fig',dest='fig',default=65,
    help='figure number')
opts, args = o.parse_args(sys.argv[1:])

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


src_spec = a.scripting.parse_srcs(opts.src,opts.cat)
cat_names = src_spec[2]
srcs = src_spec[0]
fluxes = {}
#load all given catalogs

hpm_lengths = n.array([hp.nside2npix(i) for i in 2**n.array(range(1,11))])
stack = []
for name in cat_names:
    cachename = name+'_dNds.txt'
    coveragecachename = name+'_coverage.fits'    
    try:
        if not name in no_cache and not ('all' in no_cache): 
            dNds = n.loadtxt(cachename)
            ds,dNds = dNds[:,0],dNds[:,1]
            print "loading cache" + cachename
    except(IOError):
        no_cache.append(name)
    if name in ['three_cr']: 
        area = 95*n.pi/180*2*n.pi
    elif name in ['paper_deep']:
#            area = 4*n.pi-10*n.pi/180*2*n.pi
        area = (30*n.pi/180)**2
    elif name=='nvss':
        area = (90.+40.)*180*(n.pi/180)**2
    else:
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
        #area = n.average(h)*4*n.pi/len(h)
        area = len(n.where(h>0.9)[0])*4*n.pi/len(h)
    print name," survey area = ",area, " sr"        
    if name in no_cache or 'all' in no_cache:
        print "calculating dN/dsdO for "+name
        cat = a.src.get_catalog(srcs=srcs,catalogs=[name])
        if len(cat)==0: cat=a.src.get_catalog(catalogs=[name])
        if len(cat)==0: continue
        for src in cat:
            cat[src].index = opts.global_index
        cat.update_jys(opts.freq)
        update_pos(cat)
        print name,len(cat)
        flux = [cat[src].jys for src in cat]
        med = n.median(flux)
        pos = [(cat[src].ra,cat[src].dec) for src in cat]
        dN,bins = n.histogram(n.log10(flux),bins=n.max([n.sqrt(len(flux)),10]),
            new=True)
        ds = shift_bins(bins)
        dNds = dN/area/n.diff(10**bins)
        print "saving to "+ cachename
        n.savetxt(cachename,n.vstack([ds,dNds]).transpose())
    try:cat
    except(NameError):cat = a.src.get_catalog(srcs=srcs,catalogs=[name])
    stack.append((cat[cat.keys()[0]].mfreq,ds,dNds,name,area))
    del cat

fsort = n.argsort([s[0] for s in stack])
stack = [stack[f] for f in fsort]
freqs = [s[0] for s in stack]
#axes(axisbg='black')
fig1 = figure(opts.fig)
clf()
fig2 = figure(opts.fig+1)
clf()
for i,s in enumerate(stack):
    ds = s[1]
    dNds = s[2]
    f = s[0]
    name = s[3]
    area = s[4]
    c= 1-n.abs(f-0.15)/(n.max(freqs)-0.15)
    if c == 1.0: c=0.8
    c=0
    bins = n.concatenate([ds,[ds[-1]+n.diff(ds)[-1]]])
    dN = dNds*n.diff(bins)*area
    Nlts = n.cumsum(dN[range(len(dN)-1,-1,-1)])
    Nlts = Nlts[range(len(dN)-1,-1,-1)]




    if name=='paper': 
        figure(fig1.number)
        semilogy(ds,dNds,'k',label=name +' '+ str(f)+'GHz',lw=3)
        figure(fig2.number)
        semilogy(ds,Nlts,'k-',label=name +' '+ str(f)+'GHz',lw=3)
    elif name=='paper_deep':
        figure(fig1.number)
        semilogy(ds,dNds,'k--',label=name +' '+ str(f)+'GHz',lw=3)
        figure(fig2.number)
        semilogy(ds,Nlts,'k--',label=name +' '+ str(f)+'GHz',lw=3)
    else: 
        figure(fig1.number)
        semilogy(ds,dNds,label=name+' '+str(f)+'GHz')
        figure(fig2.number)
        semilogy(ds,Nlts,'-',label=name +' '+ str(f)+'GHz')
    

    




fig1 = figure(opts.fig)
legend(loc='upper right')
xlabel("log(flux [Jy])")
ylabel("dN/dsd$\Omega$ $Jy^{-1} sr^{-1}$")
fig2 = figure(opts.fig+1)
legend(loc='upper right')
ylabel("$N(<s)/d\Omega [sr^{-1}]$")
xlabel("log(flux [Jy])")
show()
    
    


    
    

