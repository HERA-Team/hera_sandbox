#! /usr/bin/env python 

import aipy as a
import numpy as np
import sys,os,optparse

from pylab import *

o = optparse.OptionParser()
o.set_usage('gen_pol_spec.py [options] *.uv')
a.scripting.add_standard_options(o,ant=True,chan=True,cal=True,src=True)
o.add_option('--clean',dest='clean',type='float',default=1e-3,help='clean tolerance (1e-3 is recommended)')
o.add_option('--sig',dest='sig',type='int',default=-1,help='half-width of the LPF')
o.add_option('--altmin',dest='altmin',type='float',default=0,help='Minimum altitude for pointing, in degrees. When phase center is lower than this altitude, data is omitted.')
opts,args = o.parse_args(sys.argv[1:])

srcname = opts.src
uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal,uv['sdf'],uv['sfreq'],uv['nchan'])
chans = a.scripting.parse_chans(opts.chan,uv['nchan'])
freqs = np.arange(uv['sfreq'], uv['sfreq']+uv['nchan']*uv['sdf'], uv['sdf'])
freqs = freqs.take(chans)
aa.select_chans(chans)
srclist,cutoff,catalogs = a.scripting.parse_srcs(srcname,opts.cat)
cat = a.cal.get_catalog(opts.cal,srclist,cutoff,catalogs)
src = cat[srcname]
del(uv)

specs = {
         'I':np.zeros(len(chans),dtype=np.complex),
         'Q':np.zeros(len(chans),dtype=np.complex),
         'U':np.zeros(len(chans),dtype=np.complex),
         'V':np.zeros(len(chans),dtype=np.complex),
}
sumwgts = {
         'IQ':np.zeros(len(chans),dtype=np.complex),
         'UV':np.zeros(len(chans),dtype=np.complex),
        }
WGT = {}
parang = {}
for uvfile in args:
    print 'Reading',uvfile
    D = {}
    uv = a.miriad.UV(uvfile)
    a.scripting.uv_selector(uv,opts.ant)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        pol = a.miriad.pol2str[uv['pol']]
        d = d.take(chans)
        f = f.take(chans)
        if not t in D.keys():
            D[t] = {}
            aa.set_jultime(t)
            src.compute(aa)
            cat.compute(aa)
            aa.sim_cache(cat.get_crds('eq',ncrd=3))
            parang[t] = a.pol.ParAng(aa.sidereal_time()-src.ra,src.dec,aa.lat)
        if src.alt < opts.altmin * a.img.deg2rad: continue
        try:
            d = aa.phs2src(d,src,i,j,pol)
        except(a.phs.PointingError): continue
        pb = aa.passband(i,j,pol)
        wgt = aa.bm_response(i,j,pol).squeeze()
        d = np.where(f,0,d); wgt = np.where(f,0,wgt)
        if not pol in D[t].keys(): D[t][pol] = d*wgt/pb
        else: D[t][pol] += d*wgt/pb
        if not pol in WGT.keys(): WGT[pol] = wgt**2
        else: WGT[pol] += wgt**2

    for t in D:
        try:
            Q = D[t]['xx']-D[t]['yy']
            U = D[t]['xy']+D[t]['yx']
            specs['I'] += D[t]['xx']+D[t]['yy']
            specs['Q'] += Q*np.cos(2*parang[t]) - U*np.sin(2.*parang[t])
            specs['U'] += U*np.cos(2*parang[t]) + Q*np.sin(2.*parang[t])
            specs['V'] += 1.j*(D[t]['xy']-D[t]['yx'])
        except(KeyError): continue

specs['I'] /= (WGT['xx']+WGT['yy'])
specs['Q'] /= (WGT['xx']+WGT['yy'])
specs['U'] /= (WGT['xy']+WGT['yx'])
specs['V'] /= (WGT['xy']+WGT['yy'])

figure(0)
subplot(131)
title('I')
plot(freqs,np.abs(specs['I']))
subplot(132)
title('Q,U')
plot(freqs,(specs['Q']/specs['I']).real,label='Q')
plot(freqs,(specs['U']/specs['I']).real,label='U')
legend()
subplot(133)
title('V')
plot(freqs,specs['V'].real)
draw()

figure(1)
plot(freqs,np.sqrt(np.abs(specs['Q'])**2 + np.abs(specs['U'])**2)/specs['I'])
title('p')
draw()

sumwgts['IQ'] = WGT['xx']+WGT['yy']
sumwgts['UV'] = WGT['xy']+WGT['yx']

figure(2)
for pol in sumwgts:
    plot(freqs,sumwgts[pol],label=pol)
legend()
draw()

show()
