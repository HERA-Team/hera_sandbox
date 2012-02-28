#!/usr/bin/env python
"""
A script for filtering using a delay/delay-rate transform.  If a source
is specified, will remove/extract that source.  If none is specified,
will filter/extract in absolute terms.
"""

import numpy as n, pylab as p
import os, sys, optparse, re

o = optparse.OptionParser()
o.set_usage('pk_npz.py [options] *.uv')
o.set_description(__doc__)
opts, args = o.parse_args(sys.argv[1:])

def late_eor_pspec(k):
    pk = n.where(k >= .1, 2*k**-3 / 1e6, 2*(.1)**-3/1e6 * (k/.1)**-1)
    return pk

def rebin_log(x, y, nbins=10):
    '''For y=f(x), bin x into logrithmic (base 10) bins, and average y over
    these bin sizes.'''
    logx = n.log10(n.abs(x))
    hist1,bins = n.histogram(logx, bins=nbins, weights=y)
    hist2,bins = n.histogram(logx, bins=nbins)
    logx = .5 * (bins[1:] + bins[:-1])
    return 10**logx, hist1 / n.where(hist2 == 0, 1., hist2)

VERSUS_K = True
PLOT_KCUBE = True

filename_regx = re.compile(r'pk_npz__(\d+)__(\d+.\d+).npz')
pspec = {}
for filename in args:
    print 'Loading', filename
    umag,fq = filename_regx.match(filename).groups()
    umag,fq = int(umag),float(fq)
    if not pspec.has_key(umag): pspec[umag] = {}
    if not pspec[umag].has_key(fq): pspec[umag][fq] = {}
    f = n.load(filename)
    pspec[umag][fq]['k_pl'] = f['k_pl']
    pspec[umag][fq]['k_pr'] = f['k_pr']
    pspec[umag][fq]['ks'] = f['ks']
    pspec[umag][fq]['_ks'] = f['_ks']
    pspec[umag][fq]['pspec'] = f['pspec']
    pspec[umag][fq]['_pspec'] = f['_pspec']

nplots = len(pspec)
d2 = n.ceil(nplots**.5)
d1 = n.ceil(nplots / d2)
d1,d2 = int(d1),int(d2)

umags = pspec.keys(); umags.sort()
print 'Generating plots'
for cnt, umag in enumerate(umags):
    p.subplot(d1,d2,cnt + 1)
    fqs = pspec[umag].keys(); fqs.sort()
    if VERSUS_K:
        for i,fq in enumerate(fqs):
            color = 'kbgrcmy'[i%7]
            symbol = ['-','--','-.',':'][(i/7)%4]
            _ks,_pspec = pspec[umag][fq]['_ks'], pspec[umag][fq]['_pspec']
            if PLOT_KCUBE:
                # For k^3/(2pi^2)P(k), remember we've already divided by (2pi)^3
                #p.loglog(_ks, 4*n.pi*_ks**3*1e6*n.abs(_pspec), color+symbol, label=str(fq))
                p.loglog(_ks, 4*n.pi*_ks**3*1e6*n.abs(_pspec/1e-7*1.4e-5), color+symbol, label=str(fq))
            else:
                p.loglog(_ks, 1e6*n.abs(_pspec), color+symbol, label=str(fq))
        if PLOT_KCUBE:
            p.loglog(_ks, 4*n.pi*_ks**3*1e6*late_eor_pspec(_ks), 'k-')
        else:
            p.loglog(_ks, 1e6*late_eor_pspec(_ks), 'k-')
        p.xlabel(r'$k (h\ {\rm Mpc})^{-1}$')
        p.xlim(1e-2,2e0)
        #p.ylim(1e0,1e10)
        p.ylabel(r'${\rm mK}^2$')
        p.grid()
    else:
        pspec_key,ks_key = '_pspec','_ks'
        data = n.array([pspec[umag][fq][pspec_key] for fq in fqs])
        ks = n.log10(n.array([pspec[umag][fq][ks_key] for fq in fqs]))
        fqs = n.array([fq*n.ones_like(pspec[umag][fq][ks_key]) for fq in fqs])
        #p.imshow(n.log10(1e6*n.abs(data)), vmax=6, vmin=0, aspect='auto', interpolation='nearest')
        p.contourf(ks, fqs, n.log10(1e6*n.abs(data)), n.arange(-1,5,.1))
        p.ylim(fqs[0,0], fqs[-1,0])
        p.xlim(-2,0.5)
        p.xlabel(r'${\rm log}_{10}[k (h\ {\rm Mpc})^{-1}]$')
        p.ylabel(r'$\nu\ {\rm GHz}$')
    p.title(r'$|\vec u|=%d$' % (umag))

p.show()
