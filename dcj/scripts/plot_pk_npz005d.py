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
o.add_option('-m','--mode',dest='mode', default='pk',
    help='Plotting mode: pk, k3pk, snr, pk_vs_fq, k3pk_vs_fq')
o.add_option('--name',type='str',
    help='Optional name to give figure window')
opts, args = o.parse_args(sys.argv[1:])

if opts.name is None:
    opts.name = ' '.join(sys.argv)

def eor_k3pk(k):
    return 100*n.where(k >= .1, 1, (k/.1)**2)
def eor_pk(k):
    return eor_k3pk(k) * 2*n.pi**2 / k**3

def rebin_log(x, y, nbins=10):
    '''For y=f(x), bin x into logrithmic (base 10) bins, and average y over
    these bin sizes.'''
    logx = n.log10(n.abs(x))
    hist1,bins = n.histogram(logx, bins=nbins, weights=y)
    hist2,bins = n.histogram(logx, bins=nbins)
    logx = .5 * (bins[1:] + bins[:-1])
    return 10**logx, hist1 / n.where(hist2 == 0, 1., hist2)

VERSUS_K = not opts.mode.endswith('vs_fq')
PLOT_KCUBE = opts.mode.startswith('k3pk')
SNR = opts.mode.startswith('snr')
AMP = 0.2

filename_regx = re.compile(r'.*pk_npz__(\d+)__(\d+.\d+).npz')
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

prefactor = 730*.76*3.5e6 / 8 * AMP
umags = pspec.keys(); umags.sort()
print 'Generating plots'
fig = p.figure()
for cnt, umag in enumerate(umags):
    p.subplot(d1,d2,cnt + 1)
    fqs = pspec[umag].keys(); fqs.sort()
    if VERSUS_K:
        for i,fq in enumerate(fqs):
            color = 'kbgrcmy'[i%7]
            symbol = ['-','--','-.',':'][(i/7)%4]
            _ks,_pspec = pspec[umag][fq]['_ks'], pspec[umag][fq]['_pspec']
            if SNR:
                p.loglog(_ks, eor_pk(_ks)/prefactor/(1e6*n.abs(_pspec)), color+symbol, label=str(fq))
            elif PLOT_KCUBE:
                p.loglog(_ks, _ks**3/(2*n.pi**2)*prefactor*1e6*n.abs(_pspec), color+symbol, label=str(fq))
            else:
                p.loglog(_ks, 1e6*n.abs(_pspec), color+symbol, label=str(fq))
        if not SNR:
            if PLOT_KCUBE:
                p.loglog(_ks, eor_k3pk(_ks), 'k-')
            else:
                p.loglog(_ks, eor_pk(_ks)/prefactor, 'k-')
        p.xlabel(r'$k (h\ {\rm Mpc})^{-1}$')
        p.xlim(1e-2,2e0)
        if SNR:
            p.ylabel('SNR')
            p.ylim(1e-7,1e0)
        elif PLOT_KCUBE:
            p.ylabel(r'$\Delta^2(k)$')
            p.ylim(1e0,1e10)
        else:
            p.ylabel(r'$T_{\rm rms}^2 [{\rm mK}^2]$')
            p.ylim(1e-6,1e3)
        p.grid()
    else:
        pspec_key,ks_key = '_pspec','_ks'
        data = n.array([pspec[umag][fq][pspec_key] for fq in fqs])*AMP
        ks = n.log10(n.array([pspec[umag][fq][ks_key] for fq in fqs]))
        fqs = n.array([fq*n.ones_like(pspec[umag][fq][ks_key]) for fq in fqs])
        #p.imshow(n.log10(1e6*n.abs(data)), vmax=6, vmin=0, aspect='auto', interpolation='nearest')
        p.contourf(ks, fqs, n.log10(1e6*n.abs(data)), n.arange(-1,5,.1))
        p.ylim(fqs[0,0], fqs[-1,0])
        p.xlim(-2,0.5)
        p.xlabel(r'${\rm log}_{10}[k (h\ {\rm Mpc})^{-1}]$')
        p.ylabel(r'$\nu\ {\rm GHz}$')
    p.title(r'$|\vec u|=%d$' % (umag))
fig.canvas.set_window_title(opts.name)
p.show()
