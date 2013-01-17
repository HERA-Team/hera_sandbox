#! /usr/bin/env python
import numpy as n, pylab as p
import sys, scipy.optimize

def parse_line(L):
    # N | fq | flx | err | units | fq (Hz) | flx (Jy) | err (Jy) | Jy | ref | # sigma | band | ? | crd | ? | ?
    w = L.split('|')
    print w[5:8] + w[10:11]
    fq = float(w[5])
    try: jy = float(w[6])
    except(ValueError): jy = 0
    err = w[7]
    if err.startswith('+/-'): err = err[len('+/-'):]
    try: err = float(err)
    except(ValueError): err = -1
    nsig = w[10]
    if nsig.find('sigma') != -1: nsig = float(nsig.split()[0])
    else: nsig = 1
    sig = err / nsig
    print fq/1e9, jy, sig
    return fq/1e9, jy, sig

def pwrlaw(fqs, flx, ind, mfreq=.150):
    return flx * (fqs/mfreq)**ind
    
def fit_pwrlaw(fqs, dat, err, flx, ind, mfreq=.150):
    spec = pwrlaw(fqs, flx, ind, mfreq=mfreq)
    scr = n.abs(dat-spec)**2 / err**4
    wgt = 1. / err**4
    rv = n.sqrt(n.sum(scr)/n.sum(wgt))
    return rv

for f in sys.argv[1:]:
    print 'Reading', f
    if f.endswith('txt'):
        d = n.array([parse_line(L) for L in open(f).readlines() if len(L.split('|')) == 17 and not L.startswith('#')])
        fq,jy,err = d[:,0], d[:,1], d[:,2]  
        err = n.where(err <= 0, jy, err) # approx unreported errors as measured value
        #usedat = n.where(fq < .7, 1, 0)
        #fq,jy,err = fq.compress(usedat), jy.compress(usedat), err.compress(usedat)
        def fitfunc(prms): return fit_pwrlaw(fq, jy, err, prms[0], prms[1])
        rv = scipy.optimize.fmin(fitfunc, n.array([0.,0]), full_output=1, disp=0, ftol=1e-4, xtol=1e-4)
        prms,score = rv[:2]
        flx,ind = prms
        print flx, ind, score
        p.plot(fq, pwrlaw(fq, flx, ind))
        p.errorbar(fq, jy, yerr=err, fmt='.', label=f)
    else:
        npz = n.load(f)
        fq = npz['freq']
        jy = npz['spec']
        def fitfunc(prms): return fit_pwrlaw(fq, jy, 1, prms[0], prms[1])
        rv = scipy.optimize.fmin(fitfunc, n.array([0.,0]), full_output=1, disp=0, ftol=1e-4, xtol=1e-4)
        prms,score = rv[:2]
        flx,ind = prms
        print flx, ind, score
        p.plot(fq, jy, '.', label=f)
        p.plot(fq, pwrlaw(fq, flx, ind))

p.legend(loc='best')
p.gca().set_xscale('log', nonposy='clip')
p.gca().set_yscale('log', nonposy='clip')
p.xlim(.05, 2)
p.ylim(.1,1e3)
p.show()
