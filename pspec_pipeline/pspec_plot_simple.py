#! /usr/bin/env python
import matplotlib
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse, re, os

o=optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])
ONLY_POS_K = True

def dual_plot(kpl, pk, err, pkfold=None, errfold=None, umag=16., f0=.159, color='', bins=None, upperlimit=False):
    z = C.pspec.f2z(f0)
    kpr = C.pspec.dk_du(z) * umag
    k = n.sqrt(kpl**2 + kpr**2)
    k3 = n.abs(k**3 / (2*n.pi**2))
    print 'k [h Mpc^-1], P(k) [K^2], err (2sigma)'
    for _k,_pk,_err in zip(kpl,pk,err):
        print '%6.3f, %9.5f, %9.5f' % (_k, _pk.real/1e6, _err/1e6)
    print '-'*20
    pk = pk.real
    p.subplot(121)
    if not upperlimit: p.errorbar(kpl, pk, yerr=err, fmt=color+'.', capsize=0, linewidth=1.2)
    else: p.plot(kpl, pk+err, color+'-')
    p.subplot(122)
    k0 = n.abs(kpl).argmin()
    if pkfold is None:
        print 'Folding'
        pkfold = pk[k0:].copy()
        errfold = err[k0:].copy()
        pkpos,errpos = pk[k0+1:].copy(), err[k0+1:].copy()
        pkneg,errneg = pk[k0-1::-1].copy(), err[k0-1::-1].copy()
        pkfold[1:] = (pkpos/errpos**2 + pkneg/errneg**2) / (1./errpos**2 + 1./errneg**2)
        errfold[1:] = n.sqrt(1./(1./errpos**2 + 1./errneg**2))
    else:
        pkfold = pkfold.real

    if not upperlimit: p.errorbar(k[k0:], k3[k0:]*pkfold, yerr=k3[k0:]*errfold, fmt=color+'.', capsize=0,linewidth=1.2)
    else: p.plot(k[k0:], k3[k0:]*pkfold+k3[k0:]*errfold, color+'-')
    if not bins is None:
        _kpls, _k3pks, _k3errs = [], [], []
        for (dn,up) in bins:
            KRMS = False
            if dn < 0 and up > 0: KRMS = True
            ksum,kwgt = 0,0
            dsum,dwgt = 0., 0.
            for _kpl,_pk,_err in zip(kpl, pk, err):
                if dn < _kpl and _kpl <= up:
                    dsum += _pk / _err**2
                    dwgt += 1. / _err**2
                    if KRMS: ksum += _kpl**2 # krms
                    else: ksum += _kpl # kavg
                    kwgt += 1.
            if KRMS: kavg = n.sqrt(ksum/kwgt + kpr**2) # krms
            else: kavg = n.sqrt((ksum/kwgt)**2 + kpr**2) # kavg
            if dwgt == 0: continue
            _pk = dsum / dwgt
            _err = 1. / n.sqrt(dwgt)
            _k3pk = kavg**3/(2*n.pi**2) * _pk
            _k3err = kavg**3/(2*n.pi**2) * _err
            _kpls.append(dn); _kpls.append(.5*(dn+up)); _kpls.append(up)
            _k3pks.append(_k3pk); _k3pks.append(_k3pk); _k3pks.append(_k3pk)
            _k3errs.append(_k3err); _k3errs.append(_k3err); _k3errs.append(_k3err)
        kpl = n.array(_kpls)
        k3pk = n.abs(n.array(_k3pks))
        k3err = n.array(_k3errs)
    else:
        kpl = kpl
        k3pk = n.abs(k3*pk) 
        k3err = k3*err
    for _k,_k3pk,_k3err in zip(kpl,pk,err):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)
    print '-'*20
    for _k,_k3pk,_k3err in zip(k[k0:],k3[k0:]*pkfold,k3[k0:]*errfold):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)
    print '-'*20
    print "saving pspec_pk_k3pk.npz"
    print "output @ freq = ",f0
    n.savez('pspec_pk_k3pk.npz',kpl=kpl,pk=pk,err=err,k=k[k0:], k3pk=k3[k0:]*pkfold, k3err=k3[k0:]*errfold,freq=f0)

RS_VS_KPL = {} # K^2
RS_VS_KPL_FOLD = {} # K^2
dsum, dwgt = {}, {}
dsum_fold, dwgt_fold = {}, {}
afreqs=[]
chans=[]
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    afreqs=f['afreqs']
    chans=f['chans']
    RS_VS_KPL[filename] = {}
    RS_VS_KPL_FOLD[filename] = {}
    kpl,pk,err = f['kpl'], f['pk'], f['err']
    try: pkfold,errfold = f['pk_fold'],f['err_fold']
    except(KeyError): pkfold,errfold = None, None
    for _kpl, _pk, _err in zip(kpl, pk, err):
        RS_VS_KPL[filename][_kpl] = (_pk, _err)
        dsum[_kpl] = dsum.get(_kpl, 0) + _pk / _err**2
        dwgt[_kpl] = dwgt.get(_kpl, 0) + 1 / _err**2
    k0 = n.abs(kpl).argmin()
    if not pkfold is None:
        for _kpl, _pk, _err in zip(kpl[k0:], pkfold, errfold):
            RS_VS_KPL_FOLD[filename][_kpl] = (_pk, _err)
            dsum_fold[_kpl] = dsum_fold.get(_kpl, 0) + _pk / _err**2
            dwgt_fold[_kpl] = dwgt_fold.get(_kpl, 0) + 1 / _err**2
freq=n.average(afreqs)
if True:
    RS_VS_KPL['total'] = {}
    RS_VS_KPL_FOLD['total'] = {}
    for _kpl in dsum:
        RS_VS_KPL['total'][_kpl] = (dsum[_kpl] / dwgt[_kpl], 1./n.sqrt(dwgt[_kpl]))
    for _kpl in dsum_fold:
        RS_VS_KPL_FOLD['total'][_kpl] = (dsum_fold[_kpl] / dwgt_fold[_kpl], 1./n.sqrt(dwgt_fold[_kpl]))

fig=p.figure(figsize=(12,7.2))

BINS = None
colors = 'kbcm' * 10
for sep in RS_VS_KPL:
    if not sep == 'total': continue
    dsum, dwgt = {}, {}
    dsum_fold, dwgt_fold = {}, {}
    ks = RS_VS_KPL[sep].keys(); ks.sort()
    for k in ks:
        _d,_n = RS_VS_KPL[sep][k]
        dsum[k] = dsum.get(k,0) + _d / _n**2
        dwgt[k] = dwgt.get(k,0) + 1 / _n**2
    ks_fold = RS_VS_KPL_FOLD[sep].keys(); ks_fold.sort()
    for k in ks_fold:
        _d,_n = RS_VS_KPL_FOLD[sep][k]
        dsum_fold[k] = dsum_fold.get(k,0) + _d / _n**2
        dwgt_fold[k] = dwgt_fold.get(k,0) + 1 / _n**2
        
    kpl = dsum.keys(); kpl.sort()
    d = [dsum[k]/dwgt[k] for k in kpl]
    nos = [1./n.sqrt(dwgt[k]) for k in kpl]
    kpl_fold = dsum_fold.keys(); kpl_fold.sort()
    d_fold = [dsum_fold[k]/dwgt_fold[k] for k in kpl_fold]
    nos_fold = [1./n.sqrt(dwgt_fold[k]) for k in kpl_fold]
    d,kpl,nos = n.array(d, dtype=n.complex), n.array(kpl), n.array(nos)
    d_fold,kpl_fold,nos_fold = n.array(d_fold, dtype=n.complex), n.array(kpl_fold), n.array(nos_fold)
    if True: #MEDIAN STATISTIC CORRECTION
        f = 1/n.log(2)
        print 'Scaling data and noise by %f for using the median statistic.'%f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    
    if d_fold.size == 0: d_fold,nos_fold = None, None
    d_fold = n.abs(d_fold) #XXX ABSOLUTE VALUE for DELTA SQ!
    dual_plot(kpl, d, 2*nos, d_fold, 2*nos_fold, color=colors[0], bins=BINS,f0=freq) # 2-sigma error bars
    colors = colors[1:] + colors[0]

print 'Average of P(k) = ',n.median(d)

tau_h = 100 + 15. #in ns
k_h = C.pspec.dk_deta(C.pspec.f2z(.151)) * tau_h
p.subplot(121)
#p.vlines(k_h, -1e7, 1e8, linestyles='--', linewidth=1.5)
#p.vlines(-k_h, -1e7, 1e8, linestyles='--', linewidth=1.5)
p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
p.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
p.grid()

p.subplot(122)
#p.vlines(k_h, -1e7, 1e7, linestyles='--', linewidth=1.5)
p.gca().set_yscale('log', nonposy='clip')
p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$', fontsize='large')
p.xlim(0, 0.6)
p.grid()
#plot noise_curve that's saved from plot_pk_k3pk_zsa_2.py
try:
    npz = n.load('noise_curve.npz')
    xs = npz['kpls']
    ys = npz['noise']
    p.plot(xs,ys,'c--')
except: pass
p.show()

