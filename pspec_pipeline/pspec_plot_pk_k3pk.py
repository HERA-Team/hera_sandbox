#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse, re, os

o=optparse.OptionParser()
o.add_option('--flux', action='store_true', 
                help='Scale data due to flux calibration errors from pictor. Scales by factor in f option.')
o.add_option('-f', dest='flux_factor', action='store', type='float', default=.736,
            help='scaling factor for flux')
o.add_option('--beam', action='store_true',
                help='scale data by beam square instead of beam')
o.add_option('--afrf', action='store_true',
                help='scale data by factor from aggresive fringe rate filtering.')
o.add_option('-a', dest='afrf_factor', action='store', type='float', default=1.9,
            help='scaling factor for aggressive fring rate filtering.')
o.add_option('--cov', action='store_true',
            help='scale factor for signal loss in covariance removal')
o.add_option('--show', action='store_true',
            help='Show the plot')
opts,args = o.parse_args(sys.argv[1:])
print args


ONLY_POS_K = True

def dual_plot(kpl, pk, err, pkfold=None, errfold=None, umag=16., f0=.164, color='', bins=None):
    z = C.pspec.f2z(f0)
    kpr = C.pspec.dk_du(z) * umag
    k = n.sqrt(kpl**2 + kpr**2)
    k3 = n.abs(k**3 / (2*n.pi**2))
    print 'k [h Mpc^-1], P(k) [K^2], err (2sigma)'
    for _k,_pk,_err in zip(kpl,pk,err):
        print '%6.3f, %9.5f, %9.5f' % (_k, _pk.real/1e6, _err/1e6)
        #print '%6.3f, %9.5f, %9.5f' % (_k, _pk.imag/1e6, _err/1e6)
    print '-'*20
    #pk = pk.imag
    pk = pk.real
    p.subplot(121)
    #pk = n.abs(pk)
    #p.errorbar(kpl, pk, yerr=err, fmt=color+'.-')
    #p.errorbar(kpl, n.abs(pk), yerr=err, fmt='mx')
    p.errorbar(kpl, pk, yerr=err, fmt=color+'.', capsize=0)
    #p.errorbar(kpl, pk.imag, yerr=err, fmt=color+'x')
    p.subplot(122)
    k0 = n.abs(kpl).argmin()
    if pkfold is None:
        print 'Folding'
        pkfold = pk[k0:].copy()
        errfold = err[k0:].copy()
        pkpos,errpos = pk[k0+1:].copy(), err[k0+1:].copy()
        pkneg,errneg = pk[k0-1:0:-1].copy(), err[k0-1:0:-1].copy()
        pkfold[1:] = (pkpos/errpos**2 + pkneg/errneg**2) / (1./errpos**2 + 1./errneg**2)
        errfold[1:] = n.sqrt(1./(1./errpos**2 + 1./errneg**2))
    else:
        print pkfold.imag
        pkfold = pkfold.real
        #pkfold = pkfold.imag

    #p.errorbar(k, k3*pk, yerr=k3*err, fmt=color+'.', capsize=0)
    p.errorbar(k[k0:], k3[k0:]*pkfold, yerr=k3[k0:]*errfold, fmt=color+'.', capsize=0)
    if not bins is None:
        _kpls, _k3pks, _k3errs = [], [], []
        for (dn,up) in bins:
            KRMS = False
            if dn < 0 and up > 0: KRMS = True
            ksum,kwgt = 0,0
            dsum,dwgt = 0., 0.
            for _kpl,_pk,_err in zip(kpl, pk, err):
                if dn < _kpl and _kpl <= up:
                    #print _pk, _err
                    dsum += _pk / _err**2
                    dwgt += 1. / _err**2
                    #if KRMS: ksum += _kpl**2 / _err**2 # krms
                    #else: ksum += _kpl / _err**2 # kavg
                    #kwgt += 1. / _err**2
                    if KRMS: ksum += _kpl**2 # krms
                    else: ksum += _kpl # kavg
                    kwgt += 1.
            if KRMS: kavg = n.sqrt(ksum/kwgt + kpr**2) # krms
            else: kavg = n.sqrt((ksum/kwgt)**2 + kpr**2) # kavg
            if dwgt == 0: continue
            _pk = dsum / dwgt
            #print _pk
            _err = 1. / n.sqrt(dwgt)
            _k3pk = kavg**3/(2*n.pi**2) * _pk
            _k3err = kavg**3/(2*n.pi**2) * _err
            _kpls.append(dn); _kpls.append(.5*(dn+up)); _kpls.append(up)
            _k3pks.append(_k3pk); _k3pks.append(_k3pk); _k3pks.append(_k3pk)
            _k3errs.append(_k3err); _k3errs.append(_k3err); _k3errs.append(_k3err)
        kpl = n.array(_kpls)
        k3pk = n.array(_k3pks)
        k3err = n.array(_k3errs)
    else:
        kpl = kpl
        k3pk = k3*pk
        k3err = k3*err
    #p.plot(kpl, k3*pk+k3*err, color+'.-')
    #p.plot(kpl, k3pk+k3err, color+'.-')
    for _k,_k3pk,_k3err in zip(kpl,pk,err):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)
    print '-'*20
    for _k,_k3pk,_k3err in zip(k[k0:],k3[k0:]*pkfold,k3[k0:]*errfold):
        print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)
    print '-'*20
    print "saving pspec_pk_k3pk.npz"
    print "output @ freq = ",f0
    n.savez('pspec_pk_k3pk.npz',kpl=kpl,pk=pk,err=err,k=k[k0:], k3pk=k3[k0:]*pkfold, k3err=k3[k0:]*errfold,freq=f0)
    #pos = n.where(kpl >= 0, 1, 0)
    #neg = n.where(kpl <= 0, 1, 0)
    #posneg = 0.5*(k3pk.compress(pos) + k3pk.compress(neg)[::-1])
    #posneg_err = n.sqrt(1./(1./k3err.compress(pos)**2 + 1./k3err.compress(neg)[::-1]**2))
    #for _k,_k3pk,_k3err in zip(kpl.compress(pos),posneg,posneg_err):
    #    print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)
    #if ONLY_POS_K:
    #    #p.plot(n.sqrt(kpr**2+kpl.compress(pos)**2), k3pk.compress(pos)+k3err.compress(pos), 'c-')
    #    ##p.plot(n.sqrt(kpr**2+kpl.compress(pos)**2), k3pk.compress(pos)-k3err.compress(pos), 'c-')
    #    #p.plot(n.sqrt(kpr**2+kpl.compress(neg)**2), k3pk.compress(neg)+k3err.compress(neg), 'm-')
    #    p.plot(n.sqrt(kpr**2+kpl.compress(pos)**2), posneg+posneg_err, 'k-')
    #    #p.plot(n.sqrt(kpr**2+kpl.compress(neg)**2), k3pk.compress(neg)-k3err.compress(neg), 'm-')
    #    #p.plot(n.sqrt(kpr**2 + kpl**2), k3pk+k3err, color+'-')
    #else: p.plot(kpl, k3pk+k3err, color+'-')

#o = optparse.OptionParser()
#opts,args = o.parse_args(sys.argv[1:])
#args = ['data/pspec_t1_c110-149.npz']
#args = sys.argv[1:]

FG_VS_KPL_NOS = 168.74e6
FG_VS_KPL = { # K^2
#    '-0.054':   5.37262e+13,
#    '-0.027':   7.15304e+14, 
    ' 0.000':   3.50958e+15, 
#    ' 0.027':   4.12396e+14, 
#    ' 0.054':   2.60795e+13,
    #-0.0536455587089:   5.37262e+13,
    #-0.0268227793545:   7.15304e+14, 
    #0.0:                3.50958e+15, 
    #0.0268227793545:    4.12396e+14, 
    #0.0536455587089:    2.60795e+13,
}

RS_VS_KPL = {} # K^2
RS_VS_KPL_FOLD = {} # K^2
dsum, dwgt = {}, {}
dsum_fold, dwgt_fold = {}, {}
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    RS_VS_KPL[filename] = {}
    RS_VS_KPL_FOLD[filename] = {}
    kpl,pk,err = f['kpl'], f['pk'], f['err']
    try: pkfold,errfold = f['pk_fold'],f['err_fold']
    except(KeyError): pkfold,errfold = None, None
    if False: # Hacky way to get a noise bias out, if necessary
        pk -= n.median(n.concatenate([pk[:8], pk[-8:]]))
    if False: # Hacky way to estimate noise
        print 'Overriding errors for %s:' % filename
        print 'Old err:'
        print err
        err = n.std(n.concatenate([pk[:8], pk[-8:]])) * n.ones_like(pk).real
        print 'New err:'
        print err
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
#freq = f['freq']
#RS_VS_KPL = {}
if True:
    RS_VS_KPL['total'] = {}
    RS_VS_KPL_FOLD['total'] = {}
    for _kpl in dsum:
        RS_VS_KPL['total'][_kpl] = (dsum[_kpl] / dwgt[_kpl], 1./n.sqrt(dwgt[_kpl]))
    for _kpl in dsum_fold:
        RS_VS_KPL_FOLD['total'][_kpl] = (dsum_fold[_kpl] / dwgt_fold[_kpl], 1./n.sqrt(dwgt_fold[_kpl]))

BINS = None
#BINS = ((-.52,-.4),(-.4,-.2),(-.2,-.15),(-.15,-.1),(-.1,-.06),(-.06,.06),(.06,.1),(.1,.15),(.15,.2),(.2,.4),(.4,.52))
#BINS = ((-.5,-.25),(-.25,-.125),(-.125,-.0625),(-.0625,.0625),(.0625,.125),(.125,.25),(.25,.5))
#BINS = ((-.5,-.375),(-.375,-.25),(-.25,-.125),(-.125,-.0625),(-.0625,.0625),(.0625,.125),(.125,.25),(.25,.375),(.375,.5))
#BINS = ((-.5,-.375),(-.375,-.25),(-.25,-.175),(-.175,-.125),(-.125,-.1),(-.1,-.07),(-.07,-.04),(-.04,.04),(.04,.07),(.07,.1),(.1,.125),(.125,.175),(.175,.25),(.25,.375),(.375,.5))
#BINS = ((-.6,-.33),(-.33,-.22),(-.22,-.1),(-.1,.1),(.1,.22),(.22,.33),(.33,.6))
#BINS = ((-.6,-.5),(-.5,-.4),(-.4,-.3),(-.3,-.2),(-.2,-.1),(-.1,-.05),(-.05,.05),(.05,.1),(.1,.2),(.2,.3),(.3,.4),(.4,.5),(.5,.6))
#kpl = RS_VS_KPL.keys(); kpl.sort()
#d = [RS_VS_KPL[k] for k in kpl]
colors = 'kbcm' * 10
for sep in RS_VS_KPL:
    if not sep == 'total': continue
    dsum, dwgt = {}, {}
    dsum_fold, dwgt_fold = {}, {}
    ks = RS_VS_KPL[sep].keys(); ks.sort()
    for k in ks:
        _d,_n = RS_VS_KPL[sep][k]
        #print k, _d, _n, '->',
        dsum[k] = dsum.get(k,0) + _d / _n**2
        dwgt[k] = dwgt.get(k,0) + 1 / _n**2
        #print dsum[k] / dwgt[k], 1./n.sqrt(dwgt[k])
    ks_fold = RS_VS_KPL_FOLD[sep].keys(); ks_fold.sort()
    for k in ks_fold:
        _d,_n = RS_VS_KPL_FOLD[sep][k]
        #print k, _d, _n, '->',
        dsum_fold[k] = dsum_fold.get(k,0) + _d / _n**2
        dwgt_fold[k] = dwgt_fold.get(k,0) + 1 / _n**2
        #print dsum[k] / dwgt[k], 1./n.sqrt(dwgt[k])
        
    kpl = dsum.keys(); kpl.sort()
    d = [dsum[k]/dwgt[k] for k in kpl]
    nos = [1./n.sqrt(dwgt[k]) for k in kpl]
    kpl_fold = dsum_fold.keys(); kpl_fold.sort()
    d_fold = [dsum_fold[k]/dwgt_fold[k] for k in kpl_fold]
    nos_fold = [1./n.sqrt(dwgt_fold[k]) for k in kpl_fold]
    #d = [RS_VS_KPL[k][0] for k in kpl]
    #nos = [RS_VS_KPL[k][1] for k in kpl]
    if True: #if 'I' in sep: # Add foregrounds
        for cnt,k in enumerate(kpl):
            k = '%6.3f' % k
            if not FG_VS_KPL.has_key(k): continue
            d[cnt] += FG_VS_KPL[k]
            nos[cnt] = 2*n.sqrt(FG_VS_KPL_NOS*FG_VS_KPL[k])
        for cnt,k in enumerate(kpl_fold):
            k = '%6.3f' % k
            if not FG_VS_KPL.has_key(k): continue
            d_fold[cnt] += FG_VS_KPL[k]
            nos_fold[cnt] = 2*n.sqrt(FG_VS_KPL_NOS*FG_VS_KPL[k])
        
    d,kpl,nos = n.array(d, dtype=n.complex), n.array(kpl), n.array(nos)
    d_fold,kpl_fold,nos_fold = n.array(d_fold, dtype=n.complex), n.array(kpl_fold), n.array(nos_fold)
    if opts.flux:
        # PSA32 was calibrated to Pictor A @ 160 MHz = 424 Jy
        # To recalibrate to new Pic A, must multiply by square of ratio of fluxes
        # Jacobs et al 2013 says Pic A = 382 @ 150 MHz, index=-0.76, so at 160 MHz, Pic A = 364 Jy
        #f = 0.76 # psa747 calibration of Pic A = 370.6 Jy @ 160 MHz (which includes resolution effects)
        #f = 0.736 # rescale by (364/424)**2 to correct flux scale
        f = opts.flux_factor
        print 'Scaling data and noise by %f for recalibration to PicA from Jacobs et al. 2013 (PSA32 only)' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    if opts.beam:
        f = 2.35 # Use power**2 beam, which is a 1.69/0.72=2.35 penalty factor
        print 'Scaling data and noise by %f for correcting cosmo scalar to use power^2 beam' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    if opts.afrf: # For aggressive fringe-rate filtering, change beam area
        f = opts.afrf_factor
        f = 1.90 # ratio of power**2 beams for filtered * unfiltered beams: 0.306 / 0.162
        print 'Scaling data and noise by %f for beam constriction in aggressive fringe-rate filtering' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    if opts.cov: # extra penalty for signal loss in covariance diagonalization
        f = 1.2
        print 'Scaling data and noise by %f for signal loss in covariance diagonalization' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    #for _kpl,_pk,_nos in zip(kpl,d,nos): print _kpl, _pk, _nos
    print sep, colors[0]
    '''
    if True: # Hacky way to get a noise bias out, if necessary
        d -= n.median(n.concatenate([d[:8], d[-8:]]))
    if True: # Hacky way to estimate noise
        nos = n.std(n.concatenate([d[:8], d[-8:]])) * n.ones_like(d)
    '''
    if d_fold.size == 0: d_fold,nos_fold = None, None
    dual_plot(kpl, d, 2*nos, d_fold, 2*nos_fold, color=colors[0], bins=BINS)#,f0=freq) # 2-sigma error bars
    #dual_plot(kpl, d, nos, color=colors[0], bins=BINS) # 2-sigma error bars
    colors = colors[1:] + colors[0]

def mean_temp(z):
    return 28. * ((1.+z)/10.)**.5 # mK

import glob
re_z = re.compile(r'power_21cm_z(\d+\.\d+)\.dat')

for filename in glob.glob('lidz_mcquinn_k3pk/*7.3*dat'):
    print 'Reading', filename
    d = n.array([map(float, L.split()) for L in open(filename).readlines()])
    ks, pk = d[:,0], d[:,1]
    z_file = float(re_z.match(os.path.basename(filename)).groups()[0])
    z = C.pspec.f2z(.160)
    k3pk = ks**3 / (2*n.pi**2) * pk
    p.subplot(122)
    p.plot(ks, k3pk * mean_temp(z)**2, 'k-')
    
p.subplot(121)
p.gca().set_yscale('log', nonposy='clip')
p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
p.ylabel(r'$P(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$')
p.ylim(1e5,3e16)
p.grid()


p.subplot(122)
if ONLY_POS_K: p.plot([.5], [248**2], 'mv', label='GMRT2013')
else: p.plot([-.5, .5], [248**2, 248**2], 'mv', label='GMRT2013')


p.gca().set_yscale('log', nonposy='clip')
p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$')
p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$')
p.ylim(1e0,1e9)
p.xlim(0, 0.6)
p.grid()
p.savefig('pspec.png')

#p.figure(2)
#dual_plot(kpl, d, 2*nos, d_fold, 2*nos_fold, color=colors[0], bins=BINS,f0=freq) # 2-sigma error bars
#p.subplot(121)
#p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
#p.ylabel(r'$P(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$')
#p.ylim(-1e2,1e5)
#p.grid()
#p.subplot(122)
#p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$')
#p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$')
#p.ylim(-1e3,1e3)
#p.xlim(0, 0.6)
#p.grid()

f = n.load(args[0])
def posterior(kpl, pk, err, pkfold=None, errfold=None):
    k0 = n.abs(kpl).argmin()
    kpl = kpl[k0:]
    if pkfold is None:
        print 'Folding for posterior'
        pkfold = pk[k0:].copy()
        errfold = err[k0:].copy()
        pkpos,errpos = pk[k0+1:].copy(), err[k0+1:].copy()
        pkneg,errneg = pk[k0-1:0:-1].copy(), err[k0-1:0:-1].copy()
        pkfold[1:] = (pkpos/errpos**2 + pkneg/errneg**2) / (1./errpos**2 + 1./errneg**2)
        errfold[1:] = n.sqrt(1./(1./errpos**2 + 1./errneg**2))

    #ind = n.logical_and(kpl>.2, kpl<.5)
    ind = n.logical_and(kpl>.15, kpl<.5)
    #ind = n.logical_and(kpl>.12, kpl<.5)
    #print kpl,pk.real,err
    kpl = kpl[ind]
    pk= kpl**3 * pkfold[ind]/(2*n.pi**2)
    err = kpl**3 * errfold[ind]/(2*n.pi**2)
    s = n.logspace(1,3.5,100)
    data = []
    for ss in s:
        data.append(n.exp(-.5*n.sum((pk.real - ss)**2 / err**2)))
    #    print data[-1]
    data = n.array(data)
    #print data
    #print s
    #data/=n.sum(data)
    data /= n.max(data)
    p.figure(5)
    p.plot(s, data)
    p.plot(s, n.exp(-.5)*n.ones_like(s))
    p.plot(s, n.exp(-.5*2**2)*n.ones_like(s))
    p.show()

#posterior(f['kpl'], f['pk'], f['err'], f['pk_fold'], f['err_fold'])
posterior(kpl, d, nos, d_fold, nos_fold)


if opts.show:
    p.show()
