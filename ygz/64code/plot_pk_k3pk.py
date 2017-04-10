#! /usr/bin/env python
import aipy as a, numpy as np, pylab as p
import capo as C
import sys, optparse, re, os

o=optparse.OptionParser()
o.add_option('--flux', action='store_true', 
                help='Scale data due to flux calibration errors from pictor. Scales by factor in f optionp.')
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
o.add_option('--med', action='store_true',
            help='correction factor for using median statistics instead of the mean')
o.add_option('--show', action='store_true',
            help='Show the plot')
opts,args = o.parse_args(sys.argv[1:])
print args

#dspec='data/final_pspecs/v032/nocov_pspec.npz'
#args=['data/final_pspecs/v050/pspec.npz']
#args=['data/final_pspecs/v0200/pspec.npz']
#args=['data/final_pspecs/v060/pspec.npz']

def noise_level():
    tsys = 500e3 #mK
    inttime = 1886. #seconds. XXX fix with new integrationp.
    #inttime = 2930. #seconds. XXX fix with new integrationp.
    nbls=51
    ndays = 120 #effectively this many days. days of operation = 135
    nseps = 3
    folding = 2
    #ndays = 60 #effectively this many days. days of operation = 135
    nmodes = (nseps*folding*8.5*60*60/inttime)**.5 # 3 for baseline types, 8.5=total lst time used.
    pol = 2
    real = 2 #??? what is this again?
    freq = .15117 #GHz
    z = C.pspec.f2z(freq)
    X2Y = C.pspec.X2Y(z)/1e9 #h^-3 Mpc^3 / str/ Hz
    sdf = .1/203
    freqs = np.linspace(.1,.2,203)[95:115]
    B = sdf*freqs.size
    bm = np.polyval(C.pspec.DEFAULT_BEAM_POLY, freq) * 2.35 # correction for beam^2
    scalar = X2Y * bm 

    fr_correct = 1.39
    #fr_correct = 1.90

    print 'X2Y', X2Y
    print 'scalar:', scalar
    print 'BM:', bm
    print 'Tsys:', tsys
    print 

    #error bars minimum width. Consider them flat for P(k). Factor of 2 at the end is due to folding of kpl (root(2)) and root(2) in radiometer equationp.
    pk = scalar*fr_correct*( (tsys)**2 / (2*inttime*pol*real*nbls*ndays*nmodes) )

    print 'Pk noise (1 sigma):', pk
    #
#    etas = np.fft.fftshift(C.pspec.f2eta(freqs))
#    kpl = etas * C.pspec.dk_deta(z)
#    kpl_pos = kpl[kpl>=0]
    ks = np.linspace(.05837,.6,100)

    d2 = pk*ks**3 / (2*np.pi**2)
    #return ks, d2
    return pk



ONLY_POS_K = True

def dual_plot(kpl, pk, err, pkfold=None, errfold=None, umag=16., f0=.164, color='', bins=None, upperlimit=False):
    z = C.pspec.f2z(f0)
    kpr = C.pspec.dk_du(z) * umag
    k = np.sqrt(kpl**2 + kpr**2)
    k3 = np.abs(k**3 / (2*np.pi**2))
    print 'k [h Mpc^-1], P(k) [K^2], err (2sigma)'
    for _k,_pk,_err in zip(kpl,pk,err):
        print '%6.3f, %9.5f, %9.5f' % (_k, _pk.real/1e6, _err/1e6)
        #print '%6.3f, %9.5f, %9.5f' % (_k, _pk.real/1e6, _err/1e6)
    print '-'*20
    #pk = pk.real
    pk = np.abs(pk)
    p.subplot(121)
    #pk = np.abs(pk)
    #p.errorbar(kpl, pk, yerr=err, fmt=color+'.-')
    #p.errorbar(kpl, np.abs(pk), yerr=err, fmt='mx')
    if not upperlimit: p.errorbar(kpl, pk, yerr=err, fmt=color+'.', capsize=0, linewidth=1.2)
    else: p.plot(kpl, pk+err, color+'-')
    #p.errorbar(kpl, pk.real, yerr=err, fmt=color+'x')
    p.subplot(122)
    k0 = np.abs(kpl).argmin()
    if pkfold is None:
        print 'Folding'
        pkfold = pk[k0:].copy()
        errfold = err[k0:].copy()
        pkpos,errpos = pk[k0+1:].copy(), err[k0+1:].copy()
        pkneg,errneg = pk[k0-1::-1].copy(), err[k0-1::-1].copy()
#        print pkneg, pkpos
        pkfold[1:] = (pkpos/errpos**2 + pkneg/errneg**2) / (1./errpos**2 + 1./errneg**2)
        errfold[1:] = np.sqrt(1./(1./errpos**2 + 1./errneg**2))
    else:
        #print pkfold.real
        pkfold = np.abs(pkfold)
        #pkfold = pkfold.real

    #p.errorbar(k, k3*pk, yerr=k3*err, fmt=color+'.', capsize=0)
    #print color
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
                    #print _pk, _err
                    dsum += _pk / _err**2
                    dwgt += 1. / _err**2
                    #if KRMS: ksum += _kpl**2 / _err**2 # krms
                    #else: ksum += _kpl / _err**2 # kavg
                    #kwgt += 1. / _err**2
                    if KRMS: ksum += _kpl**2 # krms
                    else: ksum += _kpl # kavg
                    kwgt += 1.
            if KRMS: kavg = np.sqrt(ksum/kwgt + kpr**2) # krms
            else: kavg = np.sqrt((ksum/kwgt)**2 + kpr**2) # kavg
            if dwgt == 0: continue
            _pk = dsum / dwgt
            #print _pk
            _err = 1. / np.sqrt(dwgt)
            _k3pk = kavg**3/(2*np.pi**2) * _pk
            _k3err = kavg**3/(2*np.pi**2) * _err
            _kpls.append(dn); _kpls.append(.5*(dn+up)); _kpls.append(up)
            _k3pks.append(_k3pk); _k3pks.append(_k3pk); _k3pks.append(_k3pk)
            _k3errs.append(_k3err); _k3errs.append(_k3err); _k3errs.append(_k3err)
        kpl = np.array(_kpls)
        k3pk = np.array(_k3pks)
        k3err = np.array(_k3errs)
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
    np.savez('pspec_pk_k3pk.npz',kpl=kpl,pk=pk,err=err,k=k[k0:], k3pk=k3[k0:]*pkfold, k3err=k3[k0:]*errfold,freq=f0)
    #pos = np.where(kpl >= 0, 1, 0)
    #neg = np.where(kpl <= 0, 1, 0)
    #posneg = 0.5*(k3pk.compress(pos) + k3pk.compress(neg)[::-1])
    #posneg_err = np.sqrt(1./(1./k3err.compress(pos)**2 + 1./k3err.compress(neg)[::-1]**2))
    #for _k,_k3pk,_k3err in zip(kpl.compress(pos),posneg,posneg_err):
    #    print '%6.3f, %9.5f (%9.5f +/- %9.5f)' % (_k, _k3pk+_k3err,_k3pk,_k3err)
    #if ONLY_POS_K:
    #    #p.plot(np.sqrt(kpr**2+kpl.compress(pos)**2), k3pk.compress(pos)+k3err.compress(pos), 'c-')
    #    ##p.plot(np.sqrt(kpr**2+kpl.compress(pos)**2), k3pk.compress(pos)-k3err.compress(pos), 'c-')
    #    #p.plot(np.sqrt(kpr**2+kpl.compress(neg)**2), k3pk.compress(neg)+k3err.compress(neg), 'm-')
    #    p.plot(np.sqrt(kpr**2+kpl.compress(pos)**2), posneg+posneg_err, 'k-')
    #    #p.plot(np.sqrt(kpr**2+kpl.compress(neg)**2), k3pk.compress(neg)-k3err.compress(neg), 'm-')
    #    #p.plot(np.sqrt(kpr**2 + kpl**2), k3pk+k3err, color+'-')
    #else: p.plot(kpl, k3pk+k3err, color+'-')

#o = optparse.OptionParser()
#opts,args = o.parse_args(sys.argv[1:])
#args = ['data/pspec_t1_c110-149.npz']
#args = sys.argv[1:]

FG_VS_KPL_NOS = 168.74e6
FG_VS_KPL = { # K^2
#    '-0.054':   5.37262e+13,
    '-0.049':   1.06015e+15,
    ' 0.000':   8.44e+15, 
    ' 0.049':   1.22195e+15, 

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
    f = np.load(filename)
    RS_VS_KPL[filename] = {}
    RS_VS_KPL_FOLD[filename] = {}
    kpl,pk,err = f['kpl'], f['pk'], f['err']
    try: pkfold,errfold = f['pk_fold'],f['err_fold']
    except(KeyError): pkfold,errfold = None, None
    if False: # Hacky way to get a noise bias out, if necessary
        pk -= np.median(np.concatenate([pk[:8], pk[-8:]]))
    if False: # Hacky way to estimate noise
        print 'Overriding errors for %s:' % filename
        print 'Old err:'
        print err
        err = np.abs(np.std(np.concatenate([pk[:8], pk[-8:]])) * np.ones_like(pk))
        print 'New err:'
        print err
    for _kpl, _pk, _err in zip(kpl, pk, err):
        RS_VS_KPL[filename][_kpl] = (_pk, _err)
        dsum[_kpl] = dsum.get(_kpl, 0) + _pk / _err**2
        dwgt[_kpl] = dwgt.get(_kpl, 0) + 1 / _err**2
    k0 = np.abs(kpl).argmin()
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
        RS_VS_KPL['total'][_kpl] = (dsum[_kpl] / dwgt[_kpl], 1./np.sqrt(dwgt[_kpl]))
    for _kpl in dsum_fold:
        RS_VS_KPL_FOLD['total'][_kpl] = (dsum_fold[_kpl] / dwgt_fold[_kpl], 1./np.sqrt(dwgt_fold[_kpl]))

fig=p.figure(figsize=(12,7.2))
if False: # put in raw delay spec
    f = np.load(dspec)
    kpl,pk,err = f['kpl'], f['pk'], f['err']
    if True: #if 'I' in sep: # Add foregrounds
        for cnt,k in enumerate(kpl):
            k = '%6.3f' % k
            if not FG_VS_KPL.has_key(k): continue
            pk[cnt] += FG_VS_KPL[k]
            err[cnt] = 2*np.sqrt(FG_VS_KPL_NOS*FG_VS_KPL[k])
#    pk *= .76 # psa747 calibration of Pic A = 370.6 Jy @ 160 MHz (which includes resolution effects)
#    err *= .76
    #pk *= 2.35 # Use power**2 beam, which is a 1.69/0.72=2.35 penalty factor
    #err *= 2.35
#    pk *= 1.46 # fix scale to reflect new corrected bandwidth normalization
#    err *= 1.46
    pk *= 1.4
    err *= 1.4
    if True: # For aggressive fringe-rate filtering, change beam area
        pk *= 1.90 # ratio of power**2 beams for filtered * unfiltered beams: 0.306 / 0.162
        err *= 1.90
    print pk.shape, err.shape
    print pk + 2*err
    dual_plot(kpl, pk, 2*err, color='c', upperlimit=True) # 2-sigma error bars



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
        #print dsum[k] / dwgt[k], 1./np.sqrt(dwgt[k])
    ks_fold = RS_VS_KPL_FOLD[sep].keys(); ks_fold.sort()
    for k in ks_fold:
        _d,_n = RS_VS_KPL_FOLD[sep][k]
        #print k, _d, _n, '->',
        dsum_fold[k] = dsum_fold.get(k,0) + _d / _n**2
        dwgt_fold[k] = dwgt_fold.get(k,0) + 1 / _n**2
        #print dsum[k] / dwgt[k], 1./np.sqrt(dwgt[k])
        
    kpl = dsum.keys(); kpl.sort()
    d = [dsum[k]/dwgt[k] for k in kpl]
    nos = [1./np.sqrt(dwgt[k]) for k in kpl]
    kpl_fold = dsum_fold.keys(); kpl_fold.sort()
    d_fold = [dsum_fold[k]/dwgt_fold[k] for k in kpl_fold]
    nos_fold = [1./np.sqrt(dwgt_fold[k]) for k in kpl_fold]
    #d = [RS_VS_KPL[k][0] for k in kpl]
    #nos = [RS_VS_KPL[k][1] for k in kpl]
    if True: #if 'I' in sep: # Add foregrounds
        print 'Adding Foregrounds'
        for cnt,k in enumerate(kpl):
            k = '%6.3f' % k
            if not FG_VS_KPL.has_key(k): continue
            d[cnt] += FG_VS_KPL[k]
            nos[cnt] = 2*np.sqrt(FG_VS_KPL_NOS*FG_VS_KPL[k])
        for cnt,k in enumerate(kpl_fold):
            k = '%6.3f' % k
            if not FG_VS_KPL.has_key(k): continue
            d_fold[cnt] += FG_VS_KPL[k]
            nos_fold[cnt] = 2*np.sqrt(FG_VS_KPL_NOS*FG_VS_KPL[k])
        
    d,kpl,nos = np.array(d, dtype=np.complex), np.array(kpl), np.array(nos)
    d_fold,kpl_fold,nos_fold = np.array(d_fold, dtype=np.complex), np.array(kpl_fold), np.array(nos_fold)
    if False:
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
    if False:
        f = 2.35 # Use power**2 beam, which is a 1.69/0.72=2.35 penalty factor
        print 'Scaling data and noise by %f for correcting cosmo scalar to use power^2 beam' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    if True: # For aggressive fringe-rate filtering, change beam area
        f = opts.afrf_factor
#       f = 1.90 # ratio of power**2 beams for filtered * unfiltered beams: 0.306 / 0.162
        f = 1.39 # ratio of power**2 beams for filtered * unfiltered beams: 0.306 / 0.162
        print 'Scaling data and noise by %f for beam constriction in aggressive fringe-rate filtering.' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    if True: # extra penalty for signal loss in covariance diagonalization
        #f = 1.2
#        f = 1.15 # this is the conservative number, for the biggest modes outside the of the wedge. Get ~1.03 correction at k=.2
        f = 1.02 # this is the conservative number, for the biggest modes outside the of the wedge. Get ~1.03 correction at k=.2
        print 'Scaling data and noise by %f for signal loss from empirically estimating covariances.' % f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    if True:
        f = 1/np.log(2)
        print 'Scaling data and noise by %f for using the median statistic.'%f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f
    if True:
        f1 = 1.048
        f2 = 1.013
        f  = 1.0015
        print 'Scaling data and noise by %f(first mode outside horizon), %f(second mode outside horizon), and %f(all modes greater)  for using delay filter.'%(f1,f2,f)
        k0 = np.argmin(np.abs(kpl))
        d[k0+2] *= f1; d[k0-2] *= f1
        d[k0+3] *= f2; d[k0-3] *= f2
        d[k0+4:] *= f; d[:k0-3] *= f
        nos[k0+2] *= f1; nos[k0-2] *= f1
        nos[k0+3] *= f2; nos[k0-3] *= f2
        nos[k0+4:] *= f; nos[:k0-3] *= f

        d_fold[2] *= f1; nos_fold *= f1
        d_fold[3] *= f2; nos_fold *= f2
        d_fold[4:] *= f; nos_fold *= f
    if True:
        f=1 + 2e-9 
        print 'Scaling data and noise by %f for signal loss from bandpass calibrationp.'%f
        d *= f
        nos *= f
        d_fold *= f
        nos_fold *= f

    #for _kpl,_pk,_nos in zip(kpl,d,nos): print _kpl, _pk, _nos
    print sep, colors[0]
    '''
    if True: # Hacky way to get a noise bias out, if necessary
        d -= np.median(np.concatenate([d[:8], d[-8:]]))
    if True: # Hacky way to estimate noise
        nos = np.std(np.concatenate([d[:8], d[-8:]])) * np.ones_like(d)
    '''
    if d_fold.size == 0: d_fold,nos_fold = None, None
    #dual_plot(kpl, d, 2*nos, d_fold, 2*nos_fold, color=colors[0], bins=BINS)#,f0=freq) # 2-sigma error bars
    dual_plot(kpl, d, 2*nos, d_fold, 2*nos_fold, color=colors[0], bins=BINS,f0=.151) # 2-sigma error bars
    #dual_plot(kpl, d, nos, color=colors[0], bins=BINS) # 2-sigma error bars
    colors = colors[1:] + colors[0]

def mean_temp(z):
    return 28. * ((1.+z)/10.)**.5 # mK

import glob
re_z = re.compile(r'power_21cm_z(\d+\.\d+)\.dat')
kpl_pos = ks[k0+1:]
for filename in glob.glob('lidz_mcquinn_k3pk/*7.3*dat'):
    print 'Reading', filename
    d = np.array([map(float, L.split()) for L in open(filename).readlines()])
    ks, pk = d[:,0], d[:,1]
    z_file = float(re_z.match(os.path.basename(filename)).groups()[0])
    z = C.pspec.f2z(.151)
    k3pk = ks**3 / (2*np.pi**2) * pk
    p.subplot(122)
    p.plot(ks, k3pk * mean_temp(z)**2, 'm-')
    
tau_h = 100 + 15. #in ns
k_h = C.pspec.dk_deta(C.pspec.f2z(.151)) * tau_h
p.subplot(121)
p.vlines(k_h, -1e7, 1e8, linestyles='--', linewidth=1.5)
p.vlines(-k_h, -1e7, 1e8, linestyles='--', linewidth=1.5)
#p.vlines(k_h, 1, 1e18, linestyles='--', linewidth=1.5)
#p.vlines(-k_h, 1, 1e18, linestyles='--', linewidth=1.5)
#p.gca().set_yscale('log', nonposy='clip')
p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
p.ylabel(r'$P(k)[\ {\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$',fontsize='large')
p.ylim(-.6e7,1.75e7)
#p.ylim(1e5,5e16)
p.grid()


p.subplot(122)
#if ONLY_POS_K: p.plot([.5], [248**2], 'mv', label='GMRT2013')
#else: p.plot([-.5, .5], [248**2, 248**2], 'mv', label='GMRT2013')
p.vlines(k_h, 1, 1e17, linestyles='--', linewidth=1.5)
#theoretical_ks = np.linspace(.058,.5, 100)
#theor_errs = 1441090 * np.array(theoretical_ks)**3  / (2*np.pi**2)
#p.plot(theoretical_ks, theor_errs, 'c--')

#GMRT RESULTS
#p.plot([.1,.13,.17,.19,.27,.31,.4,.5,.6][:-1], [2e5,4e5,1e5,2e5,2e5,4e5,5e5,248**2,3e5][:-1], 'yv', label='GMRT2013')
p.plot([.5], 248**2, 'yv', label='GMRT2013', markersize=9)
#Dillon upper limit 9.82e7 mK^2 at k=0.2
p.plot([.2], 9.82e7 * .2**3/(2*np.pi**2), 'mv', label='Dillon2013', markersize=9)
#Parsons2014 upper lmit 41mk at .27
p.plot([.27], 41**2, 'gv', label='Parsons2014', markersize=9)

theo_noise = noise_level()
#print k0
#print kpl_pos[0], theo_noise[0]
#2 for the 2 sigma
p.plot(np.array(kpl_pos), 2*np.array(kpl_pos)**3*theo_noise/(2*np.pi**2), 'c--')
p.gca().set_yscale('log', nonposy='clip')
p.xlabel(r'$k\ [h\ {\rm Mpc}^{-1}]$', fontsize='large')
p.ylabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$', fontsize='large')
p.ylim(1e0,1e5)
p.xlim(0, 0.6)
p.grid()
#p.savefig('pspec_greece.png')

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

f = np.load(args[0])
def posterior(kpl, pk, err, pkfold=None, errfold=None, f0=.151, umag=16.):
    import scipy.interpolate as interp
    k0 = np.abs(kpl).argmin()
    kpl = kpl[k0:]
    z = C.pspec.f2z(f0)
    kpr = C.pspec.dk_du(z) * umag
    k = np.sqrt(kpl**2 + kpr**2)
    if pkfold is None:
        print 'Folding for posterior'
        pkfold = pk[k0:].copy()
        errfold = err[k0:].copy()
        pkpos,errpos = pk[k0+1:].copy(), err[k0+1:].copy()
        pkneg,errneg = pk[k0-1:0:-1].copy(), err[k0-1:0:-1].copy()
        pkfold[1:] = (pkpos/errpos**2 + pkneg/errneg**2) / (1./errpos**2 + 1./errneg**2)
        errfold[1:] = np.sqrt(1./(1./errpos**2 + 1./errneg**2))

    #ind = np.logical_and(kpl>.2, kpl<.5)
    ind = np.logical_and(k>.15, k<.5)
    #ind = np.logical_and(kpl>.12, kpl<.5)
    #print kpl,pk.real,err
    k = k[ind]
    pkfold = pkfold[ind]
    errfold = errfold[ind]
#    if True:
    if False:
        #remove k=.345 point
        w = np.floor(k*100)!=34
        k = k[w]
        pkfold=pkfold[w]
        errfold = errfold[w]
    pk= k**3 * pkfold/(2*np.pi**2)
    err = k**3 * errfold/(2*np.pi**2)
    err_omit = err.copy()
    err_omit[3] = 1e10 # give no weight to this point
    #s = np.logspace(1,3.5,100)
    s = np.linspace(-5000,5000,10000)
#    print s
    data = []
    data_omit = []
    for _k, _pk, _err in zip(k, pk, err):
        print _k, _pk.real, _err
    #    print '%6.3f    %9.5f     9.5f'%(_k, _pk.real, _err)
    for ss in s:
        data.append(np.exp(-.5*np.sum((np.abs(pk) - ss)**2 / err**2)))
        data_omit.append(np.exp(-.5*np.sum((np.abs(pk) - ss)**2 / err_omit**2)))
    #    print data[-1]
    data = np.array(data)
    data_omit = np.array(data_omit)
    #print data
    #print s
    #data/=np.sum(data)
    data /= np.max(data)
    data_omit /= np.max(data_omit)
    p.figure(5, figsize=(6.5,5.5))
    p.plot(s, data, 'k', linewidth=2)
#    p.plot(s, data_omit, 'k--', linewidth=1)
    #use a spline interpolator to get the 1 and 2 sigma limits.
    #spline = interp.interp1d(data,s)
    #print spline
    #print max(data), min(data)
    #print spline(.68), spline(.95)
    #p.plot(spline(np.linspace(.0,1,100)),'o')
#    p.plot(s, np.exp(-.5)*np.ones_like(s))
#    p.plot(s, np.exp(-.5*2**2)*np.ones_like(s))
    data_c = np.cumsum(data)
    data_omit_c = np.cumsum(data_omit)
    data_c /= data_c[-1]
    data_omit_c /= data_omit_c[-1]
    mean = s[np.argmax(data)]
    s1lo,s1hi = s[data_c<0.1586][-1], s[data_c>1-0.1586][0]
    s2lo,s2hi = s[data_c<0.0227][-1], s[data_c>1-0.0227][0]
    print 'Posterior: Mean, (1siglo,1sighi), (2siglo,2sighi)'
    print 'Posterior:', mean, (s1lo,s1hi), (s2lo,s2hi)
    mean_o = s[np.argmax(data_omit)]
    s1lo_o,s1hi_o = s[data_omit_c<0.1586][-1], s[data_omit_c>1-0.1586][0]
    s2lo_o,s2hi_o = s[data_omit_c<0.0227][-1], s[data_omit_c>1-0.0227][0]
    print 'Posterior (omit):', mean_o, (s1lo_o,s1hi_o), (s2lo_o,s2hi_o)
    #sig1 = []
    #sig2 = []
    #s1 = data[maxarg:] - np.exp(-.5)
    #sig1.append(np.floor(np.median(np.where(np.abs(s1)<.01)))+maxarg)
    #s1 = data[:maxarg] - np.exp(-.5)
    #sig1.append(np.floor(np.median(np.where(np.abs(s1)<.01))))

    #s2 = data[maxarg:] - np.exp(-.5*2**2)
    #sig2.append(np.floor(np.median(np.where(np.abs(s2)<.01)))+maxarg)
    #s2 = data[:maxarg] - np.exp(-.5*2**2)
    #sig2.append(np.floor(np.median(np.where(np.abs(s2)<.01))))
    
    #p.vlines(s[sig1[-1]],0,1,color=(0,107/255.,164/255.), linewidth=2)
    p.vlines(s1lo,0,1,color=(0,107/255.,164/255.), linewidth=2)
    p.vlines(s1hi,0,1,color=(0,107/255.,164/255.), linewidth=2)
    #p.vlines(s1lo_o,0,1,color=(0,107/255.,164/255.), linestyle='--', linewidth=2)
    #p.vlines(s1hi_o,0,1,color=(0,107/255.,164/255.), linestyle='--', linewidth=2)

    #p.vlines(s[sig2[-1]],0,1,color=(0,107/255.,164/255.), linewidth=2)
    # limits for data_omit
    p.vlines(s2lo,0,1,color=(1,128/255.,14/255.), linewidth=2)
    p.vlines(s2hi,0,1,color=(1,128/255.,14/255.), linewidth=2)
    #p.vlines(s2lo_o,0,1,color=(1,128/255.,14/255.), linestyle='--', linewidth=2)
    #p.vlines(s2hi_o,0,1,color=(1,128/255.,14/255.), linestyle='--', linewidth=2)

    p.xlabel(r'$k^3/2\pi^2\ P(k)\ [{\rm mK}^2]$', fontsize='large')
    p.ylabel('Posterior Distribution', fontsize='large')
    p.xlim(0,700)
    p.grid(1)
    p.subplots_adjust(left=.15, top=.95, bottom=.15, right=.95)
    p.show()

#posterior(f['kpl'], f['pk'], f['err'], f['pk_fold'], f['err_fold'])
posterior(kpl, d, 2*nos, d_fold, nos_fold)



if opts.show:
    p.show()
