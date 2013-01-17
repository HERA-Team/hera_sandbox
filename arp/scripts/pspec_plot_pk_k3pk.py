#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse

def dual_plot(kpl, pk, err, umag=16., f0=.164, color='', bins=None):
    z = C.pspec.f2z(f0)
    kpr = C.pspec.dk_du(z) * umag
    k = n.sqrt(kpl**2 + kpr**2)
    k3 = n.abs(k**3 / (2*n.pi**2))
    print 'k [h Mpc^-1], P(k) [K^2], err (2sigma)'
    for _k,_pk,_err in zip(kpl,pk,err):
        print '%6.3f, %9.5f, %9.5f' % (_k, _pk/1e6, _err/1e6)
    print '-'*20
    p.subplot(121)
    p.errorbar(kpl, pk, yerr=err, fmt=color+'.-')
    p.subplot(122)
    #p.errorbar(kpl, k3*pk, yerr=k3*err, fmt=color+'.-')
    if not bins is None:
        _kpls, _k3pks, _k3errs = [], [], []
        for (dn,up) in bins:
            kavg = n.sqrt((0.5*(up+dn))**2 + kpr**2)
            dsum,dwgt = 0., 0.
            for _kpl,_pk,_err in zip(kpl, pk, err):
                if dn < _kpl and _kpl <= up:
                    #print _pk, _err
                    dsum += _pk / _err**2
                    dwgt += 1. / _err**2
            if dwgt == 0: continue
            _pk = dsum / dwgt
            #print _pk
            _err = 1. / n.sqrt(dwgt)
            _k3pk = kavg**3/(2*n.pi**2) * _pk
            _k3err = kavg**3/(2*n.pi**2) * _err
            _kpls.append(dn); _kpls.append(up)
            _k3pks.append(_k3pk); _k3pks.append(_k3pk)
            _k3errs.append(_k3err); _k3errs.append(_k3err)
        kpl = n.array(_kpls)
        k3pk = n.array(_k3pks)
        k3err = n.array(_k3errs)
    else:
        kpl = kpl
        k3pk = k3*pk
        k3err = k3*err
    #p.plot(kpl, k3*pk+k3*err, color+'.-')
    p.plot(kpl, k3pk+k3err, color+'.-')

o = optparse.OptionParser()
opts,args = o.parse_args(sys.argv[1:])

FG_VS_KPL_NOS = 168.74e6
FG_VS_KPL = { # K^2
    '-0.054':   5.37262e+13,
    '-0.027':   7.15304e+14, 
    ' 0.000':   3.50958e+15, 
    ' 0.027':   4.12396e+14, 
    ' 0.054':   2.60795e+13,
    #-0.0536455587089:   5.37262e+13,
    #-0.0268227793545:   7.15304e+14, 
    #0.0:                3.50958e+15, 
    #0.0268227793545:    4.12396e+14, 
    #0.0536455587089:    2.60795e+13,
}

RS_VS_KPL = {} # K^2
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    RS_VS_KPL[filename] = {}
    for _kpl, _pk, _err in zip(f['kpl'], f['pk'], f['err']):
        RS_VS_KPL[filename][_kpl] = (_pk, _err)

BINS = None
BINS = ((-.52,-.4),(-.4,-.2),(-.2,-.1),(-.1,-.05),(-.05,.05),(.05,.1),(.1,.2),(.2,.4),(.4,.52))
#BINS = ((-.6,-.33),(-.33,-.22),(-.22,-.1),(-.1,.1),(.1,.22),(.22,.33),(.33,.6))
#BINS = ((-.6,-.5),(-.5,-.4),(-.4,-.3),(-.3,-.2),(-.2,-.1),(-.1,-.05),(-.05,.05),(.05,.1),(.1,.2),(.2,.3),(.3,.4),(.4,.5),(.5,.6))
#kpl = RS_VS_KPL.keys(); kpl.sort()
#d = [RS_VS_KPL[k] for k in kpl]
colors = 'kbcm'
for sep in RS_VS_KPL:
    dsum, dwgt = {}, {}
    ks = RS_VS_KPL[sep].keys(); ks.sort()
    for k in ks:
        _d,_n = RS_VS_KPL[sep][k]
        #print k, _d, _n, '->',
        dsum[k] = dsum.get(k,0) + _d / _n**2
        dwgt[k] = dwgt.get(k,0) + 1 / _n**2
        #print dsum[k] / dwgt[k], 1./n.sqrt(dwgt[k])
        
    kpl = dsum.keys(); kpl.sort()
    d = [dsum[k]/dwgt[k] for k in kpl]
    nos = [1./n.sqrt(dwgt[k]) for k in kpl]
    #d = [RS_VS_KPL[k][0] for k in kpl]
    #nos = [RS_VS_KPL[k][1] for k in kpl]
    if 'I' in sep: # Add foregrounds
        for cnt,k in enumerate(kpl):
            k = '%6.3f' % k
            if not FG_VS_KPL.has_key(k): continue
            d[cnt] += FG_VS_KPL[k]
            nos[cnt] = 2*n.sqrt(FG_VS_KPL_NOS*FG_VS_KPL[k])
    d,kpl,nos = n.array(d, dtype=n.complex), n.array(kpl), n.array(nos)
    # Currently calibrated to Pictor A @ 160 MHz = 424 Jy
    # To recalibrate to new Pic A, must multiply by square of ratio of fluxes
    #d *= 1.448 # Recalibrate to new Pic A spec from Jacobs 12/21/12?  How well do we know this?
    #d *= 0.774 # Recalibrate to Pic A from Perley et al. 1997
    #d *= 1.125 # Recalibrate to Pic A from Slee 1995
    d *= .65 # current best guess by psa747 data
    d *= 2.35 # Use power**2 beam
    #for _kpl,_pk,_nos in zip(kpl,d,nos): print _kpl, _pk, _nos
    print sep, colors[0]
    dual_plot(kpl, d, 2*nos, color=colors[0], bins=BINS) # 2-sigma error bars
    colors = colors[1:]
    
p.subplot(121)
p.gca().set_yscale('log', nonposy='clip')
p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
p.ylabel(r'$P(k)\ [{\rm mK}^2\ (h^{-1}\ {\rm Mpc})^3]$')
p.ylim(1e5,3e16)
p.grid()
p.subplot(122)
p.gca().set_yscale('log', nonposy='clip')
p.xlabel(r'$k_\parallel\ [h\ {\rm Mpc}^{-1}]$')
p.ylabel(r'$k^3/2\pi\ P(k)\ [{\rm mK}^2]$')
p.ylim(1e1,1e7)
p.grid()
p.show()
