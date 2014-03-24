#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C, scipy.interpolate
import re, sys

def mean_temp(z):
    return 28. * ((1.+z)/10.)**.5 # mK

re_z = re.compile(r'power_21cm_z(\d+\.\d+)\.dat')
npz = n.load('../fg_vs_umag_vs_fq.npz')

DIM,RES = 200, .5
SIZE = DIM / RES
aa = a.cal.get_aa('psa746_v008', n.array([.150]))
im = a.img.Img(size=DIM, res=RES)
L,M = im.get_LM((SIZE/2,SIZE/2))
lmag = L[SIZE/2]
tx,ty,tz = im.get_top((SIZE/2,SIZE/2))
tx,ty,tz = tx.flatten(), ty.flatten(), tz.flatten()
bm_resp = aa[0].bm_response((tx,ty,tz), pol='y')**2
bm_resp /= aa[0].bm_response((0,0,1), pol='y')**2
bm_resp = n.where(tz <= 0, 0, bm_resp)
bm_resp.shape = (SIZE, SIZE)

ax1 = p.subplot(211)
p.setp(ax1.get_xticklabels(), visible=False)
ax2 = p.subplot(212, sharex=ax1)

dat = {}
for cnt, filename in enumerate(sys.argv[1:]):
    print 'Reading', filename
    d = n.array([map(float, L.split()) for L in open(filename).readlines()])
    ks, pk = d[:,0], d[:,1]
    #z = float(re_z.match(filename).groups()[0])
    z = 9.0
    #k3pk = ks**3 / (2*n.pi**2) * pk * mean_temp(z)**2
    pk *= mean_temp(z)**2
    interp = scipy.interpolate.interp1d(ks, pk, kind='cubic',
        bounds_error=False, fill_value=1e-6)
    #dat[filename] = (ks, k3pk)
    dat[filename] = interp

    krange = 10**n.arange(n.log10(.065), n.log10(9.5), .01)
    umags = [16,32,64,128]
    for umag,color in zip(umags,'bgrc'):
        islope = -C.pspec.dk_deta(z) * umag / .150
        ks = islope * lmag
        kresp = bm_resp.sum(axis=0) / n.sum(bm_resp)
        pk = []
        for k in krange:
            _k = n.abs(k + ks)
            pk.append(n.sum(dat[filename](_k) * kresp))
        pk = n.array(pk)
        k3pk = krange**3 / (2*n.pi**2) * pk
        p.subplot(211)
        p.loglog(krange, k3pk, color+'-')
        if cnt == 0:
            kresp /= kresp.max()
            p.subplot(212)
            p.loglog(n.abs(ks+0.25), kresp, color+':')
            p.loglog(n.abs(ks+0.50), kresp, color+'--')
            p.loglog(n.abs(ks+1.00), kresp, color+'-')

    pk = n.array(dat[filename](krange))
    k3pk = krange**3 / (2*n.pi**2) * pk
    p.subplot(211)
    p.loglog(krange, k3pk, 'k-')
#p.xlim(2e-2,2)
#p.ylim(1e-1,2e0)
p.subplot(211)
p.ylabel(r'$\Delta^2(k)$')
p.xlim(.06, 3)
p.grid()

p.subplot(212)
p.xlabel(r'$k_\parallel\ h\ {\rm Mpc}^{-1}$')
p.ylabel('Normalized Response')
p.xlim(.06, 3)
p.ylim(.1,2)
p.grid()

p.show()
