#! /usr/bin/env python
import aipy as a, numpy as n
import pylab as p
import capo
import optparse, sys, os, random

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True)
o.add_option('-t', '--taps', type='int', default=1,
    help='Taps to use in the PFB.  Default 1, which instead uses windowed FFT')
opts,args = o.parse_args(sys.argv[1:])

PLOT = False

NTAPS = opts.taps
if NTAPS > 1: PFB = True
else: PFB = False
WINDOW = 'blackman-harris'
#WINDOW = 'none'

# XXX Currently hardcoded for PSA898
A_ = [0,16,8,24,4,20,12,28]
B_ = [i+1 for i in A_]
C_ = [i+2 for i in A_]
D_ = [i+3 for i in A_]
ANTPOS = n.array([A_, B_, C_, D_])

class CoV:
    def __init__(self, X, bls):
        self.bls = bls
        self.X = X
        self.nprms = X.shape[0] / len(bls)
        self.C = cov(X)
    def get_C(self, bl1, bl2=None):
        if bl2 is None: bl2 = bl1
        i,j = self.bls.index(bl1), self.bls.index(bl2)
        return self.C[i*self.nprms:(i+1)*self.nprms, j*self.nprms:(j+1)*self.nprms].copy()
    def get_Cinv(self, bl1, bl2=None):
        if bl2 is None: bl2 = bl1
        return n.linalg.inv(self.get_C(bl1,bl2))
    def get_x(self, bl):
        i = self.bls.index(bl)
        return self.X[i*self.nprms:(i+1)*self.nprms].copy()
    def get_Ck(self, k):
        inds = n.arange(len(self.bls)) * self.nprms + k
        return self.C.take(inds,axis=0).take(inds,axis=1)

def bl_index(bl):
    i,j = a.miriad.bl2ij(bl)
    return i * 32 + j

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

def cov2(m1,m2):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X1 = n.array(m1, ndmin=2, dtype=n.complex)
    X2 = n.array(m2, ndmin=2, dtype=n.complex)
    X1 -= X1.mean(axis=1)[(slice(None),n.newaxis)]
    X2 -= X2.mean(axis=1)[(slice(None),n.newaxis)]
    N = X1.shape[1]
    fact = float(N - 1)
    return (n.dot(X1, X2.T.conj()) / fact).squeeze()

def gen_Q(knum, n_k, dim):
    Q = n.zeros((dim,dim), dtype=n.float)
    for i in n.arange(0,Q.shape[0],n_k):
        for j in n.arange(0,Q.shape[1],n_k):
            Q[i+knum,j+knum] = 1
    if False: # Remove auto-products from estimator to avoid +noise bias
        Q[n.arange(Q.shape[0]),n.arange(Q.shape[1])] = 0
    return Q

def quick_diag_dot(A, B):
    return n.array([n.dot(A[...,i,:], B[...,:,i]) for i in range(A.shape[-2])])
def quick_trace_dot(A, B):
    return quick_diag_dot(A,B).sum()


# Get a dict of all separations and the bls that contribute
bl2sep = {}
sep2bl = {}
for ri in range(ANTPOS.shape[0]):
    for ci in range(ANTPOS.shape[1]):
        for rj in range(ANTPOS.shape[0]):
            for cj in range(ci,ANTPOS.shape[1]):
                if ri >= rj and ci == cj: continue # exclude repeat +/- listings of certain bls
                #sep = a.miriad.ij2bl(rj-ri, cj-ci)
                sep = a.miriad.ij2bl(cj-ci, rj-ri) # prefer to have column as leading term to orient E/W baselines
                i,j = ANTPOS[ri,ci], ANTPOS[rj,cj]
                bl = a.miriad.ij2bl(i,j)
                if i > j: i,j,sep = j,i,-sep
                bl2sep[bl] = sep
                sep = n.abs(sep)
                sep2bl[sep] = sep2bl.get(sep,[]) + [bl]

uv = a.miriad.UV(args[0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
del(uv)

afreqs = freqs.take(chans)
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

if PFB:
    # XXX unsure how much of a BW modification a windowed PFB needs.  I think not much...
    B = sdf * afreqs.size / NTAPS
    etas = n.fft.fftshift(capo.pspec.f2eta(afreqs[:afreqs.size/NTAPS]))
else:
    B = sdf * afreqs.size / capo.pfb.NOISE_EQUIV_BW[WINDOW]
    etas = n.fft.fftshift(capo.pspec.f2eta(afreqs))
kpl = etas * capo.pspec.dk_deta(z)
bm = n.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq)
scalar = capo.pspec.X2Y(z) * bm * B
#scalar = 1
print 'Freq:',fq
print 'z:', z
print 'B:', B
print 'scalar:', scalar

#cen_fqs = n.arange(.115,.190,.005)
#cen_fqs = n.array([.150])
#kwargs = {'cen_fqs':cen_fqs,'B':B, 'ntaps':NTAPS, 'window':WINDOW, 'bm_fqs':afreqs.clip(.120,.190)}
#window = a.dsp.gen_window(freqs.size, window=WINDOW)

T, N = {}, {}
times = []
eor_mdl = {}
for filename in args:
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, opts.ant, opts.pol)
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        if len(times) == 0 or times[-1] != t:
            if len(times) % 8 == 0:
                eor_mdl[t] = n.random.normal(size=chans.size) * n.exp(2j*n.pi*n.random.uniform(size=chans.size))
            else: eor_mdl[t] = eor_mdl[times[-1]]
            times.append(t)
        bl = a.miriad.ij2bl(i,j)
        sep = bl2sep[bl]
        if sep < 0:
            #print 'Conj:', a.miriad.bl2ij(bl)
            d,sep = n.conj(d),-sep
        d,f = d.take(chans), f.take(chans)
        w = n.logical_not(f).astype(n.float)
        Trms = d * capo.pspec.jy2T(afreqs)
        if True: # generate noise
            TSYS = 560e3 # mK
            B = 100e6 / uvi['nchan']
            NDAY = 44
            NBL = 4
            NPOL = 2
            T_INT = 43. # for just compressed data
            #T_INT = 351. # for fringe-rate filtered data
            Trms_ = n.random.normal(size=Trms.size) * n.exp(2j*n.pi*n.random.uniform(size=Trms.size))
            Trms_ *= TSYS / n.sqrt(B * T_INT * NDAY * NBL * NPOL)
            #Trms_ *= n.sqrt(n.sqrt(351./43)) # penalize for oversampling fr-filtered data
            Trms_ *= 1.14 # adjust to suboptimal flux calibration
            #Trms_ += eor_mdl[t]
            Trms_ = eor_mdl[t]
            Nrms  = Trms_ * w
        if PFB:
            _Trms = capo.pfb.pfb(Trms, window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
            _Nrms = capo.pfb.pfb(Nrms, window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
            _Wrms = capo.pfb.pfb(w   , window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
        else:
            window = a.dsp.gen_window(Trms.size, WINDOW)
            _Trms = n.fft.ifft(window * Trms)
            _Nrms = n.fft.ifft(window * Nrms)
            _Wrms = n.fft.ifft(window * w)
        gain = n.abs(_Wrms[0])
        #print 'Gain:', gain
        if gain > 0:
            _Trms.shape = (_Trms.size,1)
            C = n.zeros((_Trms.size, _Trms.size), dtype=n.complex)
            for k1 in xrange(_Wrms.size):
              for k2 in xrange(_Wrms.size):
                #C[k1,k2] = _Wrms[k2-k1]
                C[k1,k2] = _Wrms[k1-k2]
            _C = n.linalg.inv(C)
            _Trms = n.dot(_C, _Trms).squeeze()
            _Nrms = n.dot(_C, _Nrms).squeeze()
        _Trms = n.fft.fftshift(_Trms)
        _Nrms = n.fft.fftshift(_Nrms)
        if False: # swap in a simulated signal post delay transform
            _Trms = n.random.normal(size=_Trms.size) * n.exp(2j*n.pi*n.random.uniform(size=_Trms.size))
            mask = n.ones(_Trms.size); mask[15:25] = 0
            _Trms += .3*eor_mdl[times[-1]] * mask
        T[bl] = T.get(bl, []) + [_Trms]
        N[bl] = N.get(bl, []) + [_Nrms]

n_k = chans.size
bls = T.keys()
for bl in bls: T[bl],N[bl] = n.array(T[bl]),n.array(N[bl])
if False:
    print 'Fringe-rate filtering the noise to match the data'
    for bl in N:
        _N = n.fft.ifft(N[bl], axis=0)
        _N[23:] = 0 # This was calculated by hand for fr-filter with max_fr=1. and min_fr=0.
        N[bl] = n.fft.fft(_N, axis=0)
if False:
    print 'Adding extra noise into the data'
    for bl in bls: T[bl] += N[bl]
Ts = n.concatenate([T[bl] for bl in bls], axis=-1).T
Ns = n.concatenate([N[bl] for bl in bls], axis=-1).T
if False:
    print 'Switching sign of alternate integrations to decorrelate sky'
    sign = 1
    for i in xrange(Ts.shape[1]):
        if i % 8 == 0: sign = -sign
        Ts[:,i] *= sign
        Ns[:,i] *= sign
if False:
    print 'Switching sign of various baselines & modes to decorrelate sky'
    for i in xrange(Ts.shape[0]):
        if n.random.uniform() > .5:
            Ts[i] *= -1
            Ns[i] *= -1

print Ts.shape
print ' '.join(['%d_%d' % a.miriad.bl2ij(bl) for bl in bls])
#capo.arp.waterfall(cov(Ts), mode='log', drng=2); p.show()
p.subplot(131); capo.arp.waterfall(Ts, mode='log', mx=1, drng=2); p.colorbar(shrink=.5)
p.subplot(132); capo.arp.waterfall(Ns, mode='log', mx=1, drng=2); p.colorbar(shrink=.5)
p.subplot(133); capo.arp.waterfall(cov(Ts), mode='log', drng=3); p.colorbar(shrink=.5)
p.show()

for boot in xrange(20):
    if True: # pick a sample of baselines with replacement
        #bls_ = [random.choice(bls) for bl in bls]
        bls_ = random.sample(bls, len(bls))
        gp1,gp2 = bls_[:len(bls)/2],bls_[len(bls)/2:] # ensure gp1 and gp2 can't share baselines
        # ensure each group has at least 2 kinds of baselines
        gp1 = random.sample(gp1, 2) + [random.choice(gp1) for bl in gp1[:len(gp1)-2]]
        gp2 = random.sample(gp2, 2) + [random.choice(gp2) for bl in gp2[:len(gp2)-2]]
    else:
        gp1,gp2 = bls[:len(bls)/2],bls[len(bls)/2:]
    bls_ = gp1 + gp2
    print 'Bootstrap sample %d:' % boot,
    for gp in [gp1,gp2]: print '(%s)' % (','.join(['%d_%d' % a.miriad.bl2ij(bl) for bl in gp])),
    print
    Ts = n.concatenate([T[bl] for bl in bls_], axis=1).T
    Ns = n.concatenate([N[bl] for bl in bls_], axis=1).T # this noise copy processed as if it were the data
    L = len(bls_)

    _Cxtot,_Cntot = 1, 1
    #PLT1,PLT2 = 3,3
    #PLT1,PLT2 = 2,2
    PLT1,PLT2 = 2,3
    for cnt in xrange(PLT1*PLT2-1):
        #print cnt, '/', PLT1*PLT2-1
        if PLOT:
            p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(cov(Ts), mode='log', mx=0, drng=3)
            ##p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(cov(Ns), mode='log', mx=0, drng=2)
        SZ = Ts.shape[0]
        Cx,Cn = cov(Ts), cov(Ns)
        for c in [Cx,Cn]: # Normalize covariance matrices
            d = n.diag(c); d.shape = (1,SZ)
            c /= n.sqrt(d) * 2
        g = .3
        _Cx,_Cn = -g*Cx, -g*Cn
        if False: # restrict to certain modes in the covariance diagonalization
            mask = n.zeros_like(Cx)
            for k in xrange(17,24):
                for b in xrange(L):
                    mask[b*n_k+k] = mask[:,b*n_k+k] = 1
            _Cx *= mask; _Cn *= mask
        ind = n.arange(SZ)
        for b in xrange(L): # zero out redundant off-diagonals
            indb = ind[:-b*n_k]
            _Cx[indb,indb+b*n_k] = _Cx[indb+b*n_k,indb] = 0
            _Cn[indb,indb+b*n_k] = _Cn[indb+b*n_k,indb] = 0
        _Cx[ind,ind] = _Cn[ind,ind] = 0 # set these to zero temporarily to avoid noise bias into cross terms
        if True: # remove covariances common to all bl pairs.  XXX This is responsible for >75% of noise bias
            for _C in [_Cx,_Cn]:
                _C.shape = (L,n_k,L,n_k)
                sub_C = n.zeros_like(_C)
                #for i in xrange(L):
                #    for j in xrange(L):
                #        not_ij = n.array([bl for bl in xrange(L) if not bl in [i,j]])
                #        _C_avg = n.mean(n.mean(_C.take(not_ij,axis=0).take(not_ij,axis=2), axis=2), axis=0)
                #        #capo.arp.waterfall(_C_avg, mode='log', drng=2); p.colorbar(shrink=.5); p.show()
                #        sub_C[i,:,j,:] = _C_avg
                for i in xrange(L):
                    bli = bls_[i]
                    for j in xrange(L):
                        blj = bls_[j]
                        if bli in gp1 and blj in gp1: gp = gp1
                        elif bli in gp2 and blj in gp2: gp = gp2
                        else: continue # make sure we only compute average using bls in same group
                        _Csum,_Cwgt = 0,0
                        # XXX as constructed, this will explode if a group consists entirely of one bl.
                        for bli_ in gp:
                            i_ = bls_.index(bli_)
                            if i_ == i: continue
                            for blj_ in gp:
                                j_ = bls_.index(blj_)
                                if j_ == j: continue
                                _Csum += _C[i_,j_]
                                _Cwgt += 1
                        sub_C[i,j] = _Csum / _Cwgt # XXX careful if _Cwgt is 0
                _C.shape = sub_C.shape = (L*n_k,L*n_k)
                #p.clf()
                #p.subplot(131); capo.arp.waterfall(_C, mode='log', mx=0, drng=2)
                #p.subplot(132); capo.arp.waterfall(sub_C, mode='log', mx=0, drng=2)
                _C -= sub_C
                #p.subplot(133); capo.arp.waterfall(_C, mode='log', mx=0, drng=2)
                #p.show()
        if True: # divide baselines into two independent groups to avoid cross-contamination of noise
            mask = n.ones_like(Cx)
            for bl1 in xrange(len(gp1)):
                for bl2 in xrange(len(gp2)):
                    bl2 += len(gp1)
                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            #capo.arp.waterfall(mask, mode='real'); p.show()
            _Cx *= mask; _Cn *= mask
        _Cx[ind,ind] = _Cn[ind,ind] = 1
        Ts,Ns = n.dot(_Cx,Ts), n.dot(_Cn,Ns)
        _Cxtot,_Cntot = n.dot(_Cx,_Cxtot), n.dot(_Cn,_Cntot)
    if PLOT:
        p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ts), mode='log', mx=0, drng=3)
        ##p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ns), mode='log', mx=0, drng=2)
        p.show()

    #p.subplot(221); capo.arp.waterfall(_Cxtot, mode='log', drng=2)
    #p.subplot(222); capo.arp.waterfall(cov(Ts), mode='log', drng=2)
    #p.subplot(223); capo.arp.waterfall(cov(n.dot(_Cxtot,Ts)), mode='log', drng=2)
    #p.subplot(224); capo.arp.waterfall(cov(X), mode='log', drng=2)
    #p.show()
    Ts = n.concatenate([T[bl] for bl in bls_], axis=1).T
    Ns = n.concatenate([N[bl] for bl in bls_], axis=1).T # this noise copy processed as if it were the data

    pspecs,dspecs = [], []
    Cx,Cn = CoV(Ts, bls_), CoV(Ns, bls_)
    Cx_ = CoV(n.dot(_Cxtot,Ts), bls_)
    Cn1_,Cn2_ = CoV(n.dot(_Cntot,Ns), bls_), CoV(n.dot(_Cxtot,Ns), bls_)
    for cnt,bli in enumerate(bls_):
        for blj in bls_[cnt:]:
            if bli == blj: continue
            if True: # exclude intra-group pairings
                if (bli in gp1 and blj in gp1) or (bli in gp2 and blj in gp2): continue
            print a.miriad.bl2ij(bli), a.miriad.bl2ij(blj)
            # XXX behavior here is poorly defined for repeat baselines in bootstrapping
            xi,xj = Cx.get_x(bli), Cx.get_x(blj)
            xi_,xj_ = Cx_.get_x(bli), Cx_.get_x(blj)
            if True: # do an extra final removal of leakage from particular modes
                Ts = n.concatenate([xi_,xj_], axis=0)
                cx = cov(Ts)
                if PLOT:
                    p.clf()
                    p.subplot(121); capo.arp.waterfall(cx, mode='log', mx=0, drng=3)
                for cnt1 in xrange(9):
                    d = n.diag(cx); d.shape = (1,d.size); cx /= n.sqrt(d) * 2
                    g = .3
                    _cx = -g*cx
                    mask = n.zeros_like(cx)
                    for k in xrange(17,24):
                        mask[k] = mask[:,k] = 1
                        mask[k+n_k] = mask[:,k+n_k] = 1
                    ind = n.arange(n_k)
                    #mask[ind,ind] = mask[ind+n_k,ind+n_k] = 1 # don't need this b/c _cx gets set to 1 on diag
                    mask[ind,ind+n_k] = mask[ind+n_k,ind] = 0
                    _cx *= mask
                    _cx[ind,ind] = _cx[ind+n_k,ind+n_k] = 1
                    #p.subplot(132); capo.arp.waterfall(_cx, mode='log', mx=0, drng=3)
                    Ts = n.dot(_cx, Ts)
                    cx = cov(Ts)
                if PLOT:
                    p.subplot(122); capo.arp.waterfall(cx, mode='log', mx=0, drng=3)
                    p.show()
                xi_,xj_ = Ts[:n_k],Ts[n_k:]
            ni,nj = Cn.get_x(bli), Cn.get_x(blj)
            n1i_,n1j_ = Cn1_.get_x(bli), Cn1_.get_x(blj)
            n2i_,n2j_ = Cn2_.get_x(bli), Cn2_.get_x(blj)
            nij = n.sqrt(n.mean(n.abs(ni*nj.conj())**2, axis=1))
            n1ij_ = n.sqrt(n.mean(n.abs(n1i_*n1j_.conj())**2, axis=1))
            n2ij_ = n.sqrt(n.mean(n.abs(n2i_*n2j_.conj())**2, axis=1))
            #f1 = n.sqrt(n.mean(n.abs(nij)**2)/n.mean(n.abs(n1ij_)**2))
            #f2 = n.sqrt(n.mean(n.abs(nij)**2)/n.mean(n.abs(n2ij_)**2))
            f1 = n.sqrt((n.abs(nij)**2)/(n.abs(n1ij_)**2))
            f2 = n.sqrt((n.abs(nij)**2)/(n.abs(n2ij_)**2))
            print 'Rescale factor:', f1, f2
            #fudge = max(f1,f2)
            fudge = 1.

            if False:
                p.subplot(221)
                p.plot(n.average(xi_*xi_.conj(), axis=1).real, 'k')
                p.plot(n.average(xi*xi.conj(), axis=1).real, 'r')
                p.subplot(222)
                p.plot(n.average(xj_*xj_.conj(), axis=1).real, 'k')
                p.plot(n.average(xj*xj.conj(), axis=1).real, 'r')
                p.subplot(223)
                p.plot(n.average(xi_*xj_.conj(), axis=1).real, 'k')
                p.plot(n.average(xi*xj.conj(), axis=1).real, 'r')
                p.plot(fudge*n.average(xi_*xj_.conj(), axis=1).real, 'g')
                p.subplot(224)
                p.plot(n.average(n1i_*n1i_.conj(), axis=1).real, 'k')
                p.plot(n.average(n1j_*n1j_.conj(), axis=1).real, 'b')
                p.plot(n.average(n1i_*n1j_.conj(), axis=1).real, 'g')
                p.plot(n.average(ni*nj.conj(), axis=1).real, 'r')
                p.plot(fudge*n.average(n1i_*n1j_.conj(), axis=1).real, 'c')
                p.show()
            elif False:
                p.subplot(131); p.plot(n.average(xi*xj.conj(), axis=1).real)
                p.subplot(132); p.plot(n.average(xi_*xj_.conj(), axis=1).real)
            pk_avg = scalar * n.average(xi * xj.conj(), axis=1)
            pk_avg_ = scalar * n.average(xi_ * xj_.conj(), axis=1) * fudge # XXX
            dspecs.append(pk_avg)
            pspecs.append(pk_avg_)
    pspecs,dspecs = n.array(pspecs), n.array(dspecs)
    #avg_1d = n.average(pspecs, axis=0)
    avg_1d = n.average(dspecs, axis=0)
    #p.subplot(133)
    p.plot(avg_1d.real,'.')
    #p.plot(n.average(dspecs, axis=0).real/scalar)
    #p.show()
    std_1d = n.std(pspecs, axis=0) / n.sqrt(pspecs.shape[0])
    #std_1d = n.std(pspecs, axis=0) # in new noise subtraction, this remaining dither is essentially a bootstrap error, but with 5/7 of the data

    print 'Writing pspec_boot%04d.npz' % boot
    n.savez('pspec_boot%04d.npz'%boot, kpl=kpl, pk=avg_1d, err=std_1d)
p.show()

import sys; sys.exit(0)
print 'Making Fisher Matrix'
Qs,Q_C = {}, {}
for k in range(n_k):
    Qs[k] = gen_Q(k,n_k,_Cxtot.shape[0])
    Q_C[k] = n.dot(Qs[k], _Cxtot)
F = n.zeros((40,40), dtype=n.complex)

for k1 in range(n_k):
  for k2 in range(n_k):
    F[k1,k2] = 0.5 * quick_trace_dot(Q_C[k1],Q_C[k2])
_F = n.linalg.inv(F)

print 'Making Normalization/Windowing'
# Set M = F^-1/2
w,v = n.linalg.eig(_F)
M = n.dot(v, n.dot(n.diagflat(n.sqrt(w)), v.T))
# Normalize M s.t. rows of W sum to 1
W = n.dot(M,F)
norm = n.sum(W, axis=-1); norm.shape = (norm.size,1)
M /= norm
W = n.dot(M,F) # this is just to prove normalization

p.subplot(221); capo.arp.waterfall(F, mode='log', drng=2)
p.subplot(222); capo.arp.waterfall(_F, mode='log', drng=2)
p.subplot(223); capo.arp.waterfall(M, mode='log', drng=2)
p.subplot(224); capo.arp.waterfall(W, mode='log', drng=2)
p.show()

print 'Minding ps and qs'
qs = []
for k in range(n_k):
    print k
    _CQ_C = n.dot(_Cxtot, Q_C[k])
    qs.append(0.5 * quick_diag_dot(Ts.T.conj(), n.dot(_CQ_C, Ts)))
qs = n.array(qs)
print qs.shape, M.shape
ps = n.dot(M,qs)
covp = n.dot(W,M.T)

p.subplot(131); capo.arp.waterfall(qs, mode='log', drng=2)
p.subplot(132); capo.arp.waterfall(ps, mode='log', drng=2)
p.subplot(133); capo.arp.waterfall(covp, mode='log', drng=2); p.colorbar(shrink=.5)
p.show()

p.subplot(121); p.plot(scalar * n.average(qs, axis=1))
p.subplot(122); p.plot(scalar * n.average(ps, axis=1))
p.show()

