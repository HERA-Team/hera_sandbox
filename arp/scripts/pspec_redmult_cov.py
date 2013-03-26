#! /usr/bin/env python
import aipy as a, numpy as n
import pylab as p
import capo
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True)
o.add_option('-t', '--taps', type='int', default=1,
    help='Taps to use in the PFB.  Default 1, which instead uses windowed FFT')
opts,args = o.parse_args(sys.argv[1:])

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
for filename in args:
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, opts.ant, opts.pol)
    curtime = [None]
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        if len(times) == 0 or times[-1] != t: times.append(t)
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
            if t != curtime[-1] or False: # Set to true for independent (i.e. thermal) noise for each bl
                curtime.append(t)
                if len(curtime) % 8 == 2 or False: # Set False for shared (i.e. eor) noise that matches fringe rate
                    Trms_ = n.random.normal(size=Trms.size) * n.exp(2j*n.pi*n.random.uniform(size=Trms.size))
                    Trms_ *= TSYS / n.sqrt(B * T_INT * NDAY * NBL * NPOL)
                    #Trms_ *= n.sqrt(n.sqrt(351./43)) # penalize for oversampling fr-filtered data
                    Trms_ *= 1.14 # adjust to suboptimal flux calibration
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
        T[bl] = T.get(bl, []) + [_Trms]
        N[bl] = N.get(bl, []) + [_Nrms]

n_k = chans.size
bls = T.keys()
Ts = n.concatenate([T[bl] for bl in bls], axis=-1).T
Ns = n.concatenate([N[bl] for bl in bls], axis=-1).T
print Ts.shape
print ' '.join(['%d_%d' % a.miriad.bl2ij(bl) for bl in bls])
#capo.arp.waterfall(cov(Ts), mode='log', drng=2); p.show()

X = Ts.copy()
N = Ns.copy() # this noise copy processed as if it were the data
_Cxtot,_Cntot = 1., 1
PLT1,PLT2 = 3,3
for cnt in xrange(PLT1*PLT2-1):
    print cnt, '/', PLT1*PLT2-1
    p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(cov(X), mode='log', mx=0, drng=3)
    #p.subplot(PLT1,PLT2,i+1); capo.arp.waterfall(cov(N2), mode='log', mx=0, drng=2)
    SZ = X.shape[0]
    Cx,Cn = cov(X), cov(N)
    for c in [Cx,Cn]: # Normalize covariance matrices
        d = n.diag(c); d.shape = (1,SZ)
        c /= n.sqrt(d) * 2
    #g = 1. / len(bls)
    g = 1
    _Cx,_Cn = -g*Cx, -g*Cn
    ind = n.arange(SZ)
    for b in xrange(len(bls)): # zero out redundant off-diagonals
        indb = ind[:-b*n_k]
        _Cx[indb,indb+b*n_k] = _Cx[indb+b*n_k,indb] = 0
        _Cn[indb,indb+b*n_k] = _Cn[indb+b*n_k,indb] = 0
    _Cx[ind,ind] = _Cn[ind,ind] = 0 # set these to zero temporarily to avoid noise bias into cross terms
    for _C in [_Cx,_Cn]: # subtract off covariances common to all baselines pairs
        # XXX this loop may be introducing a small noise bias by injecting a common component in each
        # square of the covariance matrix, and then subtracting this common cov * same baseline
        _C.shape = (len(bls),n_k,len(bls),n_k)
        sub_C = n.zeros_like(_C)
        for i in xrange(len(bls)):
            for j in xrange(len(bls)):
                not_ij = n.array([bl for bl in xrange(len(bls)) if not bl in [i,j]])
                _C_avg = n.mean(n.mean(_C.take(not_ij,axis=0).take(not_ij,axis=2), axis=2), axis=0)
                sub_C[i,:,j,:] = _C_avg
        _C.shape = (len(bls)*n_k,len(bls)*n_k)
        sub_C.shape = (len(bls)*n_k,len(bls)*n_k)
        _C -= sub_C
    _Cx[ind,ind] = _Cn[ind,ind] = 1
    if True: # for first cycles, only use auto-products to avoid noise contamination that causes bias
        _Cx.shape = _Cn.shape = (len(bls),n_k,len(bls),n_k)
        for i in xrange(len(bls)):
            for j in xrange(len(bls)):
                # XXX need to make this symmetric for baselines at beginning and end.
                if abs(i - j) <= cnt: continue
                _Cx[i,:,j,:] = _Cn[i,:,j,:] = 0
        _Cx.shape = _Cn.shape = (len(bls)*n_k,len(bls)*n_k)
    X,N = n.dot(_Cx,X), n.dot(_Cx,N)
    _Cxtot,_Cntot = n.dot(_Cx,_Cxtot), n.dot(_Cn,_Cntot)
## Normalize to maintain amplitude assuming all k modes are independent & same amplitude
#norm1 = n.sqrt(n.sum(n.abs(_C1tot)**2, axis=1)); norm1.shape = (norm1.size,1)
#x1 /= norm1; x2 /= norm2
p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(X), mode='log', mx=0, drng=3)
#p.subplot(PLT1,PLT2,i+2); capo.arp.waterfall(cov(N2), mode='log', mx=0, drng=2)
p.show()

p.subplot(221); capo.arp.waterfall(_Cxtot, mode='log', drng=2)
p.subplot(222); capo.arp.waterfall(cov(Ts), mode='log', drng=2)
p.subplot(223); capo.arp.waterfall(cov(n.dot(_Cxtot,Ts)), mode='log', drng=2)
p.subplot(224); capo.arp.waterfall(cov(X), mode='log', drng=2)
p.show()

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

pspecs = []
Cx,Cn = CoV(Ts, bls), CoV(Ns, bls)
Cx_ = CoV(X, bls)
Cn1_,Cn2_ = CoV(N1, bls), CoV(N2, bls)
for cnt,bli in enumerate(bls):
    for blj in bls[cnt:]:
        if bli == blj: continue
        print a.miriad.bl2ij(bli), a.miriad.bl2ij(blj)
        xi,xj = Cx.get_x(bli), Cx.get_x(blj)
        xi_,xj_ = Cx_.get_x(bli), Cx_.get_x(blj)
        ni,nj = Cn.get_x(bli), Cn.get_x(blj)
        n1i_,n1j_ = Cn1_.get_x(bli), Cn1_.get_x(blj)
        n2i_,n2j_ = Cn2_.get_x(bli), Cn2_.get_x(blj)
        nij = n.sqrt(n.mean(n.abs(ni*nj.conj())**2, axis=1))
        n1ij_ = n.sqrt(n.mean(n.abs(n1i_*n1j_.conj())**2, axis=1))
        n2ij_ = n.sqrt(n.mean(n.abs(n2i_*n2j_.conj())**2, axis=1))
        f1 = n.sqrt(n.mean(n.abs(nij)**2)/n.mean(n.abs(n1ij_)**2))
        f2 = n.sqrt(n.mean(n.abs(nij)**2)/n.mean(n.abs(n2ij_)**2))
        print 'Rescale factor:', f1, f2
        fudge = max(f1,f2)

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
            p.plot(n.average(n2i_*n2i_.conj(), axis=1).real, 'k')
            p.plot(n.average(n2j_*n2j_.conj(), axis=1).real, 'b')
            p.plot(n.average(n2i_*n2j_.conj(), axis=1).real, 'g')
            p.plot(n.average(ni*nj.conj(), axis=1).real, 'r')
            p.plot(fudge*n.average(n2i_*n2j_.conj(), axis=1).real, 'c')
            p.show()
        elif True:
            p.subplot(121); p.plot(n.average(xi*xj.conj(), axis=1).real)
            p.subplot(122); p.plot(n.average(xi_*xj_.conj(), axis=1).real)
        pk_avg = scalar * n.average(xi_ * xj_.conj(), axis=1) * fudge # XXX
        pspecs.append(pk_avg)
p.show()
pspecs = n.array(pspecs)
avg_1d = n.average(pspecs, axis=0)
#std_1d = n.std(pspecs, axis=0) / n.sqrt(pspecs.shape[0])
std_1d = n.std(pspecs, axis=0) # in new noise subtraction, this remaining dither is essentially a bootstrap error, but with 5/7 of the data

print 'Writing pspec.npz'
n.savez('pspec.npz', kpl=kpl, pk=avg_1d, err=std_1d)
p.show()

