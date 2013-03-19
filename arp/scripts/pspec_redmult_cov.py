#! /usr/bin/env python
import aipy as a, numpy as n
import pylab as p
import capo
import optparse, sys, os

def quick_diag_dot(A, B):
    return n.array([n.dot(A[...,i,:], B[...,:,i]) for i in range(A.shape[-2])])
def quick_trace_dot(A, B):
    return quick_diag_dot(A,B).sum()


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

def bl_index(bl):
    i,j = a.miriad.bl2ij(bl)
    return i * 32 + j

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
etas = n.fft.fftshift(capo.pspec.f2eta(afreqs))
kpl = etas * capo.pspec.dk_deta(z)

if PFB:
    # XXX unsure how much of a BW modification a windowed PFB needs.  I think not much...
    B = sdf * afreqs.size / NTAPS
else:
    B = sdf * afreqs.size / capo.pfb.NOISE_EQUIV_BW[WINDOW]
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

T, W = {}, {}

times = []
for filename in args:
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, opts.ant, opts.pol)

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
        if PFB:
            _Trms = capo.pfb.pfb(Trms, window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
            _Wrms = capo.pfb.pfb(w   , window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
        else:
            window = a.dsp.gen_window(Trms.size, WINDOW)
            _Trms = n.fft.ifft(window * Trms)
            _Wrms = n.fft.ifft(window * w)
        gain = n.abs(_Wrms[0])
        #print 'Gain:', gain
        if True: # we think clean messes up noise statistics.  Do we need it? (Yes)
            if gain > 0:
                if False:
                    _Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-2)
                    #print info['term']
                    #_Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9)
                    _Trms = _Tcln + info['res'] / gain
                else:
                    _Trms.shape = (_Trms.size,1)
                    C = n.zeros((_Trms.size, _Trms.size), dtype=n.complex)
                    for k1 in xrange(_Wrms.size):
                      for k2 in xrange(_Wrms.size):
                        #C[k1,k2] = _Wrms[k2-k1]
                        C[k1,k2] = _Wrms[k1-k2]
                    _C = n.linalg.inv(C)
                    #p.subplot(131); capo.arp.waterfall(C, mode='log', drng=2)
                    #p.subplot(132); capo.arp.waterfall(_C, mode='log', drng=2)
                    #p.subplot(133); p.plot(_Trms, 'k')
                    _Trms_ = n.dot(_C, _Trms)
                    #p.plot(_Trms_, 'g')
                    #_Tcln, info = a.deconv.clean(_Trms.squeeze(), _Wrms, tol=1e-9)
                    #_Trms__ = _Tcln + info['res'] / gain
                    #p.plot(_Trms__, 'b')
                    #p.show()
                    _Trms = _Trms_.squeeze()
        else: _Trms /= gain
        _Trms = n.fft.fftshift(_Trms)
        _Wrms = n.fft.fftshift(_Wrms)
        if False: # swap in a simulated signal
            _Trms = n.random.normal(size=_Trms.size) * n.exp(2j*n.pi*n.random.uniform(size=_Trms.size))
            #for x in range(_Trms.size):
            #    _Trms[x] += 2 * n.exp(10j*(10+x)*t)
            _Trms[26] += 2 * n.exp(100j*t)
            #_Trms[5] += 2 * n.exp(1000j*t)

        T[bl] = T.get(bl, []) + [_Trms]
        W[bl] = W.get(bl, []) + [_Wrms]


def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

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

def gen_Q(knum, n_k):
    Q = n.zeros((n_k,n_k))
    Q[knum,knum] = 1
    return Q


n_k = chans.size
bls = T.keys()
Ts = n.concatenate([T[bl] for bl in bls], axis=-1).T
print Ts.shape
if False:
    Ts = n.zeros((Ts.shape[0], Ts.shape[1]*4), dtype=n.complex)
    for t in xrange(Ts.shape[1]):
        _Trms = n.random.normal(size=Ts.shape[0]) * n.exp(2j*n.pi*n.random.uniform(size=Ts.shape[0]))
        _Trms[26] += 2 * n.exp(.1j*t)
        _Trms[52] += 2 * n.exp(.1j*t)
        Ts[:,t] = _Trms
    

print ' '.join(['%d_%d' % a.miriad.bl2ij(bl) for bl in bls])
C = CoV(Ts, bls)
capo.arp.waterfall(C.C, mode='log', drng=2); p.show()

Q = {}
for k in range(n_k): Q[k] = gen_Q(k,n_k)

pspecs = []
for cnt,bl1 in enumerate(bls):
    for bl2 in bls[cnt:]:
        if bl2 == bl1: continue # skip auto-products
        print a.miriad.bl2ij(bl1), a.miriad.bl2ij(bl2)

        x1 = C.get_x(bl1)
        x2 = C.get_x(bl2)
        C1 = C.get_C(bl1)
        C2 = C.get_C(bl2)
        if True:
            d = n.diag(C1)
            C1 /= n.sqrt(n.multiply.outer(d,d)) * 2
            d = n.diag(C2)
            C2 /= n.sqrt(n.multiply.outer(d,d)) * 2
        SZ = C1.shape[0]

        if False: # Artificially diagonalize covariance matrix
            mask = n.eye(SZ)
            #mask[1:] += n.eye(SZ)[:-1]
            #mask[:-1] += n.eye(SZ)[1:]
            #mask[2:] += n.eye(SZ)[:-2]
            #mask[:-2] += n.eye(SZ)[2:]
            mask[17] = 1
            mask[:,17] = 1
            C1 *= mask
            C2 *= mask

        _C1 = n.linalg.inv(C1)
        _C2 = n.linalg.inv(C2)

        #norm = n.sum(_C1, axis=-1); norm.shape = (norm.size,1)
        norm = n.sum(_C1, axis=0); norm.shape = (1,norm.size)
        p.subplot(121); capo.arp.waterfall(_C1, mode='log', drng=2)
        _C1 /= norm
        p.subplot(122); capo.arp.waterfall(_C1, mode='log', drng=2); p.show()
        z1 = n.dot(_C1, x1)
        z2 = n.dot(_C2, x2)
        F1 = _C1
        F2 = _C2

        if False:
            w,v = n.linalg.eig(n.linalg.inv(F1))
            M1 = n.dot(v, n.dot(n.diagflat(n.sqrt(w)), v.T))
            w,v = n.linalg.eig(n.linalg.inv(F2))
            M2 = n.dot(v, n.dot(n.diagflat(n.sqrt(w)), v.T))
        else:
            #M1 = n.linalg.inv(F1)
            #M2 = n.linalg.inv(F2)
            M1 = n.eye(F1.shape[0])
            M2 = n.eye(F2.shape[0])
        
        W1 = n.dot(M1,F1)
        #W1 = n.dot(M1,n.abs(F1)**2)
        norm = n.sum(W1, axis=-1); norm.shape = (norm.size,1)
        #norm = n.sum(W1, axis=0); norm.shape = (1,norm.size)
        M1 /= norm
        W1 = n.dot(M1,F1)
        W2 = n.dot(M2,F2)
        #W2 = n.dot(M2,n.abs(F2)**2)
        norm = n.sum(W2, axis=-1); norm.shape = (norm.size,1)
        #norm = n.sum(W2, axis=0); norm.shape = (1,norm.size)
        M2 /= norm
        W2 = n.dot(M2,F2)

        p1 = n.dot(M1, z1)
        p2 = n.dot(M2, z2)


        x1_,x2_ = x1.copy(),x2.copy()
        #for k in xrange(17,24):
        for k in xrange(40):
        #for k in [17,23]:
            c1 = C1[:,k:k+1].copy(); c1[k] = 0
            c2 = C2[:,k:k+1].copy(); c2[k] = 0
            #p.plot(n.average(x1_*x1_.conj(), axis=1), 'k')
            x1__ = x1[k:k+1] * c1
            x2__ = x2[k:k+1] * c2
            #x1__ = x1[k:k+1] * c1 / n.sqrt(n.diag(C1) * C1[k,k])
            #p.plot(n.average(x1__ * x1__.conj(), axis=1), 'g')
            x1_ -= x1__
            #p.plot(n.average(x1_*x1_.conj(), axis=1), 'r')
            #p.show()
            x2_ -= x2__
        p.subplot(131)
        p.plot(scalar*n.average(x1_*x1_.conj(), axis=1).real, 'k')
        p.plot(scalar*n.average(x1*x1.conj(), axis=1).real, 'r')
        p.plot(scalar*n.average(z1*z1.conj(), axis=1).real, 'g')
        p.subplot(132)
        p.plot(scalar*n.average(x2_*x2_.conj(), axis=1).real, 'k')
        p.plot(scalar*n.average(x2*x2.conj(), axis=1).real, 'r')
        p.plot(scalar*n.average(p2*p2.conj(), axis=1).real, 'g')
        p.subplot(133)
        p.plot(scalar*n.average(x1_*x2_.conj(), axis=1).real, 'k')
        p.plot(scalar*n.average(x1*x2.conj(), axis=1).real, 'r')
        p.plot(scalar*n.average(p1*p2.conj(), axis=1).real, 'g')
        p.show()
        continue
        #p.subplot(231); capo.arp.waterfall(x1, mode='log', mx=.5, drng=2); p.colorbar(shrink=.5)
        #p.subplot(232); capo.arp.waterfall(x1-x1[17:18,:]*C1[:,17:18], mode='log', mx=.5, drng=2)
        #p.subplot(233); capo.arp.waterfall(x1-x1[17:18,:]*C1[:,17:18].conj(), mode='log', mx=.5, drng=2)
        #p.subplot(234); capo.arp.waterfall(x2, mode='log', drng=2), p.colorbar(shrink=.5)
        #p.subplot(235); capo.arp.waterfall(x2-x2[17:18,:]*C2[:,17:18], mode='log', mx=.5, drng=2)
        #p.subplot(236); capo.arp.waterfall(x2-x2[17:18,:]*C2[:,17:18].conj(), mode='log', mx=.5, drng=2)
        #p.show()
        #w,v = n.linalg.eig(C21)
        #print n.abs(w)

        p.subplot(221); capo.arp.waterfall(C21, mode='log', drng=2); p.colorbar(shrink=.5)
        p.subplot(222); capo.arp.waterfall(C12, mode='log', drng=2); p.colorbar(shrink=.5)
        p.subplot(223); capo.arp.waterfall(_C21, mode='log', drng=2); p.colorbar(shrink=.5)
        p.subplot(224); capo.arp.waterfall(_C12, mode='log', drng=2); p.colorbar(shrink=.5)
        p.show()
    
        #p.subplot(131); capo.arp.waterfall(C21, mode='log', drng=3); p.colorbar(shrink=.5)
        #p.subplot(132); capo.arp.waterfall(_C21, mode='log', drng=3); p.colorbar(shrink=.5)
        #p.subplot(133); capo.arp.waterfall(_C21T, mode='log', drng=3); p.colorbar(shrink=.5)
        #p.show()

        #print n.all(_C21 == _C12)
    
        p.subplot(241); capo.arp.waterfall(x1, mode='log', drng=2)
        p.subplot(242); capo.arp.waterfall(_C21, mode='log', drng=2)
        p.subplot(243); capo.arp.waterfall(z1, mode='log', drng=2)
        p.subplot(244); capo.arp.waterfall(p1, mode='log', drng=2)
        p.subplot(245); capo.arp.waterfall(x2, mode='log', drng=2)
        p.subplot(246); capo.arp.waterfall(_C12, mode='log', drng=2)
        p.subplot(247); capo.arp.waterfall(z2, mode='log', drng=2)
        p.subplot(248); capo.arp.waterfall(p2, mode='log', drng=2)
        p.show()

        p.plot(scalar*n.average(p1*p2.conj(), axis=1).real, 'k')
        p.plot(scalar*n.average(x1*x2.conj(), axis=1).real, 'r')
        p.show()

        print 'Making Fisher/Gisher Matrix'
        F,G = n.zeros((n_k,n_k), dtype=n.complex), n.zeros((n_k,n_k), dtype=n.complex)

        q, Q_C = [], {}
        Q_C12,Q_C21 = {},{}
        for k in range(n_k):
            q.append(0.5 * quick_diag_dot(z1.T.conj(), n.dot(Q[k], z2)))
            Q_C[k] = n.dot(Q[k], _C12)
            Q_C12[k] = n.dot(Q[k], _C12)
            Q_C21[k] = n.dot(Q[k], _C21)
        q = n.array(q)
        #p.subplot(131); capo.arp.waterfall(x1, mode='log', drng=2)
        #p.subplot(132); capo.arp.waterfall(z1, mode='log', drng=2)
        #p.subplot(133); capo.arp.waterfall(q, mode='log', drng=2)
        #p.show()
        if True:
            for k1 in range(n_k):
                for k2 in range(n_k):
                    #F[k1,k2] = 0.5 * quick_trace_dot(Q_C[k1],Q_C[k2])
                    F[k1,k2] = 0.5 * quick_trace_dot(Q_C21[k1],Q_C12[k2])
                    #F[k1,k2] = 0.5 * quick_trace_dot(Q_C21[k2],Q_C12[k1])
                    #F[k1,k2] = 0.5 * quick_trace_dot(Q_C12[k1],Q_C21[k2])
                    #F[k1,k2] = 0.5 * quick_trace_dot(Q_C12[k2],Q_C21[k1])
        else:
            F = cov(q)
        G = cov(q)

        print 'Making Normalization/Windowing'
        #_G = n.linalg.inv(G)
        _G = n.linalg.inv(F)
        if True:# Set M = G^-1/2
            w,v = n.linalg.eig(_G)
            M = n.dot(v, n.dot(n.diagflat(n.sqrt(w)), v.T))
        else:
            M = _G
        # Normalize M s.t. rows of W sum to 1
        W = n.dot(M,F)
        norm = n.sum(W, axis=-1); norm.shape = (norm.size,1)
        M /= norm
        #print n.diag(M)
        #capo.arp.waterfall(_G, mode='log', drng=4); p.show()
        W = n.dot(M,F)
        #print 'Wdiag', W.diagonal()
        #p.subplot(111); capo.arp.waterfall(W, mode='log', drng=2); p.show()
        pk = scalar * n.dot(M,q)
        p.subplot(131); capo.arp.waterfall(q, mode='real'); p.colorbar(shrink=.5)
        p.subplot(132); capo.arp.waterfall(M, mode='log', drng=2)
        p.subplot(133); capo.arp.waterfall(pk, mode='real'); p.colorbar(shrink=.5)
        p.show()
        
        pk_avg = n.average(pk, axis=1)
        pspecs.append(pk_avg)
        p.plot(pk_avg.real, 'k')
        p.plot(scalar*n.average(x1*x2.conj(), axis=1).real, 'r')
        # XXX still have to get Sigma
        p.show()
pspecs = n.array(pspecs)
avg_1d = n.average(pspecs, axis=0)
#std_1d = n.std(pspecs, axis=0) / n.sqrt(pspecs.shape[0])
std_1d = n.ones_like(avg_1d)

print 'Writing pspec.npz'
n.savez('pspec.npz', kpl=kpl, pk=avg_1d, err=std_1d)
p.show()

