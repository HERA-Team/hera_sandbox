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
#capo.arp.waterfall(C.C, mode='log', drng=2); p.show()

Q = {}
for k in range(n_k): Q[k] = gen_Q(k,n_k)

pspecs = []
for cnt,bl1 in enumerate(bls):
    for bl2 in bls[cnt:]:
        if bl2 == bl1: continue # skip auto-products
        print a.miriad.bl2ij(bl1), a.miriad.bl2ij(bl2)

        x1,x2 = C.get_x(bl1), C.get_x(bl2)
        SZ = x1.shape[0]
        #n1 = n.random.normal(size=SZ) * n.exp(2j*n.pi*n.random.uniform(size=SZ))
        n1 = n.exp(2j*n.pi*n.random.uniform(size=SZ)); n1.shape = (SZ,1)
        n2 = n.exp(2j*n.pi*n.random.uniform(size=SZ)); n2.shape = (SZ,1)
        _C1tot,_C2tot = n.eye(SZ), n.eye(SZ) # these aren't actually the inverses, just for norm purposes
        for i in xrange(8):
            X = n.concatenate([x1,x2], axis=0)
            #p.subplot(3,3,i+1); capo.arp.waterfall(cov(X), mode='log', mx=0, drng=2)
            p.subplot(3,3,i+1); capo.arp.waterfall(cov(X), mode='log', drng=2)
            C1,C2 = cov(x1), cov(x2)
            C12,C21 = cov2(x1,x2), cov2(x2,x1)
            #C1,C2 = C.get_C(bl1), C.get_C(bl2)
            if True:
                d1 = n.diag(C1); d1.shape = (1,SZ)
                C1 /= n.sqrt(d1) * 2
                d2 = n.diag(C2); d2.shape = (1,SZ)
                C2 /= n.sqrt(d2) * 2
                C12 /= n.sqrt(d2) * 2
                C21 /= n.sqrt(d1) * 2

            if False: _C1,_C2 = n.linalg.inv(C1), n.linalg.inv(C2)
            else:
                _C1,_C2 = -C1,-C2
                _C12,_C21 = -C12,-C21
                #_C1,_C2 = n.zeros_like(C1), n.zeros_like(C2)
                #_C12,_C21 = n.zeros_like(C12), n.zeros_like(C21)
                ind = n.arange(SZ)
                _C1[ind,ind] = _C2[ind,ind] = 1
                _C12[ind,ind] = _C21[ind,ind] = 0

            x1_,x2_ = n.dot(_C1,x1), n.dot(_C2,x2)
            x1__,x2__ = n.dot(_C12,x2), n.dot(_C21,x1)
            x1,x2 = x1_+x1__, x2_+x2__
            n1_,n2_ = n.dot(_C1,n1), n.dot(_C2,n2)
            n1__,n2__ = n.dot(_C12,n2), n.dot(_C21,n1)
            n1,n2 = n1_+n1__, n2_+n2__
            _C1tot,_C2tot = n.dot(_C1+_C12, _C1tot), n.dot(_C2+_C21, _C2tot)
        norm1 = n.sqrt(n.sum(n.abs(_C1tot)**2, axis=1)); norm1.shape = (norm1.size,1)
        norm2 = n.sqrt(n.sum(n.abs(_C2tot)**2, axis=1)); norm2.shape = (norm2.size,1)
        x1 /= norm1; x2 /= norm2
        n1 /= norm1; n2 /= norm2
        print n.abs(n1).flatten()
        print n.abs(n2).flatten()
        #print norm1.flatten()
        #print norm2.flatten()
        X = n.concatenate([x1,x2], axis=0)
        p.subplot(3,3,i+2); capo.arp.waterfall(cov(X), mode='log', mx=0, drng=2)
        p.show()
        pk_avg = scalar * n.average(x1 * x2.conj(), axis=1)

        x1orig,x2orig = C.get_x(bl1), C.get_x(bl2)
        p.subplot(131)
        p.plot(scalar*n.average(x1*x1.conj(), axis=1).real, 'k')
        p.plot(scalar*n.average(x1orig*x1orig.conj(), axis=1).real, 'r')
        p.subplot(132)
        p.plot(scalar*n.average(x2*x2.conj(), axis=1).real, 'k')
        p.plot(scalar*n.average(x2orig*x2orig.conj(), axis=1).real, 'r')
        p.subplot(133)
        p.plot(scalar*n.average(x1*x2.conj(), axis=1).real, 'k')
        p.plot(scalar*n.average(x1orig*x2orig.conj(), axis=1).real, 'r')
        p.show()
        pspecs.append(pk_avg)
pspecs = n.array(pspecs)
avg_1d = n.average(pspecs, axis=0)
std_1d = n.std(pspecs, axis=0) / n.sqrt(pspecs.shape[0])
#std_1d = n.ones_like(avg_1d)

print 'Writing pspec.npz'
n.savez('pspec.npz', kpl=kpl, pk=avg_1d, err=std_1d)
p.show()

