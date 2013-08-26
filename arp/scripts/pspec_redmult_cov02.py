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

if PFB:
    # XXX unsure how much of a BW modification a windowed PFB needs.  I think not much...
    B = sdf * afreqs.size / NTAPS
else:
    B = sdf * afreqs.size / capo.pfb.NOISE_EQUIV_BW[WINDOW]
bm = n.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq)
scalar = capo.pspec.X2Y(z) * bm * B
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
                #_Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9, maxiter=100, stop_if_div=False, verbose=False)
                #_Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9)
                _Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-2)
                #print info['term']
                #_Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9)
                _Trms = _Tcln + info['res'] / gain
        else: _Trms /= gain
        _Trms = n.fft.fftshift(_Trms)
        _Wrms = n.fft.fftshift(_Wrms)
        if True: # swap in a simulated signal
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
        #self.C = n.cov(X)
        self.C = cov(X)
        self.Cinv = {}
        self.zs = {}
    def get_C(self, bl1, bl2=None):
        if bl2 is None: bl2 = bl1
        i,j = self.bls.index(bl1), self.bls.index(bl2)
        return self.C[i*self.nprms:(i+1)*self.nprms, j*self.nprms:(j+1)*self.nprms]
    def get_Cinv(self, bl1, bl2=None):
        if bl2 is None: bl2 = bl1
        if bl1 == bl2:
            Cinv = self.Cinv.get(bl1, n.linalg.inv(self.get_C(bl1)))
            self.Cinv[bl1] = Cinv
        else: Cinv = n.linalg.inv(self.get_C(bl1,bl2))
        return Cinv
    def get_z(self, bl):
        i = self.bls.index(bl)
        z = self.zs.get(bl, n.dot(self.get_Cinv(bl), self.X[i*self.nprms:(i+1)*self.nprms]))
        self.zs[bl] = z
        return z

def gen_Q(knum, n_k):
    Q = n.zeros((n_k,n_k))
    Q[knum,knum] = 1
    return Q


n_k = chans.size
bls = T.keys()
Ts = n.concatenate([T[bl] for bl in bls], axis=-1).T
print Ts.shape

def gen_Q_old(knum, n_k):
    Q = n.zeros((Ts.shape[0],Ts.shape[0]))
    for i in n.arange(0,Q.shape[0],n_k):
        for j in n.arange(0,Q.shape[1],n_k):
            Q[i+knum,j+knum] = 1
    return Q

#p.subplot(121)
#capo.arp.waterfall(Ts, mode='log', drng=2)
#p.subplot(122)
#capo.arp.waterfall(Ts, mode='phs')
#p.show()

C = CoV(Ts, bls)
Cold = cov(Ts) # XXX old

dsum,dwgt = 0.,0.
p.subplot(121); capo.arp.waterfall(Cold, mode='log', drng=3)
p.subplot(122)
for cnt,bl0 in enumerate(bls):
  for bl1 in bls[cnt:]:
    d = C.get_C(bl0,bl1).diagonal()
    dsum += d; dwgt += 1
    p.plot(d)
dspec = dsum / dwgt
p.plot(dspec, 'k:')
p.show()

if False: # Substitute an idealized covariance matrix
    if True:
        Cold_ = sum([gen_Q_old(k, n_k) for k in xrange(n_k)])
        Cold_ *= Cold
    else:
        Cold_ = sum([dspec[k]*gen_Q_old(k, n_k) for k in xrange(n_k)])
    p.subplot(121); capo.arp.waterfall(Cold, mode='log', drng=3)
    p.subplot(122); capo.arp.waterfall(Cold_, mode='log', drng=3)
    p.show()
    Cold = Cold_

if False: # Use pseudoinverse for singular matrices
    # The problem is that Cold is singular, so inverting it generates problems
    if True: # Use Moore-Penrose pseudoinverse
        #_Cold = n.linalg.pinv(Cold) # XXX old
        U = n.zeros((Ts.shape[0], n_k))
        for i in xrange(n_k): U[i,i], U[i+n_k,i] = 1,-1
        P = n.eye(Ts.shape[0]) - n.dot(U, n.dot(n.linalg.inv(n.dot(U.T.conj(),U)), U.T.conj()))
        Cold = n.dot(P,n.dot(Cold,P))
        _Cold = n.linalg.pinv(Cold) # XXX old
    else:
        U = n.zeros((Ts.shape[0], n_k))
        for i in xrange(n_k): U[i,i], U[i+n_k,i] = 1,-1
        P = n.eye(Ts.shape[0]) - n.dot(U, n.dot(n.linalg.inv(n.dot(U.T.conj(),U)), U.T.conj()))
        p.subplot(121); capo.arp.waterfall(Cold, mode='log', drng=3)
        p.subplot(122); capo.arp.waterfall(n.dot(P,n.dot(Cold,P)), mode='log', drng=3)
        p.show()
        _Cold = n.dot(P, n.dot(n.linalg.inv(Cold + 100.*n.dot(U,U.T.conj())), P))
    p.subplot(121); capo.arp.waterfall(Cold, mode='log', drng=3)
    p.subplot(122); capo.arp.waterfall(_Cold, mode='log', drng=3)
    p.show()
else:
    _Cold = n.linalg.inv(Cold) # XXX old

if False:
    p.subplot(131); capo.arp.waterfall(Ts, mode='log', drng=2)
    p.subplot(132); capo.arp.waterfall(n.dot(_Cold,Ts), mode='log', drng=2)
    p.subplot(133); capo.arp.waterfall(cov(n.dot(_Cold,Ts)), mode='log', drng=2)
    p.show()

#Tbl0 = n.array(T[bls[0]]).T; Tbl0 -= Tbl0.mean(axis=1)[(slice(None),n.newaxis)]
#Tbl1 = n.array(T[bls[1]]).T; Tbl1 -= Tbl1.mean(axis=1)[(slice(None),n.newaxis)]
#p.subplot(131); capo.arp.waterfall(C.get_C(bls[0]), mode='log', drng=2); p.colorbar(shrink=.5)
#p.subplot(132); capo.arp.waterfall(cov(Tbl0), mode='log', drng=2); p.colorbar(shrink=.5)
#p.subplot(133); capo.arp.waterfall(n.dot(Tbl0, Tbl0.T.conj())/(Tbl0.shape[1]-1), mode='log', drng=2); p.colorbar(shrink=.5)
#p.show()
#p.subplot(221); capo.arp.waterfall(C.C, mode='log', drng=2)
#p.subplot(222); capo.arp.waterfall(C.get_C(bls[0]), mode='log', drng=2)
#p.subplot(223); capo.arp.waterfall(C.get_C(bls[1]), mode='log', drng=2)
#p.subplot(224); capo.arp.waterfall(C.get_C(bls[0],bls[1]), mode='log', drng=2)
#p.show()
#p.subplot(131); capo.arp.waterfall(n.dot(Tbl0, Tbl1.T.conj()), mode='log', drng=2)
#p.subplot(132); capo.arp.waterfall(n.dot(Tbl1, Tbl0.T.conj()), mode='log', drng=2)
#p.subplot(133); capo.arp.waterfall(C.get_C(bls[0],bls[1]), mode='log', drng=2)
#p.show()
Qs,Q_Z = {}, {}
Qs_old, Q_C = {}, {}
for k in range(n_k):
    Qs[k] = gen_Q(k,n_k)
    Q_Z[k] = {}
    for bl in bls:
        Q_Z[k][bl] = n.dot(Qs[k], C.get_Cinv(bl)) # C^-1 = Z
    # XXX old
    Qs_old[k] = gen_Q_old(k,n_k)
    #p.imshow(Qs_old[k], interpolation='nearest'); p.show()
    Q_C[k] = n.dot(Qs_old[k], _Cold)


for bl1 in T:
    z1 = C.get_z(bl1)
    _C1 = Z11 = C.get_Cinv(bl1) # Z11 = Cinv_1
    #Z1 = n.cov(z1) # XXX 
    #p.subplot(121); capo.arp.waterfall(Z11, mode='log', drng=2)
    #p.subplot(122); capo.arp.waterfall(Z1, mode='log', drng=2)
    #p.show()
    for bl2 in T:
        if bl2 == bl1: continue # skip auto-products
        z2 = C.get_z(bl2)
        _C2 = Z22 = C.get_Cinv(bl2) # Z22 = Cinv_2
        #Z2 = n.cov(z2) # XXX 
        #p.subplot(121); capo.arp.waterfall(Z22, mode='log', drng=2)
        #p.subplot(122); capo.arp.waterfall(Z2, mode='log', drng=2)
        #p.show()
        Z21 = n.dot(_C2, n.dot(C.get_C(bl2,bl1), _C1)) # Z21 = C^-1_2 C_21 C^-1_1
        #Z21alt = n.dot(z2-z2.mean(axis=1)[(slice(None),n.newaxis)],(z1-z1.mean(axis=1)[(slice(None),n.newaxis)]).T.conj())
        #p.subplot(121); capo.arp.waterfall(Z21, mode='log', drng=2)
        #p.subplot(122); capo.arp.waterfall(Z21alt, mode='log', drng=2)
        #p.show()
        #Z21 = Z21alt
        '''
        p.subplot(131); capo.arp.waterfall(Z11, mode='log', drng=2)
        p.subplot(132); capo.arp.waterfall(Z21, mode='log', drng=2)
        p.subplot(133); capo.arp.waterfall(Z22, mode='log', drng=2)
        p.show()
        '''
        print 'Making Fisher/Gisher Matrix'
        F,G = n.zeros((n_k,n_k), dtype=n.complex), n.zeros((n_k,n_k), dtype=n.complex)
        Fold = n.zeros((n_k,n_k), dtype=n.complex) # XXX old
        Q_Z21 = {}
        for k in Qs: Q_Z21[k] = n.dot(Qs[k], Z21)

        for k1 in range(n_k):
          for k2 in range(n_k):
            F[k1,k2] = 0.5 * quick_trace_dot(Q_Z[k1][bl1],Q_Z[k2][bl2])
            G[k1,k2] = 0.25 * quick_trace_dot(Q_Z21[k1],Q_Z21[k2]) \
                    + 0.25 * quick_trace_dot(Q_Z[k1][bl2],Q_Z[k2][bl1])
            Fold[k1,k2] = 0.5 * quick_trace_dot(Q_C[k1],Q_C[k2]) # XXX old
        _G = n.linalg.inv(G)
        # XXX begin old
        _Fold = n.linalg.inv(Fold)
        p.subplot(121); capo.arp.waterfall(Fold, mode='log', drng=2)
        p.subplot(122); capo.arp.waterfall(_Fold, mode='log', drng=2)
        p.show()
        w,v = n.linalg.eig(_Fold)
        Mold = n.dot(v, n.dot(n.diagflat(n.sqrt(w)), v.T))
        #Mold = _Fold
        #Mold = n.eye(Mold.shape[0]) # override with diagonal, in case inverting F is problematic
        Wold = n.dot(Mold,Fold)
        norm = n.sum(Wold, axis=-1); norm.shape = (norm.size,1)
        Mold /= norm
        Wold = n.dot(Mold,Fold)
        print Wold.diagonal()
        p.subplot(131); capo.arp.waterfall(Wold, mode='log', drng=2); p.colorbar(shrink=.5)
        p.subplot(132); capo.arp.waterfall(Mold, mode='log', drng=2); p.colorbar(shrink=.5)
        p.subplot(133); capo.arp.waterfall(Fold, mode='log', drng=2); p.colorbar(shrink=.5)
        p.show()
        # XXX end old

        print 'Making Normalization/Windowing'
        # Set M = G^-1/2
        w,v = n.linalg.eig(_G)
        M = n.dot(v, n.dot(n.diagflat(n.sqrt(w)), v.T))
        # Normalize M s.t. rows of W sum to 1
        W = n.dot(M,F)
        norm = n.sum(W, axis=-1); norm.shape = (norm.size,1)
        M /= norm
        W = n.dot(M,F)
        
        '''
        p.subplot(321); capo.arp.waterfall(C, mode='log', drng=2)
        p.subplot(322); capo.arp.waterfall(_C, mode='log', drng=2)
        p.subplot(323); capo.arp.waterfall(F, mode='log', drng=2)
        p.subplot(324); capo.arp.waterfall(_F, mode='log', drng=2)
        p.subplot(325); capo.arp.waterfall(M, mode='log', drng=2)
        p.subplot(326); capo.arp.waterfall(n.dot(M,n.dot(F,M.T.conj())), mode='log', drng=2)
        p.show()
        '''

        print 'Minding ps and qs'
        qs = []
        qold = [] # XXX old
        for k in range(n_k):
            qs.append(0.5 * quick_diag_dot(z1.T.conj(), n.dot(Qs[k], z2)))
            _CQ_C = n.dot(_Cold, n.dot(Qs_old[k], _Cold)) # XXX old
            qold.append(0.5 * quick_diag_dot(Ts.T.conj(), n.dot(_CQ_C, Ts))) # XXX old
        qold = n.array(qold) # XXX old
        pold = scalar*n.dot(Mold,qold) # XXX old
        p.subplot(121); capo.arp.waterfall(qold, mode='log', drng=2)
        p.subplot(122); capo.arp.waterfall(pold, mode='log', drng=2)
        p.show()
        qs = n.array(qs)
        ps = scalar*n.dot(M,qs)
        print ps.shape
        #p.subplot(131); capo.arp.waterfall(ps, mode='log', drng=3)
        p.subplot(121); capo.arp.waterfall(pold, mode='log', drng=3)
        #p.subplot(133); p.plot(n.average(ps, axis=1)); p.plot(n.average(pold, axis=1))
        p.subplot(122)
        p.plot(n.average(pold, axis=1))
        p.plot(dspec)
        p.show()
        # XXX still have to get Sigma

pk_a_priori = Ts[:,0:40] * n.conj(Ts[:,40:80])
print n.average(pk_a_priori, axis=0).real

print 'Making Covariance Matrix'
C = n.cov(Ts, rowvar=False)
print 'C.shape =', C.shape
_C = n.linalg.inv(C)


#_C = n.eye(_C.shape[0])

if False: # Create a fake signal with cov of original signal
    print 'Creating a fake signal with CoV of original signal'
    print Ts.shape
    #Ts_ = n.random.normal(size=Ts.shape) * n.exp(2j*n.pi*n.random.uniform(size=Ts.shape))
    Ts_ = n.exp(2j*n.pi*n.random.uniform(size=Ts.shape))
    #for i in n.arange(0,Ts.shape[1],n_k):
    #    for t in xrange(Ts.shape[0]):
    #        Ts_[t,i+26] += 1. * n.exp(2j*n.pi*float(t)/10)
    L = n.linalg.cholesky(C)
    Ts_ = n.dot(L, Ts_.T).T
    print Ts_.dtype
    #Ts_ = Ts * n.exp(2j*n.pi*n.random.uniform(size=(Ts.shape[0],1)))
    p.subplot(121); capo.arp.waterfall(Ts, mode='log', drng=1); p.colorbar(shrink=.5)
    p.subplot(122); capo.arp.waterfall(Ts_, mode='log', drng=1); p.colorbar(shrink=.5)
    p.show()
    #Ts = Ts_
    C = n.cov(Ts, rowvar=False)
    _C = n.linalg.inv(C)

print 'Making Fisher Matrix'
Qs = {}
for k in range(n_k): Qs[k] = gen_Q(k,n_k)
F = n.zeros((40,40), dtype=n.complex)
Q_C = {}
for k in Qs: Q_C[k] = n.dot(Qs[k], _C)

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
W = n.dot(M,F)

#capo.arp.waterfall(W, mode='log', drng=2); p.show()

print 'Minding ps and qs'
qs, ps = [], []
for k in range(n_k):
    print k
    _CQ_C = n.dot(_C, n.dot(Qs[k], _C))
    #qs.append(0.5 * n.dot(Ts_c, n.dot(_CQ_C, Ts_t)))
    qs.append(0.5 * quick_diag_dot(Ts.conj(), n.dot(_CQ_C, Ts.T)))
qs = n.array(qs)
print qs.shape, M.shape
ps = n.dot(M,qs)
pk = scalar * ps.T
print pk.shape

etas = n.fft.fftshift(capo.pspec.f2eta(afreqs))
kpl = etas * capo.pspec.dk_deta(z)

print 'Generating bootstrap errors'
NBOOT2 = 1
BOOTLEN2 = 100
avg_1ds, std_1ds = [], []
for h in xrange(NBOOT2):
  boot = []
  for i in xrange(BOOTLEN2):
    dsum, dwgt = 0., 0.
    if i % 10 == 0: print '%d/%d-%d/%d' % (h,NBOOT2,i,BOOTLEN2)
    N = pk.shape[0]
    if True: # do average (as opposed to median)
        for j in n.random.randint(0,N,N):
            #w = 1./std_2d[j]**2
            w = 1.
            dsum += pk[j] * w
            dwgt += w
        boot.append(dsum/dwgt)
    else: # do median instead of average
        samples = n.array([pk[j] for j in n.random.randint(0,N,N)])
        boot.append(n.median(samples, axis=0))
  boot = n.array(boot)
  '''
  p.clf(); p.plot(pk[:,26].real, pk[:,26].imag, '.')
  p.plot(boot[:,26].real, boot[:,26].imag, 'x'); p.show()
  '''
  print boot.shape
  avg_1d = n.average(boot, axis=0)
  #med_1d = n.median(boot, axis=0)
  std_1d = n.std(boot, axis=0)
  avg_1ds.append(avg_1d)
  std_1ds.append(std_1d)
avg_1ds = n.array(avg_1ds); avg_1d = n.average(avg_1ds, axis=0)
std_1ds = n.array(std_1ds); std_1d = n.average(std_1ds, axis=0)
print avg_1ds[:,26], avg_1d[26]
print std_1ds[:,26], std_1d[26]

print 'Writing pspec.npz'
'''
avg_1d = n.average(pk, axis=0)
std_1d = n.std(pk, axis=0) / n.sqrt(pk.shape[0])
for _kpl,_pk,_err in zip(kpl,avg_1d,std_1d):
    print _kpl, _pk.real, _err
'''
n.savez('pspec.npz', kpl=kpl, pk=avg_1d, err=std_1d)

'''
#p.subplot(121); capo.arp.waterfall(pk, mode='log', drng=4); p.colorbar(shrink=.5)
p.subplot(121); capo.arp.waterfall(pk, mode='real', drng=10); p.colorbar(shrink=.5)
p.subplot(122); capo.arp.waterfall(pk, mode='phs')
p.show()
#norm = n.sqrt(n.sum(n.abs(_cov)**2, axis=-1))
#norm.shape = (norm.size,1)
#_cov /= norm # XXX need to get normalization correct

import sys; sys.exit(0)
# Go through data again, this time applying inverse cov matrix
for filename in args:
    outfile = filename + '.pspec'
    print filename,'->',outfile
    if os.path.exists(outfile):
        print '    File exists.  Skipping...'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi, override={'nchan':afreqs.size, 'sfreq':afreqs[0]})
    uvo._wrhd('history', uvi['history'] + 'PSPEC: ' + ' '.join(sys.argv) + '\n')
    #uvo.add_var('k3pk_fq', 'r')
    #uvo.add_var('k3pk_wgt', 'r')

    a.scripting.uv_selector(uvi, opts.ant, opts.pol)

    _Tlist,_Wlist,curtime = {},{},None
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        if t != curtime:
            #print t
            uvo.copyvr(uvi)
            for sep,bls in sep2bl.items():
                for cnt,bl0 in enumerate(bls):
                    if not _Tlist.has_key(bl0): continue
                    for bl1 in bls[cnt:]: # this includes "auto-pspecs"
                        if not _Tlist.has_key(bl1): continue
                        pk = scalar * _Tlist[bl0] * n.conj(_Tlist[bl1])
                        uvo.write((crd,curtime,(bl_index(bl0),bl_index(bl1))), pk, n.zeros(pk.shape, dtype=n.int))
            # Clear the current pspec data and start a new integration
            _Tlist,_Wlist = {},{}
            curtime = t

        bl = a.miriad.ij2bl(i,j)
        sep = bl2sep[bl]
        if sep < 0:
            #print 'Conj:', a.miriad.bl2ij(bl)
            d,sep = n.conj(d),-sep

        d,f = d.take(chans), f.take(chans)
        w = n.logical_not(f).astype(n.float)
        Trms = d * C.pspec.jy2T(afreqs)
        if PFB:
            _Trms = C.pfb.pfb(Trms, window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
            _Wrms = C.pfb.pfb(w   , window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
        else:
            window = a.dsp.gen_window(Trms.size, WINDOW)
            _Trms = n.fft.ifft(window * Trms)
            _Wrms = n.fft.ifft(window * w)
        gain = n.abs(_Wrms[0])
        #print 'Gain:', gain
        if True: # we think clean messes up noise statistics.  Do we need it? (Yes)
            if gain > 0:
                #_Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9, maxiter=100, stop_if_div=False, verbose=False)
                #_Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9)
                _Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-2)
                #print info['term']
                #_Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9)
                _Trms = _Tcln + info['res'] / gain
        else: _Trms /= gain
        _Trms = n.fft.fftshift(_Trms)
        if True:
            _Trms.shape = (_Trms.size,1)
            print n.median(n.abs(_Trms)**2),
            _Trms = n.dot(covs[bl],_Trms).flatten()
            print n.median(n.abs(_Trms)**2)
        _Wrms = n.fft.fftshift(_Wrms)

        _Tlist[bl] = _Trms
        _Wlist[bl] = _Wrms

    # Gotta do this one last time to catch the last integration.
    for sep,bls in sep2bl.items():
        for cnt,bl0 in enumerate(bls):
            if not _Tlist.has_key(bl0): continue
            for bl1 in bls[cnt:]: # this includes "auto-pspecs"
                if not _Tlist.has_key(bl1): continue
                pk = scalar * _Tlist[bl0] * n.conj(_Tlist[bl1])
                uvo.write((crd,curtime,(bl_index(bl0),bl_index(bl1))), pk, n.zeros(pk.shape, dtype=n.int))
    del(uvi); del(uvo)
'''
