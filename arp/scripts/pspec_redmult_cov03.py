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
                #_Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9, maxiter=100, stop_if_div=False, verbose=False)
                #_Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9)
                _Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-2)
                #print info['term']
                #_Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-9)
                _Trms = _Tcln + info['res'] / gain
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
        return self.C[i*self.nprms:(i+1)*self.nprms, j*self.nprms:(j+1)*self.nprms]
    def get_Cinv(self, bl1, bl2=None):
        if bl2 is None: bl2 = bl1
        return n.linalg.inv(self.get_C(bl1,bl2))
    def get_x(self, bl):
        i = self.bls.index(bl)
        return self.X[i*self.nprms:(i+1)*self.nprms]
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
    

C = CoV(Ts, bls)
p.subplot(131); capo.arp.waterfall(C.C, mode='log', drng=2); p.colorbar(shrink=.5)
p.subplot(132); capo.arp.waterfall(C.get_Ck(17), mode='real'); p.colorbar(shrink=.5)
p.subplot(133); capo.arp.waterfall(C.get_Ck(17), mode='imag'); p.colorbar(shrink=.5)

'''
            d = n.diag(C21)
            R21 = C21/n.sqrt(n.multiply.outer(d,d))
p.subplot(234); capo.arp.waterfall(C.C/n.sqrt(n.mult, mode='log', drng=2); p.colorbar(shrink=.5)
p.subplot(235); capo.arp.waterfall(C.get_Ck(17), mode='real'); p.colorbar(shrink=.5)
p.subplot(236); capo.arp.waterfall(C.get_Ck(23), mode='real'); p.colorbar(shrink=.5)
'''
p.show()
#Cold = cov(Ts) # XXX old

#dsum,dwgt = 0.,0.
#p.subplot(121); capo.arp.waterfall(Cold, mode='log', drng=3)
#p.subplot(122)
#for cnt,bl0 in enumerate(bls):
#  for bl1 in bls[cnt:]:
#    d = C.get_C(bl0,bl1).diagonal()
#    dsum += d; dwgt += 1
#    p.plot(d)
#dspec = dsum / dwgt
#p.plot(dspec, 'k:')
#p.show()

Q = {}
for k in range(n_k): Q[k] = gen_Q(k,n_k)

for bl1 in T:
    for bl2 in T:
        #if bl2 == bl1: continue # skip auto-products
        print a.miriad.bl2ij(bl1), a.miriad.bl2ij(bl2)

        x1 = C.get_x(bl1)
        x2 = C.get_x(bl2)
        C21 = C.get_C(bl2,bl1)
        C12 = C.get_C(bl1,bl2)
        SZ = C21.shape[0]
        if False: # Artificially diagonalize covariance matrix
            mask = n.eye(SZ)
            #mask[1:] += n.eye(SZ)[:-1]
            #mask[:-1] += n.eye(SZ)[1:]
            #mask[20] = 1; mask[:,20] = 1
            #mask[23] = 1; mask[:,23] = 1
            #C21 *= mask
            #C12 *= mask
            C21 = mask
            C12 = mask
        if False: # Add some extra power to the diagonal
            amp = n.abs(C12.diagonal())
            thresh = amp.max()/2
            mask = n.where(amp >= thresh, 0, thresh)
            #C12 += 3*n.diagflat(n.random.normal(size=SZ))
            #C21 += 3*n.diagflat(n.random.normal(size=SZ))
            #C12 += n.diagflat(mask*n.exp(2j*n.pi*n.random.uniform(size=SZ)))
            #C21 += n.diagflat(mask*n.exp(2j*n.pi*n.random.uniform(size=SZ)))
            C12 += n.diagflat(mask)
            C21 += n.diagflat(mask)
        #C21 *= n.abs(C21)
        #C12 *= n.abs(C12)
        p.subplot(121); capo.arp.waterfall(C21, mode='log', drng=6); p.colorbar(shrink=.5)
        p.subplot(122); capo.arp.waterfall(C12, mode='log', drng=6); p.colorbar(shrink=.5)
        p.show()

        if True:
            _C21 = n.linalg.inv(C21)
            _C12 = n.linalg.inv(C12)
        else: # Alternate more numerically robust way to compute inverse
            d = n.diag(C21)
            R21 = C21/n.sqrt(n.multiply.outer(d,d))
            D21 = n.diagflat(n.sqrt(n.diag(C21)))
            d = n.diag(C12)
            R12 = C12/n.sqrt(n.multiply.outer(d,d))
            D12 = n.diagflat(n.sqrt(n.diag(C12)))
            #w,v = n.linalg.eig(R21)
            _D21 = n.linalg.inv(D21)
            _C21 = n.dot(_D21, n.dot(n.linalg.inv(R21), _D21))
            _D12 = n.linalg.inv(D12)
            _C12 = n.dot(_D12, n.dot(n.linalg.inv(R12), _D12))

        p.subplot(121); capo.arp.waterfall(_C21, mode='log', drng=6); p.colorbar(shrink=.5)
        p.subplot(122); capo.arp.waterfall(_C12, mode='log', drng=6); p.colorbar(shrink=.5)
        p.show()

        print n.all(_C21 == _C12)
        z1 = n.dot(_C21, x1)
        z2 = n.dot(_C12, x2)

        p.subplot(231); capo.arp.waterfall(x1, mode='log', drng=2)
        p.subplot(232); capo.arp.waterfall(_C21, mode='log', drng=2)
        p.subplot(233); capo.arp.waterfall(z1, mode='log', drng=2)
        p.subplot(234); capo.arp.waterfall(x2, mode='log', drng=2)
        p.subplot(235); capo.arp.waterfall(_C12, mode='log', drng=2)
        p.subplot(236); capo.arp.waterfall(z2, mode='log', drng=2)
        p.show()

        print 'Making Fisher/Gisher Matrix'
        F,G = n.zeros((n_k,n_k), dtype=n.complex), n.zeros((n_k,n_k), dtype=n.complex)

        q = []
        Q_C = {}
        for k in range(n_k):
            q.append(0.5 * quick_diag_dot(z1.T.conj(), n.dot(Q[k], z2)))
            Q_C[k] = n.dot(Q[k], _C12)
        q = n.array(q)
        for k1 in range(n_k):
            for k2 in range(n_k):
                F[k1,k2] = 0.5 * quick_trace_dot(Q_C[k1],Q_C[k2])
        G = cov(q)

        print 'Making Normalization/Windowing'
        #_G = n.linalg.inv(G)
        _G = n.linalg.inv(F)
        if False:# Set M = G^-1/2
            w,v = n.linalg.eig(_G)
            M = n.dot(v, n.dot(n.diagflat(n.sqrt(w)), v.T))
        else:
            M = _G
        # Normalize M s.t. rows of W sum to 1
        W = n.dot(M,F)
        norm = n.sum(W, axis=-1); norm.shape = (norm.size,1)
        M /= norm
        print n.diag(M)
        capo.arp.waterfall(_G, mode='log', drng=4); p.show()
        W = n.dot(M,F)
        print 'Wdiag', W.diagonal()
        p.subplot(111); capo.arp.waterfall(W, mode='log', drng=2); p.show()
        pk = scalar * n.dot(M,q)
        p.subplot(131); capo.arp.waterfall(q, mode='real'); p.colorbar(shrink=.5)
        p.subplot(132); capo.arp.waterfall(M, mode='log', drng=2)
        p.subplot(133); capo.arp.waterfall(pk, mode='real'); p.colorbar(shrink=.5)
        p.show()
        
        pk_avg = n.average(pk, axis=1)
        p.plot(pk_avg.real)
        p.plot(n.average(x1*x2.conj(), axis=1).real)
        #pk_true = n.zeros_like(pk_avg); pk_true[26] += 4
        pk_true = n.ones_like(pk_avg); pk_true[26] += 4
        p.plot(n.dot(W, pk_true))
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
