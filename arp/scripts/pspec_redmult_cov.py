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

T, N, W = {}, {}, {}
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
            #g = 1.
            #Trms = g*Trms + (1-g)*Trms_
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
            if False:
                _Tcln, info = a.deconv.clean(_Trms, _Wrms, tol=1e-2)
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
                _Trms = n.dot(_C, _Trms).squeeze()
                _Nrms = n.dot(_C, _Nrms).squeeze()
        _Trms = n.fft.fftshift(_Trms)
        _Nrms = n.fft.fftshift(_Nrms)
        _Wrms = n.fft.fftshift(_Wrms)
        if False: # swap in a simulated signal post delay transform
            _Trms = n.random.normal(size=_Trms.size) * n.exp(2j*n.pi*n.random.uniform(size=_Trms.size))
            #_Trms[26] += 2 * n.exp(100j*t)
        T[bl] = T.get(bl, []) + [_Trms]
        N[bl] = N.get(bl, []) + [_Nrms]
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

n_k = chans.size
bls = T.keys()
Ts = n.concatenate([T[bl] for bl in bls], axis=-1).T
Ns = n.concatenate([N[bl] for bl in bls], axis=-1).T
print Ts.shape
if False:
    p.subplot(221); capo.arp.waterfall(Ts, mode='log', mx=1, drng=2); p.colorbar(shrink=.5)
    p.subplot(222); capo.arp.waterfall(Ns, mode='log', mx=1, drng=2); p.colorbar(shrink=.5)
    p.subplot(223); capo.arp.waterfall(n.fft.ifft(Ts,axis=1), mode='log', mx=0,drng=2); p.colorbar(shrink=.5)
    p.subplot(224); capo.arp.waterfall(n.fft.ifft(Ns,axis=1), mode='log', mx=0,drng=2); p.colorbar(shrink=.5)
    p.show()

print ' '.join(['%d_%d' % a.miriad.bl2ij(bl) for bl in bls])
C = CoV(Ts, bls)
Cn = CoV(Ns, bls)
capo.arp.waterfall(C.C, mode='log', drng=2); p.show()

X = Ts.copy()
N1 = Ns.copy()
N2 = Ns.copy()
PLT1,PLT2 = 3,3
for cnt in xrange(PLT1*PLT2-1):
    p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(cov(X), mode='log', mx=0, drng=2)
    #p.subplot(PLT1,PLT2,i+1); capo.arp.waterfall(cov(N2), mode='log', mx=0, drng=2)
    SZ = X.shape[0]
    ints = n.arange(X.shape[1])
    #DEC = 8, 10, 20, 5
    #DEC = 8
    DEC = 1
    seta = n.array([xi for xi in ints if (xi/DEC)%2 == 0])
    setb = n.array([xi for xi in ints if (xi/DEC)%2 == 1])
    xa,xb = X[:,seta],X[:,setb]
    n1a,n1b = N1[:,seta],N1[:,setb]
    n2a,n2b = N2[:,seta],N2[:,setb]
    Ca,Cb = cov(xa), cov(xb)
    Cna,Cnb = cov(n2a), cov(n2b)

    # Normalize covariance matrices
    for c in [Ca,Cb,Cna,Cnb]:
        d = n.diag(c); d.shape = (1,SZ)
        c /= n.sqrt(d) * 2

    #gain1,gain2 = 1,1
    #gain1 = .1
    gain1 = 1. / len(bls)
    #gain1 = .5 / len(bls)
    #gain1,gain2 = .3, .1
    _Ca = -gain1*Ca
    _Cb = -gain1*Cb
    _Cna = -gain1*Cna
    _Cnb = -gain1*Cnb
    #_C1,_C2 = n.zeros_like(C1), n.zeros_like(C2) # don't diagonalize auto-products
    #_C12,_C21 = n.zeros_like(C12), n.zeros_like(C21) # don't diagonalize cross-products
    ind = n.arange(SZ)
    # zero out redundant off-diagonals
    for b in xrange(len(bls)):
        indb = ind[:-b*n_k]
        _Ca[indb,indb+b*n_k] = 0
        _Ca[indb+b*n_k,indb] = 0
        _Cna[indb,indb+b*n_k] = 0
        _Cna[indb+b*n_k,indb] = 0
        _Cb[indb,indb+b*n_k] = 0
        _Cb[indb+b*n_k,indb] = 0
        _Cnb[indb,indb+b*n_k] = 0
        _Cnb[indb+b*n_k,indb] = 0
        if False: # only subtract modes near the horizon
            for c in xrange(len(bls)):
                h1,h2 = 12,29
                _Ca[n_k*b:n_k*b+h1,n_k*c:n_k*c+h1] = 0
                _Ca[n_k*b:n_k*b+h1,n_k*c+h2:n_k*(c+1)] = 0
                _Ca[n_k*b+h2:n_k*(b+1),n_k*c:n_k*c+h1] = 0
                _Ca[n_k*b+h2:n_k*(b+1),n_k*c+h2:n_k*(c+1)] = 0
                _Cb[n_k*b:n_k*b+h1,n_k*c:n_k*c+h1] = 0
                _Cb[n_k*b:n_k*b+h1,n_k*c+h2:n_k*(c+1)] = 0
                _Cb[n_k*b+h2:n_k*(b+1),n_k*c:n_k*c+h1] = 0
                _Cb[n_k*b+h2:n_k*(b+1),n_k*c+h2:n_k*(c+1)] = 0
                _Cna[n_k*b:n_k*b+h1,n_k*c:n_k*c+h1] = 0
                _Cna[n_k*b:n_k*b+h1,n_k*c+h2:n_k*(c+1)] = 0
                _Cna[n_k*b+h2:n_k*(b+1),n_k*c:n_k*c+h1] = 0
                _Cna[n_k*b+h2:n_k*(b+1),n_k*c+h2:n_k*(c+1)] = 0
                _Cnb[n_k*b:n_k*b+h1,n_k*c:n_k*c+h1] = 0
                _Cnb[n_k*b:n_k*b+h1,n_k*c+h2:n_k*(c+1)] = 0
                _Cnb[n_k*b+h2:n_k*(b+1),n_k*c:n_k*c+h1] = 0
                _Cnb[n_k*b+h2:n_k*(b+1),n_k*c+h2:n_k*(c+1)] = 0
    if True: # subtract off covariances common to all baselines pairs
        _Ca[ind,ind] = _Cb[ind,ind] = 0
        _Cna[ind,ind] = _Cnb[ind,ind] = 0
        for _C in [_Ca,_Cb,_Cna,_Cnb]:
            _C.shape = (len(bls),n_k,len(bls),n_k)
            sub_C = n.zeros_like(_C)
            for i in xrange(len(bls)):
                for j in xrange(len(bls)):
                    not_ij = n.array([bl for bl in xrange(len(bls)) if not bl in [i,j]])
                    _C_avg = n.mean(n.mean(_C.take(not_ij,axis=0).take(not_ij,axis=2), axis=2), axis=0)
                    sub_C[i,:,j,:] = _C_avg
            _C.shape = (len(bls)*n_k,len(bls)*n_k)
            sub_C.shape = (len(bls)*n_k,len(bls)*n_k)
            #p.clf()
            #p.subplot(131); capo.arp.waterfall(sub_C, mode='log', drng=2)
            #p.subplot(132); capo.arp.waterfall(_C, mode='log', drng=2)
            #p.subplot(133); capo.arp.waterfall(_C - sub_C, mode='log', drng=2)
            #p.show()
            _C -= sub_C
            
    _Ca[ind,ind] = _Cb[ind,ind] = 1
    _Cna[ind,ind] = _Cnb[ind,ind] = 1
    #p.clf(); capo.arp.waterfall(_Ca, mode='log', drng=2); p.show()
    
    # Use covariance of first half of data to correct second half & vice versa
    xa,xb = n.dot(_Cb,xa), n.dot(_Ca,xb)
    n1a,n1b = n.dot(_Cb,n1a), n.dot(_Cb,n1b)
    n2a,n2b = n.dot(_Cnb,n2a), n.dot(_Cnb,n2b)
    X[:,seta],X[:,setb] = xa,xb
    N1[:,seta],N1[:,setb] = n1a,n1b
    N2[:,seta],N2[:,setb] = n2a,n2b
## Normalize to maintain amplitude assuming all k modes are independent & same amplitude
#norm1 = n.sqrt(n.sum(n.abs(_C1tot)**2, axis=1)); norm1.shape = (norm1.size,1)
#norm2 = n.sqrt(n.sum(n.abs(_C2tot)**2, axis=1)); norm2.shape = (norm2.size,1)
#x1 /= norm1; x2 /= norm2
#n1 /= norm1; n2 /= norm2
#p.subplot(121); capo.arp.waterfall(_C1tot/norm1, mode='lin', drng=1, mx=1); p.colorbar(shrink=.5)
#p.subplot(122); capo.arp.waterfall(_C2tot/norm2, mode='lin', drng=1, mx=1); p.colorbar(shrink=.5)
#p.show()
p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(X), mode='log', mx=0, drng=2)
#p.subplot(PLT1,PLT2,i+2); capo.arp.waterfall(cov(N2), mode='log', mx=0, drng=2)
p.show()

pspecs = []
C_ = CoV(X, bls)
Cn1_ = CoV(N1, bls)
Cn2_ = CoV(N2, bls)
for cnt,bli in enumerate(bls):
    for blj in bls[cnt:]:
        if bli == blj: continue
        print a.miriad.bl2ij(bli), a.miriad.bl2ij(blj)
        xi,xj = C.get_x(bli), C.get_x(blj)
        xi_,xj_ = C_.get_x(bli), C_.get_x(blj)
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
#std_1d = n.ones_like(avg_1d)

print 'Writing pspec.npz'
n.savez('pspec.npz', kpl=kpl, pk=avg_1d, err=std_1d)
p.show()

