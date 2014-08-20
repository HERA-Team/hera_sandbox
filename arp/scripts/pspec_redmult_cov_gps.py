#! /usr/bin/env python
import aipy as a, numpy as n
import capo
import optparse, sys, os, random

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True)
o.add_option('-t', '--taps', type='int', default=1,
    help='Taps to use in the PFB.  Default 1, which instead uses windowed FFT')
o.add_option('-b', '--boot', type='int', default=20,
    help='Number of bootstraps.  Default is 20')
o.add_option('--plot', action='store_true',
    help='Generate plots')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
opts,args = o.parse_args(sys.argv[1:])

PLOT = opts.plot
if PLOT: import pylab as p

NBOOT = opts.boot
NTAPS = opts.taps
if NTAPS > 1: PFB = True
else: PFB = False
WINDOW = opts.window

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
    #B = sdf * afreqs.size / capo.pfb.NOISE_EQUIV_BW[WINDOW] # this is wrong if we aren't inverting
    # the window post delay transform (or at least dividing out by the gain of the window)
    # For windowed data, the FFT divides out by the full bandwidth, B, which is
    # then squared.  Proper normalization is to multiply by B**2 / (B / NoiseEqBand) = B * NoiseEqBand
    B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW]
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
    print 'Reading', filename
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
            NDAY = 92
            NBL = 1
            NPOL = 2
            T_INT = 43. # for just compressed data
            #T_INT = 351. # for fringe-rate filtered data.  don't use if fringe-rate filtering is active below
            Trms_ = n.random.normal(size=Trms.size) * n.exp(2j*n.pi*n.random.uniform(size=Trms.size))
            Trms_ *= TSYS / n.sqrt(B * T_INT * NDAY * NBL * NPOL)
            #Trms_ *= n.sqrt(n.sqrt(351./43)) # penalize for oversampling fr-filtered data
            Trms_ *= 1.14 # adjust to suboptimal flux calibration
            #Trms_ += 10*eor_mdl[t] # add in a common signal to all baselines
            Trms_ = eor_mdl[t] # override noise with common signal
            Nrms  = Trms_ * w
        if PFB:
            _Trms = capo.pfb.pfb(Trms, window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
            _Nrms = capo.pfb.pfb(Nrms, window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
            _Wrms = capo.pfb.pfb(w   , window=WINDOW, taps=NTAPS, fft=n.fft.ifft)
        else:
            window = a.dsp.gen_window(Trms.size, WINDOW)
            #window /= np.sum(window)
            _Trms = n.fft.ifft(window * Trms)
            _Nrms = n.fft.ifft(window * Nrms)
            _Wrms = n.fft.ifft(window * w)
        #gain = n.abs(_Wrms[0])
        #print 'Gain:', gain
        #if False and gain > 0: # XXX this inverts the blackmann harris out of the data completely
        #    _Trms.shape = (_Trms.size,1)
        #    C = n.zeros((_Trms.size, _Trms.size), dtype=n.complex)
        #    for k1 in xrange(_Wrms.size):
        #      for k2 in xrange(_Wrms.size):
        #        #C[k1,k2] = _Wrms[k2-k1]
        #        C[k1,k2] = _Wrms[k1-k2]
        #    _C = n.linalg.inv(C)
        #    _Trms = n.dot(_C, _Trms).squeeze()
        #    _Nrms = n.dot(_C, _Nrms).squeeze()
        _Trms = n.fft.fftshift(_Trms)
        _Nrms = n.fft.fftshift(_Nrms)
        if False: # swap in a simulated signal post delay transform
            _Trms = n.random.normal(size=_Trms.size) * n.exp(2j*n.pi*n.random.uniform(size=_Trms.size))
            mask = n.ones(_Trms.size); mask[15:25] = 0
            _Trms += .3*eor_mdl[times[-1]] * mask
        T[bl] = T.get(bl, []) + [_Trms]
        N[bl] = N.get(bl, []) + [_Nrms]

n_k = chans.size / NTAPS
bls = T.keys()
for bl in bls:
    T[bl],N[bl] = n.array(T[bl]),n.array(N[bl])
    print T[bl].shape
if True:
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
    print Ts.shape
    for i in xrange(Ts.shape[0]):
        if n.random.uniform() > .5:
            Ts[i] *= -1
            Ns[i] *= -1

print Ts.shape
print times[300], times[500]
print ' '.join(['%d_%d' % a.miriad.bl2ij(bl) for bl in bls])
if PLOT:
    #capo.arp.waterfall(cov(Ts), mode='log', drng=2); p.show()
    p.subplot(131); capo.arp.waterfall(Ts, mode='log', mx=1, drng=2); p.colorbar(shrink=.5)
    p.subplot(132); capo.arp.waterfall(Ns, mode='log', mx=1, drng=2); p.colorbar(shrink=.5)
    p.subplot(133); capo.arp.waterfall(cov(Ts), mode='log', drng=3); p.colorbar(shrink=.5)
    p.show()
    p.subplot(121); capo.arp.waterfall(cov(Ts), mode='real', mx=.005, drng=.01); p.colorbar(shrink=.5)
    p.subplot(122); capo.arp.waterfall(cov(Ts), mode='log', drng=3); p.colorbar(shrink=.5)
    p.show()

for boot in xrange(NBOOT):
    if True: # pick a sample of baselines with replacement
        #bls_ = [random.choice(bls) for bl in bls]
        bls_ = random.sample(bls, len(bls))
        #gp1,gp2 = bls_[:len(bls)/2],bls_[len(bls)/2:] # ensure gp1 and gp2 can't share baselines
        #gp1,gp2,gp3 = bls_[:4],bls_[4:9],bls_[9:]
        gp1,gp2,gp3,gp4 = bls_[:7],bls_[7:14],bls_[14:21],bls_[21:] # for 28bl
        #gp1,gp2,gp3,gp4 = bls_[:5],bls_[5:10],bls_[10:15],bls_[15:] # for 21bl
        # ensure each group has at least 2 kinds of baselines
        gp1 = random.sample(gp1, 2) + [random.choice(gp1) for bl in gp1[:len(gp1)-2]]
        gp2 = random.sample(gp2, 2) + [random.choice(gp2) for bl in gp2[:len(gp2)-2]]
        gp3 = random.sample(gp3, 2) + [random.choice(gp3) for bl in gp3[:len(gp3)-2]]
        gp4 = random.sample(gp4, 2) + [random.choice(gp4) for bl in gp4[:len(gp4)-2]]
    else:
        bls_ = random.sample(bls, len(bls))
        gp1,gp2 = bls_[:len(bls)/2],bls_[len(bls)/2:]
    #gp2 = gp2[:len(gp1)] # XXX force gp1 and gp2 to be same size
    #gp1,gp2 = gp1+gp2,[] # XXX
    #print 'XXX', len(gp1), len(gp2)
    bls_ = gp1 + gp2 + gp3 + gp4
    print 'Bootstrap sample %d:' % boot,
    for gp in [gp1,gp2,gp3,gp4]: print '(%s)' % (','.join(['%d_%d' % a.miriad.bl2ij(bl) for bl in gp])),
    print
    Ts = n.concatenate([T[bl] for bl in bls_], axis=1).T
    Ns = n.concatenate([N[bl] for bl in bls_], axis=1).T # this noise copy processed as if it were the data
    L = len(bls_)
    #temp_noise_var = n.var(n.array([T[bl] for bl in bls_]), axis=0).T
    temp_noise_var = n.average(n.array([T[bl] for bl in bls_]), axis=0).T
    print Ts.shape, temp_noise_var.shape

    _Cxtot,_Cntot = 1, 1
    #PLT1,PLT2 = 4,4
    PLT1,PLT2 = 3,3
    #PLT1,PLT2 = 1,2
    for cnt in xrange(PLT1*PLT2-1):
        print cnt, '/', PLT1*PLT2-1
        if PLOT:
            p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(cov(Ts), mode='log', mx=-1, drng=3)
            #p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(cov(Ns), mode='log', mx=0, drng=2)
        SZ = Ts.shape[0]
        Cx,Cn = cov(Ts), cov(Ns)
        for c in [Cx,Cn]: # Normalize covariance matrices
            d = n.diag(c); d.shape = (1,SZ)
            c /= n.sqrt(d) * 2
        #g = .3 # for 1*7 baselines
        g = .2 # for 4*7 baselines
        # begin with off-diagonal covariances to subtract off 
        # (psuedo-inv for limit of small off-diagonal component)
        _Cx,_Cn = -g*Cx, -g*Cn 
        ind = n.arange(SZ)
        # XXX do we also need to zero modes adjacent to diagonal, since removing them results in bigger signal loss?
        for b in xrange(L): # for each redundant baseline, zero out diagonal from covariance diagonalization
            indb = ind[:-b*n_k]
            _Cx[indb,indb+b*n_k] = _Cx[indb+b*n_k,indb] = 0
            _Cn[indb,indb+b*n_k] = _Cn[indb+b*n_k,indb] = 0
        _Cx[ind,ind] = _Cn[ind,ind] = 0 # set these to zero temporarily to avoid noise bias into cross terms
        if True: # estimate and remove signal covariance from diagonalization process 
            # do this twice: once for the signal (Cx) and once for the noise (Cn)
            # using the statistics of the signal and noise, respectively
            for _C in [_Cx,_Cn]:
                _C.shape = (L,n_k,L,n_k)
                sub_C = n.zeros_like(_C)
                # Choose a (i,j) baseline cross-multiple panel in the covariance matrix
                for i in xrange(L):
                    bli = bls_[i] 
                    for j in xrange(L):
                        blj = bls_[j]
                        # even remove signal bias if bli == blj
                        # ensure bli and blj belong to the same group
                        if bli in gp1 and blj in gp1: gp = gp1
                        elif bli in gp2 and blj in gp2: gp = gp2
                        elif bli in gp3 and blj in gp3: gp = gp3
                        elif bli in gp4 and blj in gp4: gp = gp4
                        else: continue # make sure we only compute average using bls in same group
                        # Now average over all other panels of covariance matrix (within this group)
                        # to get the average signal covariance and subtract that off so that we don't
                        # get signal loss removing residual signal covariances.
                        _Csum,_Cwgt = 0,0
                        for i_ in xrange(L):
                            bli_ = bls_[i_]
                            if not bli_ in gp: continue # XXX need to also bail if bli_ not in gp
                            if bli == bli_: continue # only average over other bls to better isolate bl systematics
                            for j_ in xrange(L):
                                blj_ = bls_[j_]
                                if not blj_ in gp: continue # XXX need to also bail if blj_ not in gp
                                if bli_ == blj_: continue # don't average over panels with noise bias
                                if blj == blj_: continue # only average over other bls to better isolate bl systematics
                                _Csum += _C[i_,:,j_] # fixes indexing error in earlier ver
                                _Cwgt += 1
                        sub_C[i,:,j] = _Csum / _Cwgt # fixes indexing error in earlier ver XXX careful if _Cwgt is 0
                _C.shape = sub_C.shape = (L*n_k,L*n_k)
                _C -= sub_C
        if True:
            # divide bls into two independent groups to avoid cross-contamination of noise
            # this is done by setting mask=0 for all panels pairing bls between different groups
            # this masks covariances between groups, so no covariance from a bl in one group is subtracted from
            # a bl in another group
            mask = n.ones_like(Cx)
            # XXX need to clean this section up
            #for bl1 in xrange(len(gp1)):
            #    for bl2 in xrange(len(gp1)):
            #        if bls_[bl1] != bls_[bl2]: continue # zero out panels where bl1 == bl2
            #        mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
            #        mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            for bl1 in xrange(len(gp1)):
                for bl2 in xrange(len(gp2)):
                    bl2 += len(gp1)
                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            for bl1 in xrange(len(gp1)):
                for bl2 in xrange(len(gp3)):
                    bl2 += len(gp1) + len(gp2)
                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            for bl1 in xrange(len(gp1)):
                for bl2 in xrange(len(gp4)):
                    bl2 += len(gp1) + len(gp2) + len(gp3)
                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            #for bl1 in xrange(len(gp2)):
            #    bl1 += len(gp1)
            #    for bl2 in xrange(len(gp2)):
            #        bl2 += len(gp1)
            #        if bls_[bl1] != bls_[bl2]: continue # zero out panels where bl1 == bl2
            #        mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
            #        mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            for bl1 in xrange(len(gp2)):
                bl1 += len(gp1)
                for bl2 in xrange(len(gp3)):
                    bl2 += len(gp1) + len(gp2)
                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            for bl1 in xrange(len(gp2)):
                bl1 += len(gp1)
                for bl2 in xrange(len(gp4)):
                    bl2 += len(gp1) + len(gp2) + len(gp3)
                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            #for bl1 in xrange(len(gp3)):
            #    bl1 += len(gp1) + len(gp2)
            #    for bl2 in xrange(len(gp3)):
            #        bl2 += len(gp1) + len(gp2)
            #        if bls_[bl1] != bls_[bl2]: continue # zero out panels where bl1 == bl2
            #        mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
            #        mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            for bl1 in xrange(len(gp3)):
                bl1 += len(gp1) + len(gp2)
                for bl2 in xrange(len(gp4)):
                    bl2 += len(gp1) + len(gp2) + len(gp3)
                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            #for bl1 in xrange(len(gp4)):
            #    bl1 += len(gp1) + len(gp2) + len(gp3)
            #    for bl2 in xrange(len(gp4)):
            #        bl2 += len(gp1) + len(gp2) + len(gp3)
            #        if bls_[bl1] != bls_[bl2]: continue # zero out panels where bl1 == bl2
            #        mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
            #        mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            _Cx *= mask; _Cn *= mask
        _Cx[ind,ind] = _Cn[ind,ind] = 1
        Ts,Ns = n.dot(_Cx,Ts), n.dot(_Cn,Ns)
        # These are a running tally of all the diagonalization steps applied
        _Cxtot,_Cntot = n.dot(_Cx,_Cxtot), n.dot(_Cn,_Cntot)
    if PLOT:
        p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ts), mode='log', mx=-1, drng=3)
        #p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ns), mode='log', mx=0, drng=3)
        p.show()

    Ts = n.concatenate([T[bl] for bl in bls_], axis=1).T
    Ns = n.concatenate([N[bl] for bl in bls_], axis=1).T # this noise copy processed as if it were the data

    pspecs,dspecs = [], []
    nspecs,n1specs,n2specs = [], [], []
    Cx,Cn = CoV(Ts, bls_), CoV(Ns, bls_)
    Cx_ = CoV(n.dot(_Cxtot,Ts), bls_)
    # Cn1 is the noise diagonalized as if it were the signal, Cn2 is the noise with the signal diagonalization applied
    Cn1_,Cn2_ = CoV(n.dot(_Cntot,Ns), bls_), CoV(n.dot(_Cxtot,Ns), bls_)
    for cnt,bli in enumerate(bls_):
        print cnt
        for blj in bls_[cnt:]:
            #print a.miriad.bl2ij(bli), a.miriad.bl2ij(blj)
            # XXX behavior here is poorly defined for repeat baselines in bootstrapping
            xi,xj = Cx.get_x(bli), Cx.get_x(blj)
            xi_,xj_ = Cx_.get_x(bli), Cx_.get_x(blj)
            pk_avg = scalar * xi * xj.conj()
            dspecs.append(pk_avg) # do this before bli == blj check to include noise bias in dspec
            if bli == blj: continue
            if True: # exclude intra-group pairings # XXX
                if (bli in gp1 and blj in gp1) or (bli in gp2 and blj in gp2) or (bli in gp3 and blj in gp3) or (bli in gp4 and blj in gp4): continue
            if True: # do an extra final removal of leakage from particular modes
                Ts = n.concatenate([xi_,xj_], axis=0)
                cx = cov(Ts)
                if PLOT:
                    p.clf()
                    p.subplot(121); capo.arp.waterfall(cx, mode='log', mx=-1, drng=3)
                for cnt1 in xrange(9):
                    d = n.diag(cx); d.shape = (1,d.size); cx /= n.sqrt(d) * 2
                    g = .3
                    _cx = -g*cx
                    mask = n.zeros_like(cx)
                    if n_k == 20: prj_ch = xrange(8,12)
                    elif n_k == 40: prj_ch = xrange(17,24)
                    elif n_k == 80: prj_ch = xrange(34,48)
                    else: raise ValueError('Only support # channels = (20,40,80) for now')
                    for k in prj_ch:
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
                    p.subplot(122); capo.arp.waterfall(cx, mode='log', mx=-1, drng=3)
                    p.show()
                xi_,xj_ = Ts[:n_k],Ts[n_k:]
            ni,nj = Cn.get_x(bli), Cn.get_x(blj)
            n1i_,n1j_ = Cn1_.get_x(bli), Cn1_.get_x(blj)
            n2i_,n2j_ = Cn2_.get_x(bli), Cn2_.get_x(blj)
            nij = n.sqrt(n.mean(n.abs(ni*nj.conj())**2, axis=1))
            n1ij_ = n.sqrt(n.mean(n.abs(n1i_*n1j_.conj())**2, axis=1))
            n2ij_ = n.sqrt(n.mean(n.abs(n2i_*n2j_.conj())**2, axis=1))
            f1 = n.sqrt(n.mean(n.abs(nij)**2)/n.mean(n.abs(n1ij_)**2))
            f2 = n.sqrt(n.mean(n.abs(nij)**2)/n.mean(n.abs(n2ij_)**2))
            #f1 = n.sqrt((n.abs(nij)**2)/(n.abs(n1ij_)**2))
            #f2 = n.sqrt((n.abs(nij)**2)/(n.abs(n2ij_)**2))
            print 'Rescale factor:', f1, f2
            #rescale = max(f1,f2)
            rescale = 1.

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
                p.plot(rescale*n.average(xi_*xj_.conj(), axis=1).real, 'g')
                p.subplot(224)
                p.plot(n.average(n1i_*n1i_.conj(), axis=1).real, 'k')
                p.plot(n.average(n1j_*n1j_.conj(), axis=1).real, 'b')
                p.plot(n.average(n1i_*n1j_.conj(), axis=1).real, 'g')
                p.plot(n.average(ni*nj.conj(), axis=1).real, 'r')
                p.plot(rescale*n.average(n1i_*n1j_.conj(), axis=1).real, 'c')
                p.show()
            elif False:
                import pylab as p
                blistr,bljstr = str(a.miriad.bl2ij(bli)),str(a.miriad.bl2ij(blj))
                print blistr, bljstr
                #p.subplot(121); p.plot(n.average(xi*xj.conj(), axis=1).real, label='%s-%s'%(blistr,bljstr))
                p.plot(n.average(xi_*xj_.conj(), axis=1).real, label='%s-%s'%(blistr,bljstr))
            pk_avg_ = scalar * xi_ * xj_.conj() * rescale# XXX
            pspecs.append(pk_avg_)
            nspecs.append(scalar * ni * nj.conj())
            n1specs.append(scalar * n1i_ * n1j_.conj())
            n2specs.append(scalar * n2i_ * n2j_.conj())
    #p.subplot(121); p.legend(loc='best')
    #p.legend(loc='best')
    pspecs,dspecs = n.array(pspecs), n.array(dspecs)
    nspecs,n1specs,n2specs = n.array(nspecs), n.array(n1specs), n.array(n2specs)
    navg_2d = n.average(nspecs, axis=0)
    n1avg_2d = n.average(n1specs, axis=0)
    n2avg_2d = n.average(n2specs, axis=0)
    f1 = n.sqrt(n.sum(n.abs(navg_2d)**2)/n.sum(n.abs(n1avg_2d)**2))
    f2 = n.sqrt(n.sum(n.abs(navg_2d)**2)/n.sum(n.abs(n2avg_2d)**2))
    print 'Rescale factor (FINAL):', f1, f2
    
    avg_2d = n.average(pspecs, axis=0) # average over baseline cross-multiples
    std_2d = n.std(pspecs, axis=0) # get distribution as a function of time
    wgt_2d = 1. / std_2d**2 # inverse variance weighting
    avg_1d = n.sum(avg_2d * wgt_2d, axis=1) / n.sum(wgt_2d, axis=1)

    if PLOT:
        import capo as C
        p.subplot(131)
        C.arp.waterfall(avg_2d, mode='log', drng=3)
        p.subplot(132)
        C.arp.waterfall(wgt_2d, mode='log', drng=3)
        p.subplot(133)
        C.arp.waterfall(avg_2d*wgt_2d, mode='log', drng=3)
        p.show()
    
    #avg_1d = n.average(dspecs, axis=0)
    #p.subplot(133)
    if PLOT: p.plot(avg_1d.real,'.')
    #p.plot(n.average(dspecs, axis=0).real/scalar)
    #p.show()
    #std_1d = n.std(pspecs, axis=0) / n.sqrt(pspecs.shape[0]) # coarse estimate of errors.  bootstrapping will do better
    std_1d = n.std(avg_2d, axis=1) / n.sqrt(pspecs.shape[0]) # coarse estimate of errors.  bootstrapping will do better
    #std_1d = n.std(pspecs, axis=0) # in new noise subtraction, this remaining dither is essentially a bootstrap error, but with 5/7 of the data

    print 'Writing pspec_boot%04d.npz' % boot
    n.savez('pspec_boot%04d.npz'%boot, kpl=kpl, pk=avg_1d, err=std_1d, scalar=scalar, times=n.array(times),
        pk_vs_t=avg_2d, err_vs_t=std_2d, temp_noise_var=temp_noise_var, nocov_vs_t=n.average(dspecs,axis=0),
        cmd=' '.join(sys.argv))
if PLOT: p.show()

