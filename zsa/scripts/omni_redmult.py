#! /usr/bin/env python
import aipy as a, numpy as n, capo 
import optparse, sys, random
import pylab as p

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

o = optparse.OptionParser()
a.scripting.add_standard_options(o, chan=True)
o.add_option('--sep', dest='sep',
    help='Separation to use.')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
opts,args = o.parse_args(sys.argv[1:])

freqs = n.linspace(.1, .2, 203)
chans = a.scripting.parse_chans(opts.chan, freqs.size)
WINDOW = opts.window
NBOOT = 20
sep = opts.sep

afreqs = freqs.take(chans)
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)
sdf = freqs[1]-freqs[0]
n_k = afreqs.size

TRNG = (300,800)
t = n.arange(-200,200,1) * 42.8
w = a.dsp.gen_window(t.size, 'none')
sig = .000288*(fq/.1788)
cen = .001059 * (fq/.1788)
fir = n.exp(-(t**2)/(2*(1./sig)**2)).astype(n.complex) * w
fir /= fir.sum()
fir *= n.exp(2j*n.pi*cen*t) # need to flip the sign for backward baselines
#fir = n.array([1.])
p.plot(fir.real)
p.plot(fir.imag)
p.show()

B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW]
etas = n.fft.fftshift(capo.pspec.f2eta(afreqs))
kpl = etas * capo.pspec.dk_deta(z)
bm = n.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq)
scalar = capo.pspec.X2Y(z) * bm * B

print 'Freq:',fq
print 'z:', z
print 'B:', B
print 'scalar:', scalar
print 'sep:', sep

_T, _W = {}, {}
for filename in args:
    print 'Reading', filename
    npz = n.load(filename)
    d = npz[sep].take(chans, axis=1)
    try: w = npz['wgt_'+sep].take(chans, axis=1) / 1e-6
    except(KeyError): w = n.ones_like(d)
    d = n.concatenate([d[-200:],d[:1460]], axis=0)
    w = n.concatenate([w[-200:],w[:1460]], axis=0)
    #w = n.where(w > .1, 1., 0)
    for ch in xrange(d.shape[1]):
        d[:,ch] = n.convolve(w[:,ch]*d[:,ch], fir, mode='same')
        w[:,ch] = n.convolve(w[:,ch], n.abs(fir), mode='same')
        d[:,ch] /= n.where(w[:,ch] > 0, w[:,ch], 1)
    w = n.where(w > .1, 1, 0)
    #w = n.average(w, axis=1).reshape((w.shape[0],1)) * n.ones((1,w.shape[1]))
    #p.subplot(131); capo.arp.waterfall(d, mode='phs')
    #p.subplot(132); capo.arp.waterfall(w)
    d *= w
    #capo.arp.waterfall(n.fft.fft(d[300:1000], axis=0), mode='lin')
    #p.colorbar()
    #p.show()
    T = d * capo.pspec.jy2T(afreqs)
    window = a.dsp.gen_window(T.shape[1], WINDOW); window.shape = (1,window.size)
    _T[filename] = n.fft.fftshift(n.fft.ifft(window * T), axes=1)
    #_W[filename] = n.fft.fftshift(n.fft.ifft(window * w), axes=1)
    _W[filename] = n.fft.ifft(window * w)

    p.subplot(121); capo.arp.waterfall(_T[filename], drng=3)
    p.subplot(122); capo.arp.waterfall(_W[filename], drng=3)
    p.show()

files = _T.keys(); files.sort()
L = len(files)
for boot in xrange(NBOOT):
    files_ = random.sample(files, len(files))
    #gp1,gp2,gp3 = files[:3],files[3:6],files[6:]
    gp1,gp2 = files[:5],files[5:]
    gp1 = random.sample(gp1, 2) + [random.choice(gp1) for bl in gp1[:len(gp1)-2]]
    gp2 = random.sample(gp2, 2) + [random.choice(gp2) for bl in gp2[:len(gp2)-2]]
    files_ = gp1 + gp2
    print 'Bootstrap sample %d:' % boot,
    for gp in [gp1,gp2]: print '(%s)' % (','.join([f for f in gp])),
    print
    Ts = n.concatenate([_T[f] for f in files_], axis=1).T

    PLOT = False
    _Cxtot = 1
    PLT1,PLT2 = 5,5
    for cnt in xrange(PLT1*PLT2-1):
        print cnt, '/', PLT1*PLT2-1
        if PLOT:
            p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(cov(Ts), mode='log', mx=-3, drng=3)
            #p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(cov(Ns), mode='log', mx=0, drng=2)
        SZ = Ts.shape[0]
        Cx = cov(Ts)
        for c in [Cx]: # Normalize covariance matrices
            d = n.diag(c); d.shape = (1,SZ)
            c /= n.sqrt(d) * 2
        g = .3 # for 1*7 baselines
        #g = .2 # for 4*7 baselines
        # begin with off-diagonal covariances to subtract off 
        # (psuedo-inv for limit of small off-diagonal component)
        _Cx = -g*Cx
        ind = n.arange(SZ)
        for b in xrange(L): # for each redundant baseline, zero out diagonal from covariance diagonalization
            indb = ind[:-b*n_k]
            _Cx[indb,indb+b*n_k] = _Cx[indb+b*n_k,indb] = 0
        _Cx[ind,ind] = 0 # set these to zero temporarily to avoid noise bias into cross terms
        if True: # estimate and remove signal covariance from diagonalization process 
            # do this twice: once for the signal (Cx) and once for the noise (Cn)
            # using the statistics of the signal and noise, respectively
            for _C in [_Cx]:
                _C.shape = (L,n_k,L,n_k)
                sub_C = n.zeros_like(_C)
                # Choose a (i,j) baseline cross-multiple panel in the covariance matrix
                for i in xrange(L):
                    bli = files_[i]
                    for j in xrange(L):
                        blj = files_[j]
                        # even remove signal bias if bli == blj
                        # ensure bli and blj belong to the same group
                        if bli in gp1 and blj in gp1: gp = gp1
                        elif bli in gp2 and blj in gp2: gp = gp2
                        else: continue # make sure we only compute average using bls in same group
                        # to get the average signal covariance and subtract that off so that we don't
                        # get signal loss removing residual signal covariances.
                        _Csum,_Cwgt = 0,0
                        for i_ in xrange(L):
                            bli_ = files_[i_]
                            if not bli_ in gp: continue # XXX need to also bail if bli_ not in gp
                            if bli == bli_: continue # only average over other bls to better isolate bl systematics
                            for j_ in xrange(L):
                                blj_ = files_[j_]
                                if not blj_ in gp: continue # XXX need to also bail if blj_ not in gp
                                if bli_ == blj_: continue # don't average over panels with noise bias
                                if blj == blj_: continue # only average over other bls to better isolate bl systematics
                                _Csum += _C[i_,:,j_] # fixes indexing error in earlier ver
                                _Cwgt += 1
                        try: sub_C[i,:,j] = _Csum / _Cwgt # fixes indexing error in earlier ver XXX careful if _Cwgt is 0
                        except(ZeroDivisionError): continue
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
            #        if files[bl1] != files[bl2]: continue # zero out panels where bl1 == bl2
            #        mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
            #        mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            for bl1 in xrange(len(gp1)):
                for bl2 in xrange(len(gp2)):
                    bl2 += len(gp1)
                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
            _Cx *= mask
        _Cx[ind,ind] = 1
        Ts = n.dot(_Cx,Ts)
        # These are a running tally of all the diagonalization steps applied
        _Cxtot = n.dot(_Cx,_Cxtot)
    if PLOT:
        p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ts), mode='log', mx=-3, drng=3)
        #p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ns), mode='log', mx=0, drng=3)
        p.show()

    for i,f in enumerate(files_):
        new_T = Ts[i*n_k:(i+1)*n_k].T
        #p.subplot(121); capo.arp.waterfall(_T[f], drng=3)
        #p.subplot(122); capo.arp.waterfall(new_T, drng=3)
        #p.show()
        _T[f] = new_T

    #for i,filename in enumerate(files_):
    #    print i, filename
    #    p.subplot(1, len(_T), i+1)
    #    capo.arp.waterfall(_T[filename], mode='log')
    #    p.colorbar(shrink=.5)
    #p.show()

    pk_sum, pk_wgt = 0., 0.
    #for i,f1 in enumerate(files):
    for f1 in gp1:
        _T1 = _T[f1]
        _W1 = _W[f1]
        #for f2 in _T.keys()[i+1:]:
        for f2 in gp2:
            _T2 = _T[f2]
            _W2 = _W[f2]
            pk12 = scalar * _T1 * n.conj(_T2)
            wt12 = _W1 * _W2
            wt12 /= wt12.max()
            pk_sum += pk12
            #p.subplot(121)
            #p.plot(kpl, n.average(pk12[TRNG[0]:TRNG[1]], axis=0).real / n.average(wt12[TRNG[0]:TRNG[1],:1], axis=0))
            #pk_wgt += 1
            #p.clf(); p.plot(wt12[:,0]); p.show()
            pk_wgt += wt12[:,:1]

    #p.subplot(121)
    #pk = pk_sum / pk_wgt
    #capo.arp.waterfall(pk, mode='real', mx=1e9, drng=2e9)#, drng=3)
    #p.colorbar()

    pk_avg = pk_sum[TRNG[0]:TRNG[1]].sum(axis=0) / pk_wgt[TRNG[0]:TRNG[1]].sum(axis=0)
    print pk_avg.shape
    
    p.subplot(121)
    p.plot(kpl, pk_avg.real, 'k^-')

    p.subplot(122)
    p.plot(n.abs(kpl), 2* n.abs(kpl)**3/(2*n.pi**2) * pk_avg.real, 'k.-')
p.show()

_T = n.concatenate([_T[f] for f in _T.keys()], axis=1)

_T = _T[TRNG[0]:TRNG[1]].T
capo.arp.waterfall(cov(_T), mode='log')
p.show()


