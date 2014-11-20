#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo
import glob, optparse, sys, random

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
o.add_option('-b', '--nboot', type='int', default=20,
    help='Number of bootstraps.  Default is 20')
o.add_option('--plot', action='store_true',
    help='Generate plots')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
opts,args = o.parse_args(sys.argv[1:])

random.seed(0)
POL = 'I'
LST_STATS = False
DELAY = False
NGPS = 5
INJECT_SIG = False
SAMPLE_WITH_REPLACEMENT = True
NOISE = .0
PLOT = opts.plot

def get_data(filenames, antstr, polstr, verbose=False):
    # XXX could have this only pull channels of interest to save memory
    lsts, dat, flg = [], {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            lst = uv['lst']
            if len(lsts) == 0 or lst != lsts[-1]: lsts.append(lst)
            bl = a.miriad.ij2bl(i,j)
            if not dat.has_key(bl): dat[bl],flg[bl] = [],[]
            dat[bl].append(d)
            flg[bl].append(f)
            #if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            #pol = a.miriad.pol2str[uv['pol']]
            #if not dat[bl].has_key(pol):
            #    dat[bl][pol],flg[bl][pol] = [],[]
            #dat[bl][pol].append(d)
            #flg[bl][pol].append(f)
    return n.array(lsts), dat, flg

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

def noise(size):
    return n.random.normal(size=size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=size))

def get_Q(mode, n_k):
    if not DELAY:
        _m = n.zeros((n_k,), dtype=n.complex)
        _m[mode] = 1.
        m = n.fft.fft(n.fft.ifftshift(_m)) * a.dsp.gen_window(nchan, WINDOW)
        Q = n.einsum('i,j', m, m.conj())
        return Q
    else:
        # XXX need to have this depend on window
        Q = n.zeros_like(C)
        Q[mode,mode] = 1
        return Q

dsets = {
    'even': glob.glob('even/sep0,1/*242.[3456]*uvAL'),
    'odd' : glob.glob('odd/sep0,1/*243.[3456]*uvAL'),
}

WINDOW = opts.window
uv = a.miriad.UV(dsets.values()[0][0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
del(uv)

afreqs = freqs.take(chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

aa = a.cal.get_aa(opts.cal, n.array([.150]))
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
jy2T = capo.pspec.jy2T(afreqs)
window = a.dsp.gen_window(nchan, WINDOW)
#if not WINDOW == 'none': window.shape=(1,nchan)
if not WINDOW == 'none': window.shape=(nchan,1)

#B = sdf * afreqs.size / capo.pfb.NOISE_EQUIV_BW[WINDOW] # this is wrong if we aren't inverting
# the window post delay transform (or at least dividing out by the gain of the window)
# For windowed data, the FFT divides out by the full bandwidth, B, which is
# then squared.  Proper normalization is to multiply by B**2 / (B / NoiseEqBand) = B * NoiseEqBand
# XXX NEED TO FIGURE OUT BW NORMALIZATION
B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW] # normalization. See above.
etas = n.fft.fftshift(capo.pspec.f2eta(afreqs)) #create etas (fourier dual to frequency)
#etas = capo.pspec.f2eta(afreqs) #create etas (fourier dual to frequency)
kpl = etas * capo.pspec.dk_deta(z) #111
print kpl

if True:
    bm = n.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq) * 2.35 # correction for beam^2
    scalar = capo.pspec.X2Y(z) * bm * B
else: scalar = 1
if not DELAY:
    # XXX this is a hack
    if WINDOW == 'hamming': scalar /= 3.67
    elif WINDOW == 'blackman-harris': scalar /= 5.72
print 'Freq:',fq
print 'z:', z
print 'B:', B
print 'scalar:', scalar
sys.stdout.flush()

# acquire the data
#antstr = '41_49,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_47,3_25,1_58,35_61,42_63,2_33'
#antstr = '9_58,42_63'#41_49,3_10,9_58,22_61,20_63'#,2_43,21_53,31_45,41_47,3_25,1_58,35_61,42_63,2_33'
#antstr = '9_58,42_63,2_33,22_61,35_61,28_34,39_46,23_30'
antstr = 'cross'
lsts,data,flgs = {},{},{}
for k in dsets:
    lsts[k],data[k],flgs[k] = get_data(dsets[k], antstr=antstr, polstr=POL, verbose=True)

if LST_STATS:
    # collect some metadata from the lst binning process
    cnt, var = {}, {}
    for filename in dsets.values[0]:
        print 'Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, '41_49', POL)
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            bl = '%d,%d,%d' % (i,j,uv['pol'])
            cnt[bl] = cnt.get(bl, []) + [uv['cnt']]
            var[bl] = var.get(bl, []) + [uv['var']]
    cnt = n.array(cnt.values()[0]) # all baselines should be the same
    var = n.array(var.values()[0]) # all baselines should be the same
else: cnt,var = n.ones_like(lsts.values()[0]), n.ones_like(lsts.values()[0])


# Align data sets in LST
lstmax = max([lsts[k][0] for k in lsts])
for k in lsts:
    print k
    for i in xrange(len(lsts[k])):
        # allow for small numerical differences (which shouldn't exist!)
        if lsts[k][i] >= lstmax - .001: break
    lsts[k] = lsts[k][i:]
    for bl in data[k]:
        data[k][bl],flgs[k][bl] = data[k][bl][i:],flgs[k][bl][i:]
j = min([len(lsts[k]) for k in lsts])
for k in lsts:
    lsts[k] = lsts[k][:j]
    for bl in data[k]:
        data[k][bl],flgs[k][bl] = n.array(data[k][bl][:j]),n.array(flgs[k][bl][:j])
lsts = lsts.values()[0]

d1,d2 = {},{}
x = {}
for k in data:
    x[k] = {}
    for bl in data[k]:
        d = data[k][bl][:,chans] * jy2T
        if conj[bl]: d = n.conj(d)
        x[k][bl] = n.transpose(d, [1,0]) # swap time and freq axes
        
bls_master = x.values()[0].keys()
nbls = len(bls_master)
print 'Baselines:', nbls

if INJECT_SIG > 0.: # Create a fake EoR signal to inject
    print 'INJECTING SIMULATED SIGNAL'
    eor = noise(x.values()[bls_master[0]].shape) * INJECT_SIG
    fringe_filter = n.ones((44,))
    # Maintain amplitude of original noise
    fringe_filter /= n.sqrt(n.sum(fringe_filter))
    for ch in xrange(eor.shape[0]):
        eor[ch] = n.convolve(eor[ch], fringe_filter, mode='same')
    _eor = n.fft.ifft(eor, axis=0)
    #wgt = n.exp(-n.fft.ifftshift(kpl)**2/(2*.3**2))
    wgt = n.zeros(_eor.shape[0]); wgt[0] = 1
    wgt.shape = wgt.shape + (1,)
    #_eor *= wgt
    #_eor = n.fft.ifft(eor, axis=0); _eor[4:-3] = 0
    #eor = n.fft.fft(_eor, axis=0)
    eor *= wgt
    for k in x:
        for bl in x[k]: x[k][bl] += eor
    if PLOT:
        capo.arp.waterfall(eor, mode='real'); p.colorbar(); p.show()

Q = {} # Create the Q's that extract power spectrum modes
for i in xrange(nchan):
    print 'Q',i
    Q[i] = get_Q(i, nchan)

# Compute baseline auto-covariances and apply inverse to data
I,_I,_Ix = {},{},{}
C,_C,_Cx = {},{},{}
for k in x:
    I[k],_I[k],_Ix[k] = {},{},{}
    C[k],_C[k],_Cx[k] = {},{},{}
    for bl in x[k]:
        C[k][bl] = cov(x[k][bl])
        I[k][bl] = n.identity(C[k][bl].shape[0])
        U,S,V = n.linalg.svd(C[k][bl].conj())
        _S = 1./S
        #_S[10:] = _S[10-1] # XXX
        _C[k][bl] = n.einsum('ij,j,jk', V.T, _S, U.T)
        _I[k][bl] = n.identity(_C[k][bl].shape[0])
        _Cx[k][bl] = n.dot(_C[k][bl], x[k][bl])
        _Ix[k][bl] = x[k][bl].copy()
        if PLOT and False:
            #p.plot(S); p.show()
            p.subplot(311); capo.arp.waterfall(x[k][bl], mode='real')
            p.subplot(323); capo.arp.waterfall(C[k][bl])
            p.subplot(324); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
            p.subplot(313); capo.arp.waterfall(_Cx[k][bl], mode='real')
            p.show()
        


for boot in xrange(opts.nboot):
    print '%d / %d' % (boot+1,opts.nboot)
    bls = bls_master[:]
    if not SAMPLE_WITH_REPLACEMENT:
        random.shuffle(bls)
        bls = bls[:-5] # XXX
    else: # sample with replacement
        bls = [random.choice(bls) for bl in bls]
    nbls = len(bls)
    gps = [bls[i::NGPS] for i in range(NGPS)]
    gps = [[random.choice(gp) for bl in gp] for gp in gps]
    #gps = [[bl for bl in gp] for gp in gps]
    bls = [bl for gp in gps for bl in gp]
    print '\n'.join([','.join(['%d_%d'%a.miriad.bl2ij(bl) for bl in gp]) for gp in gps])
    _Iz,_Isum,_IsumQ = {},{},{}
    _Cz,_Csum,_CsumQ = {},{},{}
    for k in _Cx:
        _Iz[k],_Isum[k],_IsumQ[k] = {},{},{}
        _Cz[k],_Csum[k],_CsumQ[k] = {},{},{}
        for i,gp in enumerate(gps):
            _Iz[k][i] = sum([_Ix[k][bl] for bl in gp])
            _Cz[k][i] = sum([_Cx[k][bl] for bl in gp])
            _Isum[k][i] = sum([_I[k][bl] for bl in gp])
            _Csum[k][i] = sum([_C[k][bl] for bl in gp])
            _IsumQ[k][i] = {}
            _CsumQ[k][i] = {}
            if DELAY: # this is much faster
                _Iz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Iz[k][i], axis=0), axes=0)
                _Cz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Cz[k][i], axis=0), axes=0)
                # XXX need to take fft of _Csum, _Isum here
            for ch in xrange(nchan): # XXX this loop makes computation go as nchan^3
                _IsumQ[k][i][ch] = n.dot(_Isum[k][i], Q[ch])
                _CsumQ[k][i][ch] = n.dot(_Csum[k][i], Q[ch])
        if PLOT:
            _Csumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            _Isumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            for i in xrange(len(gps)): _Isumk[i,:,i,:] = _Isum[k][i]
            _Isumk.shape = (NGPS*nchan, NGPS*nchan)
            _Isum[k] = _Isumk
            for i in xrange(len(gps)): _Csumk[i,:,i,:] = _Csum[k][i]
            _Csumk.shape = (NGPS*nchan, NGPS*nchan)
            _Csum[k] = _Csumk
            _Czk = n.array([_Cz[k][i] for i in _Cz[k]])
            _Czk = n.reshape(_Czk, (_Czk.shape[0]*_Czk.shape[1], _Czk.shape[2]))
            p.subplot(211); capo.arp.waterfall(_Czk, mode='real')
            p.subplot(223); capo.arp.waterfall(_Csum[k])
            p.subplot(224); capo.arp.waterfall(cov(_Czk))
            p.show()

    FI = n.zeros((nchan,nchan), dtype=n.complex)
    FC = n.zeros((nchan,nchan), dtype=n.complex)
    qI = n.zeros((nchan,_Iz.values()[0].values()[0].shape[1]), dtype=n.complex)
    qC = n.zeros((nchan,_Cz.values()[0].values()[0].shape[1]), dtype=n.complex)
    for k1 in _Cz:
        for k2 in _Cz:
            if k1 == k2: continue
            for bl1 in _Cz[k1]:
                for bl2 in _Cz[k2]:
                    if bl1 == bl2: continue
                    print k1, k2, bl1, bl2
                    if PLOT and False:
                        p.subplot(231); capo.arp.waterfall(C[m], drng=3)
                        p.subplot(232); capo.arp.waterfall(_C[m], drng=3)
                        p.subplot(233); capo.arp.waterfall(n.dot(C[m],_C[m]), drng=3)
                        p.subplot(234); p.semilogy(S)
                        p.subplot(236); capo.arp.waterfall(V, drng=3)
                        p.show()
                        p.subplot(311); capo.arp.waterfall(x[m], mode='real', mx=5, drng=10); p.colorbar(shrink=.5)
                        p.subplot(312); capo.arp.waterfall(_Cx, mode='real'); p.colorbar(shrink=.5)
                        p.subplot(313); capo.arp.waterfall(_Ix, mode='real'); p.colorbar(shrink=.5)
                        p.show()

                    #bC = n.zeros((nchan,1), dtype=n.complex)
                    #bI = n.zeros((nchan,1), dtype=n.complex)
                    #_CN = n.dot(_C,N)
                    if False: # use ffts to do q estimation fast
                        qI += n.conj(_Iz[k1][bl1]) * _Iz[k2][bl2]
                        qC += n.conj(_Cz[k1][bl1]) * _Cz[k2][bl2]
                    else: # brute force with Q to ensure normalization
                        _qI = n.array([_Iz[k1][bl1].conj() * n.dot(Q[i], _Iz[k2][bl2]) for i in xrange(nchan)])
                        qI += n.sum(_qI, axis=1)
                        _qC = n.array([_Cz[k1][bl1].conj() * n.dot(Q[i], _Cz[k2][bl2]) for i in xrange(nchan)])
                        qC += n.sum(_qC, axis=1)
                    for i in xrange(nchan):
                        for j in xrange(nchan):
                            FI[i,j] += n.einsum('ij,ji', _IsumQ[k1][bl1][i], _IsumQ[k2][bl2][j])
                            FC[i,j] += n.einsum('ij,ji', _CsumQ[k1][bl1][i], _CsumQ[k2][bl2][j])
                        #bI [i,0] = n.einsum('ij,ji', E[i], N)
                        #bC[i,0] = n.einsum('ij,ji', _CE[i], _CN)

        if PLOT:
            p.subplot(121); capo.arp.waterfall(FC, drng=4)
            p.subplot(122); capo.arp.waterfall(FI, drng=4)
            p.show()

        print 'Psuedoinverse of FC'
        
        # Other choices for M
        #U,S,V = n.linalg.svd(FC.conj())
        #_S = n.sqrt(1./S)
        # _S = 1./S
        # _S = n.ones_like(S)
        #MC = n.dot(n.transpose(V), n.dot(n.diag(_S), n.transpose(U)))
        #order = n.array([10,11,9,12,8,13,7,14,6,15,5,16,4,17,3,18,2,19,1,20,0])

        # Cholesky decomposition
        order = n.array([10,11,9,12,8,20,0,13,7,14,6,15,5,16,4,17,3,18,2,19,1])
        iorder = n.argsort(order)
        FC_o = n.take(n.take(FC,order, axis=0), order, axis=1)
        L_o = n.linalg.cholesky(FC_o)
        U,S,V = n.linalg.svd(L_o.conj())
        _S = 1./S
        MC_o = n.dot(n.transpose(V), n.dot(n.diag(_S), n.transpose(U)))
        MC = n.take(n.take(MC_o,iorder, axis=0), iorder, axis=1)
        MI  = n.identity(nchan, dtype=n.complex128)
        
        print 'Normalizing M/W'
        WI = n.dot(MI, FI)
        norm  = WI.sum(axis=-1); norm.shape += (1,)
        #norm  = WI.max(axis=-1); norm.shape += (1,) # XXX
        MI /= norm; WI = n.dot(MI, FI)
        WC = n.dot(MC, FC)
        norm  = WC.sum(axis=-1); norm.shape += (1,)
        #norm  = WC.max(axis=-1); norm.shape += (1,) # XXX
        MC /= norm; WC = n.dot(MC, FC)

        print 'Generating ps'
        pC = n.dot(MC, qC) * scalar
        #pC[m] *= 1.81 # signal loss, high-SNR XXX
        #pC[m] *= 1.25 # signal loss, low-SNR XXX
        pI = n.dot(MI, qI) * scalar

        if PLOT:
            p.subplot(411); capo.arp.waterfall(qC, mode='real'); p.colorbar(shrink=.5)
            p.subplot(412); capo.arp.waterfall(pC, mode='real'); p.colorbar(shrink=.5)
            p.subplot(413); capo.arp.waterfall(qI, mode='real'); p.colorbar(shrink=.5)
            p.subplot(414); capo.arp.waterfall(pI, mode='real'); p.colorbar(shrink=.5)
            p.show()

    if PLOT:
        p.plot(kpl, n.average(pC.real, axis=1), 'b.-')
        p.plot(kpl, n.average(pI.real, axis=1), 'k.-')
        p.show()

    print 'Writing pspec_boot%04d.npz' % boot
    n.savez('pspec_boot%04d.npz'%boot, kpl=kpl, scalar=scalar, times=n.array(lsts),
        pk_vs_t=pC, err_vs_t=1./cnt, temp_noise_var=var, nocov_vs_t=pI,
        cmd=' '.join(sys.argv))


