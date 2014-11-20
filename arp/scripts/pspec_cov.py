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
MASK = True
DELAY = False
NGPS = 51
LST_STATS = False
INJECT_SIG = False
NOISE = .0
JACKNIFE = True

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

def tile_panels(panel, bls, zero_block_diag=False, zero_repeat_bls=False):
    nbls = len(bls)
    n_k,n_k = panel.shape
    C = n.zeros((nbls,n_k,nbls,n_k), dtype=panel.dtype)
    for i in xrange(nbls):
        for j in xrange(nbls):
            if zero_repeat_bls and bls[i] == bls[j]: continue
            C[i,:,j,:] = panel.copy()
    if zero_block_diag:
        for i in xrange(nbls): C[i,:,i,:] = 0
    C.shape = (nbls*n_k, nbls*n_k)
    return C

def get_Q(mode, n_k, bls, zero_block_diag=False, zero_repeat_bls=False):
    if not DELAY:
        _m = n.zeros((n_k,), dtype=n.complex)
        _m[mode] = 1.
        m = n.fft.fft(n.fft.ifftshift(_m)) * a.dsp.gen_window(nchan, WINDOW)
        E = n.einsum('i,j', m, m.conj())
        return tile_panels(E, bls, zero_block_diag=zero_block_diag, zero_repeat_bls=zero_repeat_bls)
    else:
        # XXX need to integrate across waveband, not just mode at center point
        E = n.zeros_like(C)
        E[mode,mode] = 1
        return E

def auto_cov(d, n_k):
    nbls = d.shape[0]/n_k
    orig_shape = d.shape
    C = n.zeros((nbls,n_k,nbls,n_k), dtype=d.dtype)
    d.shape = (nbls,n_k,d.shape[1])
    for i in xrange(nbls): C[i,:,i,:] = cov(d[i])
    d.shape = orig_shape
    C.shape = (nbls*n_k, nbls*n_k)
    return C


files1 = glob.glob('even/sep0,1/*242.[3456]*uvAL') # XXX
files2 = glob.glob('odd/sep0,1/*243.[3456]*uvAL') # XXX
WINDOW = opts.window
uv = a.miriad.UV(args[0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
del(uv)

afreqs = freqs.take(chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

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
bm = n.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq) * 2.35 # correction for beam^2

if True: scalar = capo.pspec.X2Y(z) * bm * B
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
times1,data1,flgs1 = capo.arp.get_dict_of_uv_data(files1, antstr=antstr, polstr='I', verbose=True)
times2,data2,flgs2 = capo.arp.get_dict_of_uv_data(files2, antstr=antstr, polstr='I', verbose=True)

if LST_STATS:
    # collect some metadata from the lst binning process
    cnt, var = {}, {}
    for filename in args:
        print 'Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            bl = '%d,%d,%d' % (i,j,uv['pol'])
            cnt[bl] = cnt.get(bl, []) + [uv['cnt']]
            var[bl] = var.get(bl, []) + [uv['var']]
    cnt = n.array(cnt.values()[0]) # all baselines should be the same
    var = n.array(var.values()[0]) # all baselines should be the same
else: cnt,var = n.ones_like(times), n.ones_like(times)

aa = a.cal.get_aa(opts.cal, n.array([.150]))
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
jy2T = capo.pspec.jy2T(afreqs)
window = a.dsp.gen_window(nchan, WINDOW)
if not WINDOW == 'none': window.shape=(1,nchan)
d1,d2 = {},{}
for data,d in zip([data1,data2],[d1,d2]):
    for k in data:
        data[k]['I'] = data[k]['I'][:,chans] * jy2T
        if conj[k]: data[k]['I'] = n.conj(data[k]['I'])
        if DELAY: _d[k] = n.fft.fftshift(n.fft.ifft(window*data[k]['I']), axes=1)
        #else: d[k] = data[k]['I']
        else: d[k] = window*data[k]['I']
        d[k] = n.transpose(d[k], [1,0])
bls_master = d1.keys()
nbls = len(bls_master)
print 'Baselines:', nbls

# Create a fake EoR signal to inject
eor = noise(_d[bls_master[0]].shape) * .5
fringe_filter = n.ones((44,))
for ch in xrange(eor.shape[0]):
    eor[ch] = n.convolve(eor[ch], fringe_filter, mode='same')
_eor = n.fft.ifft(eor, axis=0); _eor[4:-3] = 0
eor = n.fft.fft(_eor, axis=0)

# Create the Q's that extract power spectrum modes
Q = {}
for i in xrange(nchan):
    print 'Q',i
    E[i] = get_Q(i, nchan, bls)
    #E[i] = get_Q(i, nchan, bls, zero_repeat_bls=NOAUTOS) * (1 - mask)


for boot in xrange(opts.nboot):
    print '%d / %d' % (boot+1,opts.nboot)
    bls = bls_master[:]
    random.shuffle(bls)
    bls = bls[:-5] # XXX
    nbls = len(bls)
    NGPS = nbls # XXX
    gps = [bls[i::NGPS] for i in range(NGPS)]
    gps = [[random.choice(gp) for bl in gp] for gp in gps]
    #gps = [[bl for bl in gp] for gp in gps]
    bls = [bl for gp in gps for bl in gp]
    print '\n'.join([','.join(['%d_%d'%a.miriad.bl2ij(bl) for bl in gp]) for gp in gps])

    if MASK: # mask covariance matrix
        mask = n.ones((nbls,nchan,nbls,nchan))
        for i in n.cumsum([len(gp) for gp in gps])[:-1]:
            mask[:i,:,i:,:] = 0
            mask[i:,:,:i,:] = 0
        mask.shape = (nbls*nchan,nbls*nchan)
        for i in xrange(NCHAN): Q[i] *= (1-mask) # XXX maybe try to avoid doing this every loop
    else: mask = 1

    if not INJECT_SIG:
        x1,x2 = n.array([d1[k] for k in bls]), n.array([d2[k] for k in bls])
    else:
        print 'INJECTING SIMULATED SIGNAL'
        x1,x2 = n.array([d1[k]+eor for k in bls]), n.array([d2[k]+eor for k in bls])
    x1 = n.reshape(x1, (x1.shape[0]*x1.shape[1], x1.shape[2]))
    x2 = n.reshape(x2, (x2.shape[0]*x2.shape[1], x2.shape[2]))
    x = {'s':0.5*(x1+x2), 'd':0.5*(x1-x2)}
    if PLOT:
        p.subplot(411); capo.arp.waterfall(x1, mode='real', mx=10, drng=20)
        p.subplot(412); capo.arp.waterfall(x2, mode='real', mx=10, drng=20)
        p.subplot(413); capo.arp.waterfall(x['s'], mode='real', mx=10, drng=20)
        p.subplot(414); capo.arp.waterfall(x['d'], mode='real', mx=10, drng=20)
        p.show()
    C,_C = {}, {}
    FC, FI = {}, {}
    MC, MI = {}, {}
    WC, WI = {}, {}
    pC, pI = {}, {}
    for m in 'sd':
        print 'Mode:', m
        C[m] = auto_cov(x[m] + NOISE*noise(_data.shape), nchan)
        #Cs *= mask; Cd *= mask # XXX probably not going to use this anymore

        print 'Psuedoinverse of C'
        U,S,V = n.linalg.svd(C[m].conj())
        #_S = n.where(S > 1e-2, 1./S, 0)
        _S = 1./S
        _C[m] = n.einsum('ij,j,jk', V.T, _S, U.T)
        _Cx = n.dot(_C[m], x[m])
        if PLOT:
            p.subplot(231); capo.arp.waterfall(C[m], drng=3)
            p.subplot(232); capo.arp.waterfall(_C[m], drng=3)
            p.subplot(233); capo.arp.waterfall(n.dot(C[m],_C[m]), drng=3)
            p.subplot(234); p.semilogy(S)
            p.subplot(236); capo.arp.waterfall(V, drng=3)
            p.show()
            p.subplot(211); capo.arp.waterfall(x[m], mode='real', mx=5, drng=10); p.colorbar(shrink=.5)
            p.subplot(212); capo.arp.waterfall(_Cx, mode='real'); p.colorbar(shrink=.5)
            p.show()

        #bC = n.zeros((nchan,1), dtype=n.complex)
        #bI = n.zeros((nchan,1), dtype=n.complex)
        #_CN = n.dot(_C,N)
        FC[m] = n.zeros((nchan,nchan), dtype=n.complex)
        FI[m]  = n.zeros((nchan,nchan), dtype=n.complex)
        for i in xrange(nchan):
            print 'Fisher row:', i
            for j in xrange(nchan):
                _CQ[i] = _CQ.get(i, n.dot(_C[m],Q[i])) # XXX ensure this has been calculated
                _CQ[j] = _CQ.get(j, n.dot(_C[m],Q[j])) # XXX ensure this has been calculated
                FC[m][i,j] = n.einsum('ij,ji', _CQ[i], _CQ[j])
                FI[m][i,j] = n.einsum('ij,ji', Q[i], Q[j])
            #bI [i,0] = n.einsum('ij,ji', E[i], N)
            #bC[i,0] = n.einsum('ij,ji', _CE[i], _CN)

        if PLOT:
            p.subplot(121); capo.arp.waterfall(FC[m], drng=4)
            p.subplot(123); capo.arp.waterfall(FI[m], drng=4)
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
        FC_o = n.take(n.take(FC[m],order, axis=0), order, axis=1)
        L_o = n.linalg.cholesky(FC_o)
        U,S,V = n.linalg.svd(L_o.conj())
        _S = 1./S
        MC_o = n.dot(n.transpose(V), n.dot(n.diag(_S), n.transpose(U)))
        MC[m] = n.take(n.take(MC_o,iorder, axis=0), iorder, axis=1)
        MI[m]  = n.identity(nchan, dtype=n.complex128)
        
        print 'Normalizing M/W'
        WI[m] = n.dot(MI[m], FI[m])
        norm  = WI[m].sum(axis=-1); norm.shape += (1,)
        MI[m] /= norm; WI[m] = n.dot(MI[m], FI[m])
        WC[m] = n.dot(MC[m], FC[m])
        norm  = WC[m].sum(axis=-1); norm.shape += (1,)
        MC[m] /= norm; WC[m] = n.dot(MC[m], FC[m])

        print 'Generating qs'
        qC = n.array([_Cx.conj() * n.dot(Q[i], _Cx) for i in xrange(nchan)])
        qC = n.sum(qC, axis=1)
        #qC -= bC # subtract noise bias
        qI = n.array([x[m].conj() * n.dot(Q[i], x[m]) for i in xrange(nchan)])
        qI = n.sum(qI, axis=1)
        #qI -= bI # subtract noise bias

        print 'Generating ps'
        pC[m] = n.dot(MC[m], qC) * scalar
        pI[m] = n.dot(MI[m], qI) * scalar

        if PLOT:
            p.subplot(411); capo.arp.waterfall(qC, mode='real'); p.colorbar(shrink=.5)
            p.subplot(412); capo.arp.waterfall(pC[m], mode='real'); p.colorbar(shrink=.5)
            p.subplot(413); capo.arp.waterfall(qI, mode='real'); p.colorbar(shrink=.5)
            p.subplot(414); capo.arp.waterfall(pI[m], mode='real'); p.colorbar(shrink=.5)
            p.show()

    if PLOT:
        p.plot(kpl, n.average(pC['s'].real, axis=1), 'b.-')
        p.plot(kpl, n.average(pC['d'].real, axis=1), 'g.-')
        p.plot(kpl, n.average((pC['s']-pC['d']).real, axis=1), 'r.-')
        p.plot(kpl, n.average(pI['s'].real, axis=1), 'k.-')
        p.plot(kpl, n.average(pI['d'].real, axis=1), 'c.-')
        p.plot(kpl, n.average((pI['s']-pI['d']).real, axis=1), 'm.-')
        p.show()

    print 'Writing pspec_boot%04d.npz' % boot
    n.savez('pspec_boot%04d.npz'%boot, kpl=kpl, scalar=scalar, times=n.array(times1),
        pk_vs_t=pC['s']-pC['d'], err_vs_t=1./cnt, temp_noise_var=var, nocov_vs_t=pI['s']-pI['d'],
        cmd=' '.join(sys.argv))


