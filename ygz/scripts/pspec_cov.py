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

MASK = True
DELAY = False
CHOLESKY = True
#NGPS = 4
NGPS = 2
AVG_COV = True
LST_STATS = True

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

def same_group(bli, blj, gps):
    for gp in gps:
        if bli in gp and blj in gp: return gp
    return None

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

#def ndiag(size, num):
#    C = n.identity(size)
#    ind = n.arange(size)
#    for i in xrange(1,num):
#        for j in xrange(size):
#            C[j,(j+i)%size] = 1
#            C[(j+i)%size,j] = 1
#    return C

def cov_average(C, bls, nchan):
    '''Return bl levels, autos, crosses'''
    nbls = len(bls)
    original_shape = C.shape
    C.shape = (nbls,nchan,nbls,nchan)
    dsum, dwgt = {}, {}
    level = dict(zip(bls, [n.trace(C[i,:,i,:])/nchan for i in range(nbls)]))
    for i in xrange(nbls):
        for j in xrange(nbls):
            dsum[i == j] = dsum.get(i == j, 0) + C[i,:,j,:]
            dwgt[i == j] = dwgt.get(i == j, 0) + 1
    C.shape = original_shape
    return level, dsum[True]/dwgt[True], dsum[False]/dwgt[False]

def get_E(mode, n_k, bls, zero_block_diag=False, zero_repeat_bls=False):
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
#antstr = '41_49,3_10,9_58,22_61,20_63'#,2_43,21_53,31_45,41_47,3_25,1_58,35_61,42_63,2_33'
antstr = 'cross'
times,data,flgs = capo.arp.get_dict_of_uv_data(args, antstr=antstr, polstr='I', verbose=True)

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
_d = {}
for k in data:
    data[k]['I'] = data[k]['I'][:,chans] * jy2T
    if conj[k]: data[k]['I'] = n.conj(data[k]['I'])
    if DELAY: _d[k] = n.fft.fftshift(n.fft.ifft(window*data[k]['I']), axes=1)
    else: _d[k] = data[k]['I']
    _d[k] = n.transpose(_d[k], [1,0])
bls_master = _d.keys()
nbls = len(bls_master)
print nbls

_data = n.array([_d[k] for k in bls_master])
_data = n.reshape(_data, (_data.shape[0]*_data.shape[1], _data.shape[2]))
C = cov(_data)
level,auto,cross = cov_average(C, bls_master, nchan)

for boot in xrange(opts.nboot):
    print '%d / %d' % (boot+1,opts.nboot)
    bls = bls_master[:]
    random.shuffle(bls)
    #bls = bls[:-10] # XXX
    bls = bls[:-5] # XXX
    nbls = len(bls)
    gps = [bls[i::NGPS] for i in range(NGPS)]
    #gps = [[random.choice(gp) for bl in gp] for gp in gps]
    gps = [[bl for bl in gp] for gp in gps]
    bls = [bl for gp in gps for bl in gp]
    print '\n'.join([','.join(['%d_%d'%a.miriad.bl2ij(bl) for bl in gp]) for gp in gps])

    if True: _data = n.array([_d[k] for k in bls])
    else:
        print 'OVERRIDING WITH SIMULATED SIGNAL'
        eor = noise(_d[bls[0]].shape)
        _data = n.array([eor for k in bls])
    _data = n.reshape(_data, (_data.shape[0]*_data.shape[1], _data.shape[2]))
    print _data.shape

    if False: # add in additional noise
        print 'ADDING NOISE'
        _data += 1.0 * noise(_data.shape)

    C = cov(_data)
    #N = C - cross

    if MASK: # mask covariance matrix
        mask = n.ones((nbls,nchan,nbls,nchan))
        for i in n.cumsum([len(gp) for gp in gps])[:-1]:
            mask[:i,:,i:,:] = 0
            mask[i:,:,:i,:] = 0
        mask.shape = (nbls*nchan,nbls*nchan)
        #C = C * n.identity(C.shape[0])
        #C = C * n.identity(C.shape[0]) + Cavg * tile_panels(n.identity(nchan), len(bls))
        #C += Cavg * tile_panels(n.identity(nchan), len(bls))
        #C *= tile_panels(n.identity(nchan), len(bls))
        #C = C * ndiag(C.shape[0], 2) + Cavg * tile_panels(ndiag(nchan,2), len(bls))
        #C = C * n.identity(C.shape[0]) + Cavg
    else: mask = 1

    U,S,V = n.linalg.svd(C.conj())
    if opts.plot:
        p.subplot(224); p.semilogy(S)

    if AVG_COV:
        C.shape = (nbls,nchan,nbls,nchan)
        for i in xrange(nbls):
          for j in xrange(nbls):
            if bls[i] == bls[j]: C[i,:,j,:] = level[bls[i]] * auto
            else: C[i,:,j,:] = n.sqrt(level[bls[i]] * level[bls[j]]) * cross
        C.shape = (nbls*nchan,nbls*nchan)

    C *= mask
    print 'Psuedoinverse of C'
    U,S,V = n.linalg.svd(C.conj())
    print S
    _S = n.where(S > 1e-4, 1./S, 0) # for fringe rate filtered noise.
    #_S = n.where(S > 1e-3, 1./S, 0) #original ARP threshold
    #_S = 1./S
    #_S = n.concatenate([1./S[:100], n.zeros_like(S[100:])])
    _C = n.einsum('ij,j,jk', V.T, _S, U.T)
    _Cx = n.dot(_C, _data)
    if opts.plot:
        p.subplot(221); capo.arp.waterfall(C, drng=3)
        p.subplot(222); capo.arp.waterfall(_C, drng=3)
        p.subplot(223); capo.arp.waterfall(n.dot(C,_C), drng=3)
        p.subplot(224); p.semilogy(S)
        p.show()
        p.subplot(211); capo.arp.waterfall(_data, mode='real'); p.colorbar(shrink=.5)
        p.subplot(212); capo.arp.waterfall(  _Cx, mode='real'); p.colorbar(shrink=.5)
        p.show()

    E, _CE = {}, {}
    #bC = n.zeros((nchan,1), dtype=n.complex)
    #b  = n.zeros((nchan,1), dtype=n.complex)
    #_CN = n.dot(_C,N)
    for i in xrange(nchan):
        print 'E',i
        #E[i] = get_E(i, nchan, bls, zero_repeat_bls=NOAUTOS) * (1 - mask)
        E[i] = get_E(i, nchan, bls) * (1 - mask)
        _CE[i] = n.dot(_C,E[i])
        #b [i,0] = n.einsum('ij,ji', E[i], N)
        #bC[i,0] = n.einsum('ij,ji', _CE[i], _CN)
        #capo.arp.waterfall(E[i], mode='real'); p.show()
    FC = n.zeros((nchan,nchan), dtype=n.complex)
    F  = n.zeros((nchan,nchan), dtype=n.complex)
    for i in xrange(nchan):
        print i
        for j in xrange(nchan):
            FC[i,j] = n.einsum('ij,ji', _CE[i], _CE[j])
            F[i,j] = n.einsum('ij,ji', E[i], E[j])

    if opts.plot:
        p.subplot(121); capo.arp.waterfall(FC, drng=4)
        p.subplot(122); capo.arp.waterfall(F , drng=4)
        p.show()

    if False:
        MC = n.identity(nchan, dtype=n.complex128)
        #MC = n.linalg.inv(FC)
    elif False:
        print 'Psuedoinverse of FC'
        U,S,V = n.linalg.svd(FC.conj())
        print S
        #_S = n.where(S > 100, 1./S, 0)
        #_S = n.where(S > 1e-16, n.sqrt(1./S), 0)
        _S = n.sqrt(1./S)
        #_S = 1./S
        MC = n.dot(n.transpose(V), n.dot(n.diag(_S), n.transpose(U)))
    else:
        print 'Psuedoinverse of FC'
        #order = n.array([10,11,9,12,8,13,7,14,6,15,5,16,4,17,3,18,2,19,1,20,0])
        order = n.array([10,11,9,12,8,20,0,13,7,14,6,15,5,16,4,17,3,18,2,19,1])
        iorder = n.argsort(order)
        FC_o = n.take(FC,order, axis=0)
        FC_o = n.take(FC_o,order, axis=1)
        if CHOLESKY:
            L_o = n.linalg.cholesky(FC_o)
            U,S,V = n.linalg.svd(L_o.conj())
            print S
            #_S = n.where(S > 1, 1./S, 0)
            _S = 1./S
            MC_o = n.dot(n.transpose(V), n.dot(n.diag(_S), n.transpose(U)))
        else:
            U,S,V = n.linalg.svd(FC_o.conj())
            print S
            _S = n.where(S > 2e-25, n.sqrt(1./S), 0)
            MC_o = n.dot(n.transpose(V), n.dot(n.diag(_S), n.transpose(U)))
        MC = n.take(MC_o,iorder, axis=0)
        MC = n.take(MC,iorder, axis=1)
        
    print 'Normalizing Ms and Ws'
    M  = n.identity(nchan, dtype=n.complex128)
    W  = n.dot(M , F )
    norm  = W.sum(axis=-1); norm.shape += (1,); M /= norm; W  = n.dot(M , F )
    WC = n.dot(MC, FC)
    normC = WC.sum(axis=-1); normC.shape += (1,); MC /= normC; WC = n.dot(MC, FC)
    print 'Generating qs'
    qCa = n.array([_Cx.conj() * n.dot(E[i], _Cx) for i in xrange(nchan)])
    qCa = n.sum(qCa, axis=1)
    #qCa -= bC # subtract noise bias
    qa = n.array([_data.conj() * n.dot(E[i], _data) for i in xrange(nchan)])
    qa = n.sum(qa, axis=1)
    #qa -= b # subtract noise bias
    print 'Generating ps'
    pCa = n.dot(MC, qCa) * scalar
    pa = n.dot(M, qa) * scalar

    if opts.plot:
        p.subplot(411); capo.arp.waterfall(qCa, mode='real'); p.colorbar(shrink=.5)
        p.subplot(412); capo.arp.waterfall(pCa, mode='real'); p.colorbar(shrink=.5)
        p.subplot(413); capo.arp.waterfall(qa , mode='real'); p.colorbar(shrink=.5)
        p.subplot(414); capo.arp.waterfall(pa , mode='real'); p.colorbar(shrink=.5)
        p.show()
        
        p.plot(kpl, n.average(pCa.real, axis=1), 'b.-')
        #p.plot(kpl, n.dot(MC,bC)[:,0], 'b:')
        p.plot(kpl, n.average( pa.real, axis=1), 'k.-')
        #p.plot(kpl, n.dot(M,b)[:,0], 'k:')
        p.show()
    print 'Writing pspec_boot%04d.npz' % boot
    n.savez('pspec_boot%04d.npz'%boot, kpl=kpl, scalar=scalar, times=n.array(times),
        pk_vs_t=pCa, err_vs_t=1./cnt, temp_noise_var=var, nocov_vs_t=pa,
        cmd=' '.join(sys.argv))


