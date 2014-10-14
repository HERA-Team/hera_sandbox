#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo
import glob, optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
#o.add_option('-b', '--boot', type='int', default=20,
#    help='Number of bootstraps.  Default is 20')
#o.add_option('--plot', action='store_true',
#    help='Generate plots')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
opts,args = o.parse_args(sys.argv[1:])

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

def noise(size):
    return n.random.normal(size=size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=size))

def same_group(bli, blj, gps):
    for gp in gps:
        if bli in gp and blj in gp: return gp
    return None

def cov_average(C, bls, n_k, gps):
    '''gps is the list of all groups in all chunks.'''
    nbls = len(bls)
    original_shape = C.shape
    C.shape = (nbls,n_k,nbls,n_k)
    Cavg = n.zeros_like(C)
    #choose a (i,j) baseline cross-multiple panel in the covariance matrix
    for i in xrange(nbls):
        bli = bls[i]
        for j in xrange(nbls):
            blj = bls[j]
            #make sure we only compute average using baselines in the same group
            gp = same_group(bli, blj, gps)
            if gp is None: continue
            #Now average over all other panels of covariance matrix (within this group)
            #to get the average signal covariance and subtract that off so that we don't 
            #get signal loss removing residual signal covariances.
            _Csum,_Cwgt = 0,0
            for i_ in xrange(nbls):
                bli_ = bls[i_]
                if not bli_ in gp: continue # make sure averaging over baseline in same group.
                if bli == bli_: continue #only average over other bls to better isolate bl systematics.
                for j_ in xrange(nbls):
                    blj_ = bls[j_]
                    if not blj_ in gp: continue #make sure averaging over baselines in same gp. 
                    if bli_ == blj_: continue #don't average panels with noise bias (same bls).
                    if blj == blj_: continue #only average over other bls to bettwe isolate bl systematics.
                    _Csum += C[i_,:,j_]
                    _Cwgt += 1
            try:
                Cavg[i,:,j] = _Csum/_Cwgt
            except(ZeroDivisionError): #catches zero weights
                print 'weights are zero for %d_%d'%(i_,j_)
                sys.stdout.flush()
                Cavg[i,:,j] = _Csum
    C.shape = Cavg.shape = original_shape
    return Cavg

def tile_panels(panel, nbls):
    n_k,n_k = panel.shape
    C = n.zeros((nbls,n_k,nbls,n_k), dtype=panel.dtype)
    for i in xrange(nbls):
        for j in xrange(nbls):
            C[i,:,j,:] = panel.copy()
    C.shape = (nbls*n_k, nbls*n_k)
    return C


def cov_average(C, bls, n_k, gps):
    nbls = len(bls)
    original_shape = C.shape
    C.shape = (nbls,n_k,nbls,n_k)
    Cavg = n.zeros_like(C)
    dsum, dwgt = 0, 0
    for i in xrange(nbls):
        for j in xrange(nbls):
            if i == j: continue
            dsum += C[i,:,j,:]
            dwgt += 1
    C.shape = original_shape
    davg = dsum / dwgt
    return tile_panels(davg, nbls)

def get_E(mode, n_k, nbls):
    # XXX need to integrate across waveband, not just mode at center point
    _m = n.zeros((n_k,), dtype=n.complex)
    _m[mode] = 1.
    m = n.fft.fft(n.fft.ifftshift(_m)) * a.dsp.gen_window(nchan, WINDOW)
    E = n.einsum('i,j', m, m.conj())
    return tile_panels(E, nbls)
    #E = tile_panels(E, nbls)
    #E.shape = (nbls,n_k,nbls,n_k)
    #for i in xrange(nbls): E[i,:,i] = 0 # XXX
    #E.shape = (nbls*n_k,nbls*n_k)
    #return E

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
B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW] # normalization. See above.
etas = n.fft.fftshift(capo.pspec.f2eta(afreqs)) #create etas (fourier dual to frequency)
#etas = capo.pspec.f2eta(afreqs) #create etas (fourier dual to frequency)
kpl = etas * capo.pspec.dk_deta(z) #111
print kpl
bm = n.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq)

if True: scalar = capo.pspec.X2Y(z) * bm * B
else: scalar = 1
print 'Freq:',fq
print 'z:', z
print 'B:', B
print 'scalar:', scalar
sys.stdout.flush()


#files = glob.glob('*.uvA')
#files = glob.glob('*.uvL')
files = args
#antstr = '41_49,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_47,3_25,1_58,35_61,42_63,2_33'
#antstr = '41_49,3_10,9_58,22_61,20_63'#,2_43,21_53,31_45,41_47,3_25,1_58,35_61,42_63,2_33'
antstr = 'cross'
t,d,f = capo.arp.get_dict_of_uv_data(files, antstr=antstr, polstr='I', verbose=True)
aa = a.cal.get_aa(opts.cal, n.array([.150]))
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
window = a.dsp.gen_window(nchan, WINDOW)
if not WINDOW == 'none': window.shape=(1,nchan)
_d = {}
for k in d:
    d[k]['I'] = d[k]['I'][:,chans]
    #d[k]['I'] = d[k]['I'][500:1000,chans]
    if conj[k]: d[k]['I'] = n.conj(d[k]['I'])
    #_d[k] = n.fft.fftshift(n.fft.ifft(window*d[k]['I']), axes=1)
    _d[k] = d[k]['I']
    _d[k] = n.transpose(_d[k], [1,0])
bls = _d.keys()
print len(bls)

if True: _data = n.array([_d[k] for k in bls])
else:
    print 'OVERRIDING WITH SIMULATED SIGNAL'
    eor = noise(_d[bls[0]].shape)
    _data = n.array([eor for k in bls])
_data = n.reshape(_data, (_data.shape[0]*_data.shape[1], _data.shape[2]))

if False:
    print 'ADDING NOISE'
    _data += 1.0 * noise(_data.shape)

_data *= n.sqrt(scalar)

#x1,x2 = _d[bls[0]], _d[bls[1]]
#C12 = cov2(x1,x2)
#C21 = cov2(x2,x1)
#C = 0.5 * (C12 + C21)
if True:
    C = cov(_data)
else:
    'OVERRIDING COVARIANCE MATRIX'
    C = tile_panels(n.identity(nchan), len(bls))
    C += n.identity(C.shape[0])
#C *= tile_panels(n.identity(nchan), len(bls)) # XXX
#C *= n.identity(C.shape[0]) # XXX
Cavg = cov_average(C, bls, nchan, [bls])
N = C - Cavg
C = C * n.identity(C.shape[0]) + Cavg * tile_panels(n.identity(nchan), len(bls))
#C = C * n.identity(C.shape[0]) + Cavg
if False:
    _C = n.linalg.inv(C)
else:
    print 'Psuedoinverse of C'
    U,S,V = n.linalg.svd(C.conj())
    print S
    #_S = n.where(S > .1, 1./S, 0)
    _S = n.where(S > 1e7, 1./S, 0)
    _C = n.einsum('ij,j,jk', V.T, _S, U.T)
p.subplot(221); capo.arp.waterfall(C, drng=3)
p.subplot(222); capo.arp.waterfall(_C, drng=3)
p.subplot(223); capo.arp.waterfall(n.dot(C,_C), drng=3)
p.subplot(224); capo.arp.waterfall(N, drng=3)
p.show()
#_Cx1 = n.dot(_C, x1)
#_Cx2 = n.dot(_C, x2)
_Cx = n.dot(_C, _data)

p.subplot(211); capo.arp.waterfall(_data, mode='real'); p.colorbar(shrink=.5)
p.subplot(212); capo.arp.waterfall(  _Cx, mode='real'); p.colorbar(shrink=.5)
p.show()

#def get_E(mode):
#    E = n.zeros_like(C)
#    E[mode,mode] = 1
#    return E


E, _CE = {}, {}
bC = n.zeros((nchan,1), dtype=n.complex)
b  = n.zeros((nchan,1), dtype=n.complex)
_CN = n.dot(_C,N)
for i in xrange(nchan):
    print 'E',i
    #E[i] = get_E(i)
    E[i] = get_E(i, nchan, len(bls))
    _CE[i] = n.dot(_C,E[i])
    b [i,0] = n.einsum('ij,ji', E[i], N)
    bC[i,0] = n.einsum('ij,ji', _CE[i], _CN)
    #capo.arp.waterfall(E[i], mode='real')
    #p.show()
FC = n.zeros((nchan,nchan), dtype=n.complex)
F  = n.zeros((nchan,nchan), dtype=n.complex)
for i in xrange(nchan):
    print i
    for j in xrange(nchan):
        FC[i,j] = n.einsum('ij,ji', _CE[i], _CE[j])
        F[i,j] = n.einsum('ij,ji', E[i], E[j])

p.subplot(121); capo.arp.waterfall(FC, drng=4)
p.subplot(122); capo.arp.waterfall(F , drng=4)
p.show()

if False:
    MC = n.identity(nchan, dtype=n.complex128)
    #MC = n.linalg.inv(FC)
else:
    print 'Psuedoinverse of FC'
    U,S,V = n.linalg.svd(FC.conj())
    print S
    #_S = n.where(S > 100, 1./S, 0)
    #_S = n.where(S > 1e-15, n.sqrt(1./S), 0)
    _S = n.sqrt(1./S)
    #_S = 1./S
    MC = n.dot(n.transpose(V), n.dot(n.diag(_S), n.transpose(U)))
M  = n.identity(nchan, dtype=n.complex128)
#M2 = n.linalg.inv(F2)
WC = n.dot(MC, FC)
W  = n.dot(M , F )
norm1 = WC.sum(axis=-1); norm1.shape += (1,)
MC /= norm1
norm  = W.sum(axis=-1); norm.shape += (1,)
M /= norm
qCa = n.array([_Cx.conj() * n.dot(E[i], _Cx) for i in xrange(nchan)])
qCa = n.sum(qCa, axis=1)
print qCa.shape, bC.shape
qCa -= bC
pCa = n.dot(MC, qCa)
qa = n.array([_data.conj() * n.dot(E[i], _data) for i in xrange(nchan)])
qa = n.sum(qa, axis=1)
qa -= b
pa = n.dot(M, qa)

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


'''
C = cov(_data)
Cavg = cov_average(C, bls, nchan, [bls])
N = C - Cavg
if False: _N = n.linalg.inv(N)
else:
    U,S,V = n.linalg.svd(N.conj())
    print 'Psuedoinverse of N'
    print S
    _S = n.where(S > 1e-1, 1./S, 0)
    _N = n.dot(n.transpose(V), n.dot(n.diag(_S), n.transpose(U)))
if True: _C = n.linalg.inv(C)
else:
    U,S,V = n.linalg.svd(C.conj())
    print 'Psuedoinverse of C'
    print S
    _S = n.where(S > 1e0, 1./S, 0)
    _C = n.dot(n.transpose(V), n.dot(n.diag(_S), n.transpose(U)))

_Cd = n.dot(_C, _data)
_Nd = n.dot(_N, _data)

def get_E(mode):
    E = n.zeros_like(C)
    ind = n.arange(len(bls)) * nchan
    for i in ind: E[i+mode,ind+mode] = 1
    ind = n.arange(E.shape[0])
    #E[ind,ind] = 0
    return E

E = {}
_CE = {}
_NE = {}
for i in xrange(nchan):
    E[i] = get_E(i)
    _CE[i] = n.dot(_C,E[i])
    _NE[i] = n.dot(_N,E[i])
F1 = n.zeros((nchan,nchan), dtype=n.complex)
F2 = n.zeros((nchan,nchan), dtype=n.complex)
F3 = n.zeros((nchan,nchan), dtype=n.complex)
for i in xrange(nchan):
    print i
    for j in xrange(nchan):
        F1[i,j] = n.trace(n.dot(_NE[i], _NE[j]))
        F2[i,j] = n.trace(n.dot(_CE[i], _CE[j]))
        F3[i,j] = n.trace(n.dot(E[i], E[j]))

M1 = n.identity(nchan, dtype=n.complex128)
M2 = n.identity(nchan, dtype=n.complex128)
M3 = n.identity(nchan, dtype=n.complex128)
M1 = n.linalg.inv(F1)
M2 = n.linalg.inv(F2)
M3 = n.linalg.inv(F3)
W1 = n.dot(M1, F1)
W2 = n.dot(M2, F2)
W3 = n.dot(M3, F3)
norm1 = W1.sum(axis=-1); norm1.shape += (1,)
M1 /= norm1
norm2 = W2.sum(axis=-1); norm2.shape += (1,)
M2 /= norm2
norm3 = W3.sum(axis=-1); norm3.shape += (1,)
M3 /= norm3

#qa = n.array([n.dot(_d.T.conj(), n.dot(E[i], _d)) for i in xrange(nchan)])
#q1a = n.array([_d.conj() * n.dot(E[i], _d) for i in xrange(nchan)])
q1a = n.array([_Nd.conj() * n.dot(E[i], _Nd) for i in xrange(nchan)])
q1a = n.sum(q1a, axis=1)
p1a = n.dot(M1, q1a)
q2a = n.array([_Cd.conj() * n.dot(E[i], _Cd) for i in xrange(nchan)])
q2a = n.sum(q2a, axis=1)
p2a = n.dot(M2, q2a)
q3a = n.array([_data.conj() * n.dot(E[i], _data) for i in xrange(nchan)])
q3a = n.sum(q3a, axis=1)
p3a = n.dot(M3, q3a)


#C2 = cov(_d)
'''


