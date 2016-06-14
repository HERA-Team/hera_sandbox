#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo
import capo.frf_conv as fringe
import glob, optparse, sys, random
import capo.zsa as zsa
import scipy

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
o.add_option('-b', '--nboot', type='int', default=20,
    help='Number of bootstraps.  Default is 20')
o.add_option('-i', '--inject_sig', type='float', default=1.,
    help='Inject signal amplitude.  Default is 1.')
o.add_option('--plot', action='store_true',
    help='Generate plots')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
o.add_option('--bl_scale', type='float',default=1.0,
    help='Scale to change baseline length used to make Fringe Rateswhile filtering eor')
o.add_option('--fr_width_scale', type='float',default=1.0,
    help='Artificially inflate the width of the FRP by scaling factor')
o.add_option('--output', type='string', default='',
    help='output directory for pspec_boot files (default "")')
opts,args = o.parse_args(sys.argv[1:])

#Basic Parameters
random.seed(0)
#n.random.seed(1235813)
POL = opts.pol #'I'
DELAY = False
NGPS = 5
INJECT_SIG = opts.inject_sig
FRF_WIDTH = 401
NOISE = 13
ntimes = 609
PLOT = opts.plot
days = ['even','odd']
SEP='0,1'

## Define some Functions

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

def noise(size):
    size = list(size)
    size[-1] *= 10
    size = tuple(size)
    return n.random.normal(size=size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=size))

def gen_flags(size,f):
    size = list(size)
    size[-1] *= 10
    f_size = list(f.shape)
    start= int( ( size[-1] - f_size[-1] )/2.)
    end = int( ( size[-1] + f_size[-1])/2. )
    diff = end-start - f_size[-1]
    if diff != 0: end += diff
    ind = n.arange(start, end,1)
    flags = n.zeros(size)
    for ch in xrange(size[0]): flags[ch,ind] = n.copy(f[ch])
    return n.transpose(flags, [1,0])

def clip_array(d,size,axis=0):
    if size is None:
        print 'Must give valid size for clipping'
        print 'Returning origian larray'
        return d
    d_size = d.shape[axis]
    if d_size < size:
        print 'Array Must have large size than desired clip'
        print 'Clipping failed, returning original array'
        return d
    start= int( (d_size - size)/2.)
    end = int( (d_size + size)/2. )
    diff = size - ( end - start )
    if diff != 0: end += diff
    d = n.swapaxes(d,axis,0)
    _d = d[start:end]
    _d = n.swapaxes(_d,0, axis)
    return _d

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

#Get uv file info
WINDOW = opts.window
freqs = n.linspace(0.1,0.2,num=203)
sdf = n.diff(freqs)[0]
chans = a.scripting.parse_chans(opts.chan, 203)
inttime=42.9499
print 'inttime', inttime

afreqs = freqs.take(chans)
nchan = chans.size
times= n.arange(0,ntimes)*inttime
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

#aa = a.cal.get_aa(opts.cal, n.array([.150]))
aa = a.cal.get_aa(opts.cal, freqs)
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
jy2T = capo.pspec.jy2T(afreqs)
window = a.dsp.gen_window(nchan, WINDOW)
if not WINDOW == 'none': window.shape=(nchan,1)
sep2ij, blconj, bl2sep = zsa.grid2ij(aa.ant_layout)

ijs = sep2ij[SEP].split(',')
#all_bls= [a.miriad.ij2bl(  *map(int, x.split('_'))) for x in ijs]
all_bls= [ a.miriad.ij2bl(*map( int,x.split('_'))) for x in ijs]

# XXX NEED TO FIGURE OUT BW NORMALIZATION
B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW] #proper normalization
etas = n.fft.fftshift(capo.pspec.f2eta(afreqs)) #create etas (fourier dual to frequency)
#etas = capo.pspec.f2eta(afreqs) #create etas (fourier dual to frequency)
kpl = etas * capo.pspec.dk_deta(z) #111
#print kpl
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


if True: #this one is the exact one
    sep = bl2sep[all_bls[0]]
    ij_array =  sep2ij[sep].split(',')
    while True:
        ij = map( int, ij_array.pop().split('_') )
        bl = a.miriad.ij2bl(*ij )
        if not blconj[bl]: break
    print 'Using Baseline for FRP:',bl
    bins = fringe.gen_frbins(inttime)
    frp, bins = fringe.aa_to_fr_profile(aa, ij, len(afreqs)/2, bins=bins, bl_scale = opts.bl_scale)

    timebins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(), fq0=aa.get_freqs()[len(afreqs)/2], bl_scale=opts.bl_scale, fr_width_scale= opt.fr_width_scale)

    if blconj[a.miriad.ij2bl(ij[0],ij[1])]: fir = {(ij[0],ij[1],POL):n.conj(firs)} #conjugate fir if needed
    else: fir = {(ij[0],ij[1],POL):firs}

#Extract frequency range of data for each boot
for boot in xrange(opts.nboot):
    print '%d / %d' % (boot+1,opts.nboot)
    x = {}
    f = {}
    for k in days:
        x[k] = {}
        f[k] = {}
        for bl in all_bls:
            d = noise((nchan,ntimes)).T * NOISE * jy2T
            flg = n.zeros_like(d)
            if conj[bl]: d = n.conj(d) #conjugate if necessary
            x[k][bl] = n.transpose(d, [1,0]) #swap time and freq axes
            f[k][bl] = n.transpose(flg, [1,0])
    #eor = x.pop('eor'); days = x.keys() #make up for putting eor in list above
    bls_master = x.values()[0].keys()
    nbls = len(bls_master)
    print 'Baselines:', nbls
    if INJECT_SIG > 0.: #Create a fake EoR signal to inject
        print 'INJECTING SIMULATED SIGNAL'
        eor1 = noise((nchan,ntimes)).T * INJECT_SIG * jy2T
        wij_ = n.zeros_like(eor1)
        dij,wij = eor1 , n.logical_not(wij_)
        _d,_w,_,_ = fringe.apply_frf(aa,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)

            ### OLD CODE TO FRF ###
            #for cnt,ch in enumerate(chans):
            #    eor1[cnt] = n.convolve(eor1[cnt], n.conj(firs[cnt]), mode='same') #conjugate firs!!!
        eor2 = n.transpose(_d, [1,0])
    #    eor2 = clip_array(eor2,x[k][bls_master[0]].shape)
        eor = n.copy(eor2)
     #   eor1 = clip_array(eor1,x[k][bls_master[0]].shape)
        for k in days:
            for bl in x[k]:
                _x = x[k][bl].copy() + eor.copy() #add injected signal to data
                wij = n.logical_not(f[k][bl].T)
                _x_,_,_,_ =fringe.apply_frf(aa,_x.T,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)
                _x_ = n.transpose(_x_,[1,0])
                x[k][bl] = clip_array(_x_,ntimes,axis=1)

        eor = clip_array(eor,ntimes,axis=1)
    #Power spectrum stuff
    #Q = {} # Create the Q's that extract power spectrum modes
    #for i in xrange(nchan):
    #    Q[i] = get_Q(i, nchan)
    Q = [get_Q(i,nchan) for i in xrange(nchan)] #get Q matrix (does FT from freq to delay)

    #Compute baseline auto-covariances and apply inverse to data
    I,_I,_Ix = {},{},{}
    C,_C,_Cx = {},{},{}
    for k in days:
        I[k],_I[k],_Ix[k] = {},{},{}
        C[k],_C[k],_Cx[k] = {},{},{}
        for bl in x[k]:
            C[k][bl] = cov(x[k][bl])
            I[k][bl] = n.identity(C[k][bl].shape[0])
            #C[k][bl] = C[k][bl] + 1*I[k][bl] #C+IN noise
            U,S,V = n.linalg.svd(C[k][bl].conj()) #singular value decomposition
            _C[k][bl] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
            _I[k][bl] = n.identity(_C[k][bl].shape[0])
            #_Cx[k][bl] = n.dot(_C[k][bl], x[k][bl])
            #_Ix[k][bl] = x[k][bl].copy()
            _Cx[k][bl] = n.dot(_C[k][bl], eor) # XXX
            _Ix[k][bl] = eor.copy() # XXX
            if PLOT and True:
                #p.plot(S); p.show()
                print a.miriad.bl2ij(bl), k
                p.suptitle('{0}_{1}'.format(a.miriad.bl2ij(bl)[0],a.miriad.bl2ij(bl)[1]))
                p.subplot(311); capo.arp.waterfall(x[k][bl], mode='real')
                p.subplot(334); capo.arp.waterfall(C[k][bl])
                p.subplot(335); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
                p.subplot(336); capo.arp.waterfall(_C[k][bl])
                p.subplot(313); capo.arp.waterfall(_Cx[k][bl], mode='real')
                p.show()
        
    #Make boots
    bls = bls_master[:]
    if True: #shuffle and group baselines for bootstrapping
        random.shuffle(bls)
        #bls = bls[:-5] # XXX
        gps = [bls[i::NGPS] for i in range(NGPS)]
        gps = [[random.choice(gp) for bl in gp] for gp in gps]
    else: gps = [bls[i::NGPS] for i in range(NGPS)]
    bls = [bl for gp in gps for bl in gp]
    #print '\n'.join([','.join(['%d_%d'%a.miriad.bl2ij(bl) for bl in gp]) for gp in gps])
    _Iz,_Isum,_IsumQ = {},{},{}
    _Cz,_Csum,_CsumQ = {},{},{}
    print "   Getting C"
    for k in days:
        _Iz[k],_Isum[k],_IsumQ[k] = {},{},{}
        _Cz[k],_Csum[k],_CsumQ[k] = {},{},{}
        for i,gp in enumerate(gps): #sum things up over the groups
            _Iz[k][i] = sum([_Ix[k][bl] for bl in gp])
            _Cz[k][i] = sum([_Cx[k][bl] for bl in gp])
            _Isum[k][i] = sum([_I[k][bl] for bl in gp])
            _Csum[k][i] = sum([_C[k][bl] for bl in gp])
            _IsumQ[k][i] = {}
            _CsumQ[k][i] = {}
            if DELAY: #this is much faster
                _Iz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Iz[k][i], axis=0), axes=0)
                _Cz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Cz[k][i], axis=0), axes=0)
                # XXX need to take fft of _Csum, _Isum here
            for ch in xrange(nchan): # XXX this loop makes computation go as nchan^3
                _IsumQ[k][i][ch] = n.dot(_Isum[k][i], Q[ch])
                _CsumQ[k][i][ch] = n.dot(_Csum[k][i], Q[ch])
        if PLOT:
            NGPS = len(gps)
            _Csumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            _Isumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            for i in xrange(len(gps)): _Isumk[i,:,i,:] = _Isum[k][i]
            _Isumk.shape = (NGPS*nchan, NGPS*nchan)
            #_Isum[k] = _Isumk
            for i in xrange(len(gps)): _Csumk[i,:,i,:] = _Csum[k][i]
            _Csumk.shape = (NGPS*nchan, NGPS*nchan)
            #_Csum[k] = _Csumk
            _Czk = n.array([_Cz[k][i] for i in _Cz[k]])
            print _Czk.shape
            _Czk = n.reshape(_Czk, (_Czk.shape[0]*_Czk.shape[1], _Czk.shape[2]))
            p.subplot(211); capo.arp.waterfall(_Czk, mode='real')
            p.subplot(223); capo.arp.waterfall(_Csumk)
            p.subplot(224); capo.arp.waterfall(cov(_Czk))
            p.show()
    print "   Getting F and q"
    FI = n.zeros((nchan,nchan), dtype=n.complex)
    FC = n.zeros((nchan,nchan), dtype=n.complex)
    qI = n.zeros((nchan,_Iz.values()[0].values()[0].shape[1]), dtype=n.complex)
    qC = n.zeros((nchan,_Cz.values()[0].values()[0].shape[1]), dtype=n.complex)
    Q_Iz = {}
    Q_Cz = {}
    for cnt1,k1 in enumerate(days):
        for k2 in days[cnt1:]: #loop over even with eve, even with odd, etc.
            if not Q_Iz.has_key(k2): Q_Iz[k2] = {}
            if not Q_Cz.has_key(k2): Q_Cz[k2] = {}
            for bl1 in _Cz[k1]:
                for bl2 in _Cz[k2]:
                    if k1 == k2 or bl1 == bl2: continue
                    #if k1 == k2 and bl1 == bl2: continue # this results in a significant bias
                    #if bl1 == bl2: continue # also a significant noise bias
                    #print k1, k2, bl1, bl2
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
                    if False: #use ffts to do q estimation fast
                        qI += n.conj(_Iz[k1][bl1]) * _Iz[k2][bl2]
                        qC += n.conj(_Cz[k1][bl1]) * _Cz[k2][bl2]
                    else: #brute force with Q to ensure normalization
                        #_qI = n.array([_Iz[k1][bl1].conj() * n.dot(Q[i], _Iz[k2][bl2]) for i in xrange(nchan)])
                        #_qC = n.array([_Cz[k1][bl1].conj() * n.dot(Q[i], _Cz[k2][bl2]) for i in xrange(nchan)])
                        if not Q_Iz[k2].has_key(bl2): Q_Iz[k2][bl2] = [n.dot(Q[i], _Iz[k2][bl2]) for i in xrange(nchan)]
                        if not Q_Cz[k2].has_key(bl2): Q_Cz[k2][bl2] = [n.dot(Q[i], _Cz[k2][bl2]) for i in xrange(nchan)]
                        _qI = n.array([_Iz[k1][bl1].conj() * Q_Iz[k2][bl2][i] for i in xrange(nchan)])
                        qI += n.sum(_qI, axis=1)
                        _qC = n.array([_Cz[k1][bl1].conj() * Q_Cz[k2][bl2][i] for i in xrange(nchan)]) #C^-1 Q C^-1
                        qC += n.sum(_qC, axis=1)
                        #qC += [ n.diag(n.dot(_Cz[k1][bl1].T.conj() , Q_Cz[k2][bl2][i])) for i in xrange(nchan)] #C^-1 Q C^-1
                        #qC += n.sum(_qC, axis=1)
                    if DELAY: #by taking FFT of CsumQ above, each channel is already i,j separated
                        FI += n.conj(_IsumQ[k1][bl1]) * _IsumQ[k2][bl2]
                        FC += n.conj(_CsumQ[k1][bl1]) * _CsumQ[k2][bl2]
                    else:
                        for i in xrange(nchan):
                            for j in xrange(nchan):
                                FI[i,j] += n.einsum('ij,ji', _IsumQ[k1][bl1][i], _IsumQ[k2][bl2][j])
                                FC[i,j] += n.einsum('ij,ji', _CsumQ[k1][bl1][i], _CsumQ[k2][bl2][j]) #C^-1 Q C^-1 Q

    if PLOT:
        p.subplot(121); capo.arp.waterfall(FC, drng=4)
        p.subplot(122); capo.arp.waterfall(FI, drng=4)
        p.show()

    #print 'Psuedoinverse of FC'
    
    #Other choices for M
    #U,S,V = n.linalg.svd(FC.conj())
    #_S = n.sqrt(1./S)
    # _S = 1./S
    # _S = n.ones_like(S)
    #MC = n.dot(n.transpose(V), n.dot(n.diag(_S), n.transpose(U)))
    #order = n.array([10,11,9,12,8,13,7,14,6,15,5,16,4,17,3,18,2,19,1,20,0])

    print "   Getting M"
    #Cholesky decomposition
    order = n.array([10,11,9,12,8,20,0,13,7,14,6,15,5,16,4,17,3,18,2,19,1])
    iorder = n.argsort(order)
    FC_o = n.take(n.take(FC,order, axis=0), order, axis=1)
    L_o = n.linalg.cholesky(FC_o)
    #_,L_o,U = scipy.linalg.lu(FC_o)
    U,S,V = n.linalg.svd(L_o.conj())
    MC_o = n.dot(n.transpose(V), n.dot(n.diag(1./S), n.transpose(U)))
    MC = n.take(n.take(MC_o,iorder, axis=0), iorder, axis=1)
    MI  = n.identity(nchan, dtype=n.complex128)
    
    print "   Getting W"
    #print 'Normalizing M/W'
    WI = n.dot(MI, FI)
    norm  = WI.sum(axis=-1); norm.shape += (1,)
    #norm  = WI.max(axis=-1); norm.shape += (1,) # XXX
    MI /= norm; WI = n.dot(MI, FI)
    WC = n.dot(MC, FC)
    norm  = WC.sum(axis=-1); norm.shape += (1,)
    #norm  = WC.max(axis=-1); norm.shape += (1,) # XXX
    MC /= norm; WC = n.dot(MC, FC)

    print '   Generating ps'
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

    print 'pI=', n.average(pI.real), 'pC=', n.average(pC.real), 'pI/pC=', n.average(pI.real)/n.average(pC.real)

    if PLOT:
        p.plot(kpl, n.average(pC.real, axis=1), 'b.-')
        p.plot(kpl, n.average(pI.real, axis=1), 'k.-')
        p.show()

    print '   Writing pspec_bootsigloss%04d.npz' % boot

    if len(opts.output) > 0: outpath = opts.output+'/pspec_boot%04d.npz' % boot
    else: outpath = 'pspec_boot%04d.npz' % boot

    n.savez(outpath, kpl=kpl, scalar=scalar, times=times,
        pk_vs_t=pC, nocov_vs_t=pI, cmd=' '.join(sys.argv))


