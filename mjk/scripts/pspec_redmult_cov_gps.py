#! /usr/bin/env python
import aipy as a, numpy as n
import capo
import optparse, sys, os, random

def miriadbl2str(mbl):
    return "%d_%d"%a.miriad.bl2ij(mbl)


o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
o.add_option('-t', '--taps', type='int', default=1,
    help='Taps to use in the PFB.  Default 1, which instead uses windowed FFT')
o.add_option('-b', '--boot', type='int', default=20,
    help='Number of bootstraps.  Default is 20')
o.add_option('--plot', action='store_true',
    help='Generate plots')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
o.add_option('--gain', type='float', default=.3,
    help='gain parameter in approximation. .3 for 32, .1 for 64')
o.add_option('--usebls', action='store_true',
    help='use the baselines give in the command line. Default is use all of the given separations.')
o.add_option('--output', type='string', default='',
    help='output directory for pspec_boot files (default "")')
o.add_option('--noproj', action='store_true', 
    help='Skip the projecting out of modes inside the horizion.\
          This can lead to removign significant sky signal. Therfore \
          turning it off does not remove sky.')
o.add_option('--niters', type='string', default='', 
    help='tuple for number of steps in covariance removal')
opts,args = o.parse_args(sys.argv[1:])


PLOT = opts.plot
if PLOT: import pylab as p

NBOOT = opts.boot
NTAPS = opts.taps
if NTAPS > 1: PFB = True
else: PFB = False
WINDOW = opts.window
aa = a.cal.get_aa(opts.cal, .1, .1, 1) #bogus params in get_aa
ANTPOS = aa.ant_layout
del(aa)

# XXX Currently hardcoded for PSA898
#A_ = [0,16,8,24,4,20,12,28]
#B_ = [i+1 for i in A_]
#C_ = [i+2 for i in A_]
#D_ = [i+3 for i in A_]
#ANTPOS = n.array([A_, B_, C_, D_])

#Take half of the 64 antenna array. For testing purposes.
#1/13/14 : take the whole array.
#ANTPOS = n.array(
#        [[49,41,47,19,29,28,34,51],
#         [10, 3,25,48,24,55,27,57],
#         [ 9,58, 1, 4,17,13,56,59],
#         [22,61,35,18, 5,32,30,23]])
#         [20,63,42,37,40,14,54,50],
#         [43, 2,33, 6,52, 7,12,38],
#         [53,21,15,16,62,44, 0,26],
#         [31,45, 8,11,36,60,39,46]])

class CoV:
    '''Covariance class.
        input : Data matrix and bls in set. (X, bls)

        bls   = bls
        X     = data matrix
        nprms = number of channels/kmodes/prms
        C     = covariance of data matrix.
    '''
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
if False: #turning off the auto sep selection in preperation for deletion
    # Get a dict of all separations and the bls that contribute.0001
    #creates a dictionary of separations given a certain baseline
    #and vice versa, gives baselines for a given separation (returns list).
    bl2sep = {}
    sep2bl = {}
    for ri in range(ANTPOS.shape[0]):
        for ci in range(ANTPOS.shape[1]):
            for rj in range(ANTPOS.shape[0]):
                for cj in range(ci,ANTPOS.shape[1]):
                    if ri >= rj and ci == cj: continue # exclude repeat +/- listings of certain bls
                    #sep = a.miriad.ij2bl(rj-ri, cj-ci)
                    sep = (cj-ci, rj-ri) #(dx,dy) in row spacing units
                    i,j = ANTPOS[ri,ci], ANTPOS[rj,cj]
                    bl = a.miriad.ij2bl(i,j)
                    if i > j:
                        i,j = j,i
                        sep = (sep[0]*-1,sep[1]*-1)
                    bl2sep[bl] = sep
                    sep2bl[sep] = sep2bl.get(sep,[]) + [bl]
    #choose unit seperations corresponding to the bls I put in.
    if len(opts.ant.split('_'))>1: #if there are baselines requested
        #get a list of miriad format bl ints
        input_bls = [a.miriad.ij2bl(int(l.split('_')[0]),int(l.split('_')[1])) for l in opts.ant.split(',')]
        print input_bls
        myseps = list(set([bl2sep[bl] for bl in input_bls]))#get a list of the seps, one entry per
        print "based on input baselines, I am including the following seperations"
        print myseps
        mybls = []
        for sep in myseps:
            mybls += sep2bl[sep]
            revsep = (-sep[0],-sep[1])
            mybls += sep2bl[revsep] #don't forget the reverse seps. they count as the same!
        print "found %d baselines"%len(mybls)
        opts.ant = ','.join([miriadbl2str(bl) for bl in mybls])
        print opts.ant
        sys.stdout.flush()
    #WARNING: The default is to do _all_ seps in the data.

# Get a dict of all separations and the bls that contribute.0000
#creates a dictionary of separations given a certain baseline
#and vice versa, gives baselines for a given separation (returns list).
bl2sep = {}
sep2bl = {}
for ri in range(ANTPOS.shape[0]):
    for ci in range(ANTPOS.shape[1]):
        for rj in range(ANTPOS.shape[0]):
            for cj in range(ci,ANTPOS.shape[1]):
                if ri >= rj and ci == cj: continue # exclude repeat +/- listings of certain bls
                #sep = a.miriad.ij2bl(rj-ri, cj-ci)
                sep = (cj-ci, rj-ri) #(dx,dy) in row spacing units
                i,j = ANTPOS[ri,ci], ANTPOS[rj,cj]
                bl = a.miriad.ij2bl(i,j)
                if i > j: 
                    i,j = j,i
                    sep = (sep[0]*-1,sep[1]*-1)
                bl2sep[bl] = sep
                sep2bl[sep] = sep2bl.get(sep,[]) + [bl]
#choose unit seperations corresponding to the bls I put in.

if not opts.usebls:
    if len(opts.ant.split('_'))>1: #if there are baselines requested
        #get a list of miriad format bl ints
        input_bls = [a.miriad.ij2bl(int(l.split('_')[0]),int(l.split('_')[1])) for l in opts.ant.split(',')]
        print input_bls
        myseps = list(set([bl2sep[bl] for bl in input_bls]))#get a list of the seps, one entry per 
        print "based on input baselines, I am including the following seperations"
        print myseps
        mybls = []
        for sep in myseps:
            mybls += sep2bl[sep]
            revsep = (-sep[0],-sep[1])
            mybls += sep2bl[revsep] #don't forget the reverse seps. they count as the same!
        print "found %d baselines"%len(mybls)
        opts.ant = ','.join([miriadbl2str(bl) for bl in mybls])
        print opts.ant
#WARNING: The default is to do _all_ seps in the data.

else:
    print 'Using the following baselines'
    print opts.ant.split(',')
    print 'There are %d of them'%len(opts.ant.split(','))


#checking that our grid indexing is working
print [a.miriad.bl2ij(bl) for bl in sep2bl[(1,0)]]
print len(sep2bl[(0,1)])

uv = a.miriad.UV(args[0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
del(uv)

#creat active frequencies, average active frequency and convert to redshift.
afreqs = freqs.take(chans)
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

#currently not used.
if PFB:
    # XXX unsure how much of a BW modification a windowed PFB needs.  I think not much...
    B = sdf * afreqs.size / NTAPS
    etas = n.fft.fftshift(capo.pspec.f2eta(afreqs[:afreqs.size/NTAPS]))
else:
    #B = sdf * afreqs.size / capo.pfb.NOISE_EQUIV_BW[WINDOW] # this is wrong if we aren't inverting
    # the window post delay transform (or at least dividing out by the gain of the window)
    # For windowed data, the FFT divides out by the full bandwidth, B, which is
    # then squared.  Proper normalization is to multiply by B**2 / (B / NoiseEqBand) = B * NoiseEqBand
    B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW] # normalization. See above.
    etas = n.fft.fftshift(capo.pspec.f2eta(afreqs)) #creat etas (these are fourier dual to frequency)
kpl = etas * capo.pspec.dk_deta(z) #111
bm = n.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq)
scalar = capo.pspec.X2Y(z) * bm * B # cosmological volume scalar
#scalar = 1
print 'Freq:',fq
print 'z:', z
print 'B:', B
print 'scalar:', scalar
sys.stdout.flush()

#cen_fqs = n.arange(.115,.190,.005)
#cen_fqs = n.array([.150])
#kwargs = {'cen_fqs':cen_fqs,'B':B, 'ntaps':NTAPS, 'window':WINDOW, 'bm_fqs':afreqs.clip(.120,.190)}
#window = a.dsp.gen_window(freqs.size, window=WINDOW)

#T is a dictionary of the visibilities and N is the dictionary for noise estimates.
T, N, W = {}, {}, {}
times = []
eor_mdl = {}
for filename in args:
    print 'Reading', filename
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, opts.ant, opts.pol)
    for (crd,t,(i,j)),d,f in uvi.all(raw=True):
        if len(times) == 0 or times[-1] != t:
            if len(times) % 8 == 0:
                #For every 8th integration make random normal noise that is our eor model. Note for every block of 8, this noise is the same. 222
                eor_mdl[t] = n.random.normal(size=chans.size) * n.exp(2j*n.pi*n.random.uniform(size=chans.size))
            else: eor_mdl[t] = eor_mdl[times[-1]]
            times.append(t)
        #For 32 array inside 64 array. skip bls not in the subarray, but may be in the data set.
        if not (( i in ANTPOS ) and ( j in ANTPOS )) : continue
        bl = a.miriad.ij2bl(i,j)
        #sep = bl2sep[bl]
        #print i,j,':',sep
        #if n.abs(sep[0]) != 1 or sep[1] !=0:continue
        #print '-->',i,j
        sys.stdout.flush()
        #if sep[0] < 0:
            #print 'Conj:', a.miriad.bl2ij(bl)
        #    d,sep = n.conj(d),-1*n.array(sep)
        #take active data and convert from janksy's to temperature units. Current data is in janskys.
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
            _Wrms = n.fft.ifft(w)
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
        _Wrms = n.fft.fftshift(_Wrms)
        if False: # swap in a simulated signal post delay transform
            _Trms = n.random.normal(size=_Trms.size) * n.exp(2j*n.pi*n.random.uniform(size=_Trms.size))
            mask = n.ones(_Trms.size); mask[15:25] = 0
            _Trms += .3*eor_mdl[times[-1]] * mask
        #Makes list of visibilities for each baseline for all times. number of integrations by number of channels.
        T[bl] = T.get(bl, []) + [_Trms]
        N[bl] = N.get(bl, []) + [_Nrms]
        W[bl] = W.get(bl, []) + [_Wrms]

#444
n_k = chans.size / NTAPS
bls = T.keys()
for bl in bls:
    print 'bl shape (nints, nchan):'
    T[bl],N[bl] = n.array(T[bl]),n.array(N[bl])
    print '\t',T[bl].shape
if True:
    print 'Fringe-rate filtering the noise to match the data'
    for bl in N:
        _N = n.fft.ifft(N[bl], axis=0) # iffts along the time axis
        _N[23:] = 0 # This was calculated by hand for fr-filter with max_fr=1. and min_fr=0. #555
        N[bl] = n.fft.fft(_N, axis=0)
if False:
    print 'Adding extra noise into the data'
    for bl in bls: T[bl] += N[bl]
#666
Ts = n.concatenate([T[bl] for bl in bls], axis=-1).T
Ns = n.concatenate([N[bl] for bl in bls], axis=-1).T
Ws = n.concatenate([W[bl] for bl in bls],  axis=-1).T
#print bl flagging

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
print Ns.shape
#print times[300], times[500]
print ' '.join(['%d_%d' % a.miriad.bl2ij(bl) for bl in bls])
sys.stdout.flush()
if PLOT:
    #capo.arp.waterfall(cov(Ts), mode='log', drng=2); p.show()
    p.subplot(141); capo.arp.waterfall(Ts, mode='log', mx=1, drng=2); p.colorbar(shrink=.5)
    p.title('Vis in K. bls X ints.', fontsize = 8)
    p.subplot(142); capo.arp.waterfall(Ns, mode='log')#, mx=1, drng=2); p.colorbar(shrink=.5)
    p.title('FRF eor_model.', fontsize=8)
    p.subplot(143); capo.arp.waterfall(Ws); p.colorbar(shrink=0.5)
    p.title('Weights in samples bls x ints', fontsize=8)
    p.subplot(144); capo.arp.waterfall(cov(Ts), mode='log', drng=3); p.colorbar(shrink=.5)
    print cov(Ts).shape
    p.title('cov(Ts)', fontsize=8)
    p.show()
    p.subplot(121); capo.arp.waterfall(cov(Ts), mode='real', mx=.005, drng=.01); p.colorbar(shrink=.5)
    p.title('cov(Ts) real part', fontsize=8)
    #p.subplot(122); capo.arp.waterfall(cov(Ts), mode='log', drng=3); p.colorbar(shrink=.5)
    p.subplot(122); capo.arp.waterfall(cov(Ts), mode='log' ); p.colorbar(shrink=.5)
    p.title('cov(Ts) log', fontsize=8)
    p.tight_layout()
    p.show()

###########################################
#def subtract_average(C, bls, n_k):
def subtract_average(C, bls, n_k, gps):
    '''gps is the list of all groups in all chunks.'''
    print gps
    nbls = len(bls)
    C.shape = (nbls,n_k,nbls,n_k)
    sub_C = n.zeros_like(C) #array to be subtracted from C.
    #choose a (i,j) baseline cross-multiple panel in the covariance matrix
    for i in xrange(nbls):  
        bli = bls[i]
        for j in xrange(nbls):
            blj = bls[j]
            #ensure bli and blj belong to the same group
            gp = same_group(bli, blj, gps)
            if gp is None: continue
            #make sure we only compute average using baselines in the same group
            #Now average over all other panels of covariance matrix (within this group)
            #to get the average signal covariance and subtract that off so that we don't 
            #get signal loss removing residual signal covariances.
            
            _Csum,_Cwgt = 0,0
            for i_ in xrange(nbls):
                bli_ = bls_[i_]
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
                sub_C[i,:,j] = _Csum/_Cwgt
            except(ZeroDivisionError): #catches zero weights
                print 'weights are zero for %d_%d'%(i_,j_)
                sys.stdout.flush()
                sub_C[i,:,j] = _Csum

    C.shape = sub_C.shape = (nbls*n_k, nbls*n_k)
    C -= sub_C
    return C, sub_C

def make_groups(bls,ngps=4):
    #makes groups based on list of bls and number of groups you want.
    #defaults to 4 groups.
    bls_ = random.sample(bls, len(bls))
    nbls = len(bls)
    nblspg = nbls/ngps
    print 'Breaking %d bls into groups of %d'%(nbls, nblspg)
    sys.stdout.flush()
    gps = []
    for gpi in xrange(ngps):
        gps.append(bls_[gpi*nblspg:(gpi+1)*nblspg])
    leftover = nbls - (nblspg*ngps)  
    print gps
    print nbls
    print leftover
    while leftover > 0:
        rc = random.choice(range(ngps))
        gps[rc].append(bls_[-1*leftover])
        leftover-=1
        print leftover
    for gpi in xrange(ngps):
        gps[gpi] = random.sample(gps[gpi],3) + [random.choice(gps[gpi]) for bl in gps[gpi][:len(gps[gpi])-3]]
    return gps

def same_group(bli, blj, gps):
    for gp in gps:
        if bli in gp and blj in gp:
            return gp
    return None

############################################
for boot in xrange(NBOOT):
    #777
    if True: # pick a sample of baselines with replacement
        bls_ = random.sample(bls, len(bls))
        nbls = len(bls)
        #nchunks = 2 
        nchunks = 1
        chunks = []
        #Split baselines into chunks that will be diagonalized separately.
        for nc in xrange(nchunks):
            chunks.append(bls_[nc*(nbls/nchunks):(nc+1)*(nbls/nchunks)])
        #chunks does not have repeated baselines.
        chgps = []
        #chgps = chunkgroups.
        #each element of chgps is the ngps used in covariance diag for each of the 
        #chunks. eg2 [ [[gpi1,gpi2,gpi3]], [[gpj1,gpj2,gpj3]] ] for chunks i,j
        #baselines may be repeated within a chunk. 
        for chunk in chunks:
            chgps.append(make_groups(chunk))
        print 'chgps: ', chgps         
        
        #gp1,gp2,gp3,gp4 = bls_[:14],bls_[14:28],bls_[28:42],bls_[42:] # for 56bl -> 14 * 4
        #gp1,gp2,gp3,gp4 = bls_[:5],bls_[5:10],bls_[10:15],bls_[15:] # for 21bl
        # ensure each group has at least 3 kinds of baselines. Otherwise get 0 divide.
        #gp1 = random.sample(gp1, 3) + [random.choice(gp1) for bl in gp1[:len(gp1)-3]]
        #gp2 = random.sample(gp2, 3) + [random.choice(gp2) for bl in gp2[:len(gp2)-3]]
        #gp3 = random.sample(gp3, 3) + [random.choice(gp3) for bl in gp3[:len(gp3)-3]]
        #gp4 = random.sample(gp4, 3) + [random.choice(gp4) for bl in gp4[:len(gp4)-3]]
    #else:
    #    bls_ = random.sample(bls, len(bls))
    #    gp1,gp2 = bls_[:len(bls)/2],bls_[len(bls)/2:]
    #gp2 = gp2[:len(gp1)] # XXX force gp1 and gp2 to be same size
    #gp1,gp2 = gp1+gp2,[] # XXX
    #print 'XXX', len(gp1), len(gp2)

    #get bls used within a chunk
    chbls_ = []    
    for i in xrange(nchunks):
        cb = []
        for cbc in chgps[i]:cb+=cbc
        chbls_.append(cb)
    #get list of all baselines. Record keeping.
    bls_all = []
    for i in xrange(nchunks):
        bls_all.append(chbls_[i]) 
    bls_all = bls_all[0]
    print 'bls_all: ', bls_all
    all_gps = []
    for i in xrange(nchunks):
        all_gps += chgps[i]
    print 'all_gps: ', all_gps

    print 'Bootstrap sample %d:' % boot,
    for chunk in chbls_: print '(%s)' % (','.join(['%d_%d' % a.miriad.bl2ij(bl) for bl in chunk])),
    print
    #again Ts is the number of bls*channels X number of integrations. Note this rearragnes the order of bls in Ts to be that of bls_. May be repititions.
    Ts = n.concatenate([T[bl] for bl in bls_all], axis=1).T
    Ns = n.concatenate([N[bl] for bl in bls_all], axis=1).T
    #Ns = n.concatenate([N[bl] for bl in bls_], axis=1).T # this noise copy processed as if it were the data
    #CTs.append(n.concatenate([T[bl] for chbl in chbls_ for bl in chbl], axis=1).T)
    #CNs.append(n.concatenate([N[bl] for chbl in chbls_ for bl in chbl], axis=1).T)
    L = len(bls_all)#total size of both chunks
    temp_noise_var = n.var(n.array([T[bl] for bl in bls_all]), axis=0).T
    #
    #temp_noise_var = n.average(n.array([T[bl] for bl in bls_]), axis=0).T
    #print Ts.shape, temp_noise_var.shape
    sys.stdout.flush()

    _Cxtot,_Cntot = 1, 1
    #PLT1,PLT2 = 4,4
    if opts.niters:
        #override number if iters.
        PLT1,PLT2 = map(int, opts.niters.split(','))
    else:
        PLT1,PLT2 = int(3*n.sqrt(0.3/opts.gain)),int(3*n.sqrt(0.3/opts.gain))#scale the number of steps by the gain? -dcj
    #PLT1,PLT2 = 2,2
    #PLT1,PLT2 = 1,2
    #888
    for cnt in xrange(PLT1*PLT2-1):
        print cnt, '/', PLT1*PLT2-1
        SZ = Ts.shape[0]
        #MCx,MCn = CoV(Ts, bls_), CoV(Ns, bls_)
        MCx,MCn = CoV(Ts, bls_all), CoV(Ns, bls_all)
        #Cx,Cn = cov(Ts), cov(Ns)
        Cx,Cn = MCx.C, MCn.C
        #to turn 64 data into 2 "32" data sets by masking off diagonals.
        #Big_Mask = n.ones_like(Cx, n.dtype=n.float32)
        
        #one_half = ( len(gp1) + len(gp2) ) * n_k #nbls in half times nbins.
        #other_half = ( len(gp3) + len(gp4) ) * n_k #nbls in other have times nbins.
        #Big_Mask[:one_half, one_half+1:] = 0.
        #Big_Mask[one_half+1:, :one_half] = 0.
        #Cx *= Big_Mask
        #Cn *= Big_Mask
        
        if PLOT:
            if cnt%1==0:
#                p.figure(7)
#                capo.arp.waterfall(_Cx, mode='log', mx=8,drng=4); p.colorbar(shrink=.5)
                p.figure(7)
                p.clf()
                capo.arp.waterfall(Cx*scalar, mode='log', mx=8,drng=4); p.colorbar(shrink=.5)
                #p.savefig('temp/cov%d'%cnt)
                p.show()
                #p.figure(7)
                #p.clf()
                #capo.arp.waterfall(Cx[3*n_k:4*n_k,1*n_k:2*n_k]*scalar, mode='log', mx=8,drng=4); p.colorbar(shrink=.5)
                #p.savefig('temp/cov_1pair%d'%cnt)
#                p.figure(8)
#                capo.arp.waterfall(avg_Cx, mode='log', mx=8,drng=4); p.colorbar(shrink=.5)
#                p.figure(9)
#                capo.arp.waterfall(_Cx -  avg_Cx, mode='log', mx=8,drng=4); p.colorbar(shrink=.5)
            #capo.arp.waterfall(cov(Ts), mode='log', mx=-1,  drng=4); p.colorbar(shrink=.5)
            #p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(cov(Ts), mode='log', mx=0,  drng=3); p.colorbar(shrink=.5)
            #p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(cov(Ns), mode='log', mx=0, drng=2)
            print "max(cov(Ts))",n.max(Cx)
            sys.stdout.flush()
        print "max(cov(Ts))",n.max(Cx)
        #999
        #SZ is number of rows in Ts. i.e. #bls * #channels.
        dx = n.diag(Cx); dx.shape = (1,SZ); Cx /= dx
        dn = n.diag(Cn); dn.shape = (1,SZ); Cn /= dn

        #g = .3 # for 1*7 baselines
        g = opts.gain # for 4*7 baselines
        print 'gain factor = ', g
        # begin with off-diagonal covariances to subtract off
        # (psuedo-inv for limit of small off-diagonal component)
#        _Cx,_Cn = -g*_Cx, -g*_Cn
        _Cx,_Cn = -g*Cx, -g*Cn
        ind = n.arange(SZ)
        # XXX do we also need to zero modes adjacent to diagonal, since removing them results in bigger signal loss?
        for b in xrange(L): # for each redundant baseline, zero out diagonal from covariance diagonalization. Sets each diagonal of each bl-bl covariance to 0.
            indb = ind[:-b*n_k]
            _Cx[indb,indb+b*n_k] = 0 
            _Cx[indb+b*n_k,indb] = 0
            _Cn[indb,indb+b*n_k] = 0 
            _Cn[indb+b*n_k,indb] = 0
        _Cx[ind,ind] = 0 
        _Cn[ind,ind] = 0 # set these to zero temporarily to avoid noise bias into cross terms


        #subtract average
        _Cx, avg_Cx = subtract_average(_Cx, bls_all, n_k, all_gps)
        _Cn, avg_Cn = subtract_average(_Cn, bls_all, n_k, all_gps)


######This is the subtract the average part which has been moved to a function above.######

        #if True: # estimate and remove signal covariance from diagonalization process
            # do this twice: once for the signal (Cx) and once for the noise (Cn)
            # using the statistics of the signal and noise, respectively
            #for _C in [_Cx,_Cn]:


#            for _C in [_Cx]:#,_Cn]:
#                #remember L is the total number of baselines and n_k is the number of kbins (i.e. number of channels). Note shape of covariance is n_k*#bls X n_k*#bls.
#                _C.shape = (L,n_k,L,n_k)
#                sub_C = n.zeros_like(_C)
#                # Choose a (i,j) baseline cross-multiple panel in the covariance matrix
#                for i in xrange(L):
#                    bli = bls_[i]
#                    for j in xrange(L):
#                        blj = bls_[j]
#                        # even remove signal bias if bli == blj
#                        # ensure bli and blj belong to the same group
#                        if bli in gp1 and blj in gp1: gp = gp1
#                        elif bli in gp2 and blj in gp2: gp = gp2
#                        elif bli in gp3 and blj in gp3: gp = gp3
#                        elif bli in gp4 and blj in gp4: gp = gp4
#                        else: continue # make sure we only compute average using bls in same group
#                        # Now average over all other panels of covariance matrix (within this group)
#                        # to get the average signal covariance and subtract that off so that we don't
#                        # get signal loss removing residual signal covariances.
#                        #AAA, Why are we not checking if the baselines are in the same group as the one we want to subtract from? i.e. bli_ and blj_ in gp{i}
#                        #CHANGE TO GP #Seems like it needs to be for i_,j_ in gp:
#                        _Csum,_Cwgt = 0,0
#                        for i_ in xrange(L):
#                            #check if i_ in gp
#                            bli_ = bls_[i_]
#                            if not bli_ in gp: continue # make sure averaging over baseline in the same group.
#                            if bli == bli_: continue # only average over other bls to better isolate bl systematics
#                            for j_ in xrange(L):
#                                blj_ = bls_[j_]
#                                if not blj_ in gp: continue # make sure averaging over baseline in the same group.
#                                if bli_ == blj_: continue # don't average over panels with noise bias
#                                if blj == blj_: continue # only average over other bls to better isolate bl systematics
#                                _Csum += _C[i_,:,j_] # fixes indexing error in earlier ver
#                                _Cwgt += 1
#                        try:
#                            sub_C[i,:,j] = _Csum / _Cwgt # fixes indexing error in earlier ver
#                        except(ZeroDivisionError): #catches zero weights
#                           # print gp
#                           # print i,j,bli,blj
#                           # print i_,j_,bli_,blj_
#                           # print _Cwgt
#                            print 'weights are zero for %d_%d'%(i_,j_)
#                            sys.stdout.flush()
#                            sub_C[i,:,j] = _Csum
#                _C.shape = sub_C.shape = (L*n_k,L*n_k)
#                _C -= sub_C


#######This is the subtract the average part which has been moved to a function above.END######


        if PLOT:
            if cnt%1==0:
                p.figure(100)
                p.clf()
                #capo.arp.waterfall(avg_Cx*scalar*dx/g, mode='log', mx=8, drng=4); p.colorbar(shrink=.5)
                capo.arp.waterfall(avg_Cx, mode='log'); p.colorbar(shrink=.5)
                #p.savefig('temp/avg_cov%d'%cnt)
                p.figure(10)
                p.clf()
                capo.arp.waterfall(_Cx*scalar*dx/g, mode='log', mx=8, drng=4); p.colorbar(shrink=.5)
                #p.savefig('temp/cov_minus_avg_cov%d'%cnt)

               # p.figure(100)
               # p.clf()
               # CCCX = avg_Cx*scalar*dx/g
               # _CCCX = _Cx*scalar*dx/g
               # capo.arp.waterfall(CCCX[3*n_k:4*n_k,1*n_k:2*n_k], mode='log', mx=8, drng=4); p.colorbar(shrink=.5)
               # p.savefig('temp/avg_cov1pair%d'%cnt)
               # p.figure(10)
               # p.clf()
               # capo.arp.waterfall(_CCCX[3*n_k:4*n_k,1*n_k:2*n_k], mode='log', mx=8, drng=4); p.colorbar(shrink=.5)
               # p.savefig('temp/cov_minus_avg_cov1pair%d'%cnt)
                p.show()
            #p.subplot(131);capo.arp.waterfall(sub_C, mode='log',mx=-1, drng=4); p.colorbar(shrink=.5)
                #p.subplot(132);capo.arp.waterfall(_Cx, mode='log',mx=-1, drng=4); p.colorbar(shrink=.5)
                #p.subplot(133);capo.arp.waterfall(-g*Cx, mode='log',mx=-1, drng=4); p.colorbar(shrink=.5)
#                        p.show()

        if True:
            # divide bls into two independent groups to avoid cross-contamination of noise
            # this is done by setting mask=0 for all panels pairing bls between different groups
            # this masks covariances between groups, so no covariance from a bl in one group is subtracted from
            # a bl in another group
            mask = n.ones_like(Cx)
            # XXX need to clean this section up
            for i, bli in enumerate(bls_all):
                for j, blj in enumerate(bls_all):
                    if not same_group(bli, blj, all_gps) is None: continue
                    mask[i*n_k:(i+1)*n_k,j*n_k:(j+1)*n_k] = 0
            
#            #for bl1 in xrange(len(gp1)):
#            #    for bl2 in xrange(len(gp1)):
#            #        if bls_[bl1] != bls_[bl2]: continue # zero out panels where bl1 == bl2
#            #        mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
#            #        mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
#            for bl1 in xrange(len(gp1)):
#                for bl2 in xrange(len(gp2)):
#                    bl2 += len(gp1)
#                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
#                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
#            for bl1 in xrange(len(gp1)):
#                for bl2 in xrange(len(gp3)):
#                    bl2 += len(gp1) + len(gp2)
#                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
#                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
#            for bl1 in xrange(len(gp1)):
#                for bl2 in xrange(len(gp4)):
#                    bl2 += len(gp1) + len(gp2) + len(gp3)
#                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
#                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
#            #for bl1 in xrange(len(gp2)):
#            #    bl1 += len(gp1)
#            #    for bl2 in xrange(len(gp2)):
#            #        bl2 += len(gp1)
#            #        if bls_[bl1] != bls_[bl2]: continue # zero out panels where bl1 == bl2
#            #        mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
#            #        mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
#            for bl1 in xrange(len(gp2)):
#                bl1 += len(gp1)
#                for bl2 in xrange(len(gp3)):
#                    bl2 += len(gp1) + len(gp2)
#                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
#                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
#            for bl1 in xrange(len(gp2)):
#                bl1 += len(gp1)
#                for bl2 in xrange(len(gp4)):
#                    bl2 += len(gp1) + len(gp2) + len(gp3)
#                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
#                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
#            #for bl1 in xrange(len(gp3)):
#            #    bl1 += len(gp1) + len(gp2)
#            #    for bl2 in xrange(len(gp3)):
#            #        bl2 += len(gp1) + len(gp2)
#            #        if bls_[bl1] != bls_[bl2]: continue # zero out panels where bl1 == bl2
#            #        mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
#            #        mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
#            for bl1 in xrange(len(gp3)):
#                bl1 += len(gp1) + len(gp2)
#                for bl2 in xrange(len(gp4)):
#                    bl2 += len(gp1) + len(gp2) + len(gp3)
#                    mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
#                    mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
#            #for bl1 in xrange(len(gp4)):
#            #    bl1 += len(gp1) + len(gp2) + len(gp3)
#            #    for bl2 in xrange(len(gp4)):
#            #        bl2 += len(gp1) + len(gp2) + len(gp3)
#            #        if bls_[bl1] != bls_[bl2]: continue # zero out panels where bl1 == bl2
#            #        mask[bl1*n_k:(bl1+1)*n_k,bl2*n_k:(bl2+1)*n_k] = 0
#            #        mask[bl2*n_k:(bl2+1)*n_k,bl1*n_k:(bl1+1)*n_k] = 0
#            #BBB All of the above for loops mask the gp-gp baseline pairs. Get a diagonal matrix of covariances within the group.
            _Cx *= mask; _Cn *= mask
        #make diagonal 1 after applying mask.
        _Cx[ind,ind] = _Cn[ind,ind] = 1
        Ts,Ns = n.dot(_Cx,Ts), n.dot(_Cn,Ns)
        # These are a running tally of all the diagonalization steps applied
        _Cxtot,_Cntot = n.dot(_Cx,_Cxtot), n.dot(_Cn,_Cntot)
#        p.figure(cnt)
#        p.subplot(111); capo.arp.waterfall(Ts, mode='log', drng=3);p.colorbar(shrink=.5)
    if PLOT:
        capo.arp.waterfall(cov(Ts), mode='log', mx=-1, drng=4);p.colorbar(shrink=.5)
        #p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ts), mode='log', mx=0, drng=3);p.colorbar(shrink=.5)
        #p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ns), mode='log', mx=0, drng=3)
        p.show()

#    p.show()
#    import IPython
#    IPython.embed()
#    exit()
    #Ts = n.concatenate([T[bl] for bl in bls_], axis=1).T
    #Ns = n.concatenate([N[bl] for bl in bls_], axis=1).T # this noise copy processed as if it were the data
    Ts = n.concatenate([T[bl] for bl in bls_all], axis=1).T
    Ns = n.concatenate([N[bl] for bl in bls_all], axis=1).T # this noise copy processed as if it were the data

    pspecs,dspecs = [], []
    nspecs,n1specs,n2specs = [], [], []
    #Cx,Cn = CoV(Ts, bls_), CoV(Ns, bls_)
    #Cx_ = CoV(n.dot(_Cxtot,Ts), bls_)
    Cx,Cn = CoV(Ts, bls_), CoV(Ns, bls_all)
    Cx_ = CoV(n.dot(_Cxtot,Ts), bls_all)
    # Cn1 is the noise diagonalized as if it were the signal, Cn2 is the noise with the signal diagonalization applied
    #Cn1_,Cn2_ = CoV(n.dot(_Cntot,Ns), bls_), CoV(n.dot(_Cxtot,Ns), bls_)
    Cn1_,Cn2_ = CoV(n.dot(_Cntot,Ns), bls_all), CoV(n.dot(_Cxtot,Ns), bls_all)
    bls_ = bls_all #THis is new. Use with the chunnk stuff.
    for cnt,bli in enumerate(bls_):
        print cnt
        for blj in bls_[cnt:]:
            #print a.miriad.bl2ij(bli), a.miriad.bl2ij(blj)
            # XXX behavior here is poorly defined for repeat baselines in bootstrapping
            xi,xj = Cx.get_x(bli), Cx.get_x(blj) # No covariance applied.
            xi_,xj_ = Cx_.get_x(bli), Cx_.get_x(blj) #covariance applied.
            pk_avg = scalar * xi * xj.conj() # make a power spectrum from bli*blj^*.
            dspecs.append(pk_avg) # do this before bli == blj check to include noise bias in dspec
            if bli == blj: continue
            #if True: # exclude intra-group pairings # XXX
            #    if (bli in gp1 and blj in gp1) or (bli in gp2 and blj in gp2) or (bli in gp3 and blj in gp3) or (bli in gp4 and blj in gp4): continue
            if not same_group(bli, blj, all_gps) is None: continue
            if not opts.noproj: # do an extra final removal of leakage from particular modes
                print 'Projecting'
                Ts = n.concatenate([xi_,xj_], axis=0)
                cx = cov(Ts)
                #if PLOT:
                if False:
                    p.clf()
                    p.subplot(121); capo.arp.waterfall(cx, mode='log', mx=0, drng=3)
                for cnt1 in xrange(9):
                    d = n.diag(cx); d.shape = (1,d.size); cx /= n.sqrt(d) * 2
                    g = .3
                    _cx = -g*cx
                    mask = n.zeros_like(cx)
                    #CCC what is prj_ch
                    if n_k == 20: prj_ch = xrange(8,12)
                    elif n_k == 40: prj_ch = xrange(17,24)
                    elif n_k == 80: prj_ch = xrange(34,48)
                    else: #raise ValueError('Only support # channels = (20,40,80) for now')
                        prj_ch = xrange(int(n_k/2-n_k/10),int(n_k/2+n_k/10))
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
                    p.subplot(122); capo.arp.waterfall(cx, mode='log', mx=0, drng=3)
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
            sys.stdout.flush()
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
    print pspecs.shape
#    import pylab as p
#    p.figure(45)
    #capo.arp.waterfall(n.mean(pspecs,axis=-1), mode='log');p.colorbar(shrink=.5)
#    capo.arp.waterfall(pspecs[:,:,30], mode='log');p.colorbar(shrink=.5)
#    p.show()
    nspecs,n1specs,n2specs = n.array(nspecs), n.array(n1specs), n.array(n2specs)
    navg_2d = n.average(nspecs, axis=0)
    n1avg_2d = n.average(n1specs, axis=0)
    n2avg_2d = n.average(n2specs, axis=0)
    f1 = n.sqrt(n.sum(n.abs(navg_2d)**2)/n.sum(n.abs(n1avg_2d)**2))
    f2 = n.sqrt(n.sum(n.abs(navg_2d)**2)/n.sum(n.abs(n2avg_2d)**2))
    print 'Rescale factor (FINAL):', f1, f2
    sys.stdout.flush()

    avg_2d = n.average(pspecs, axis=0) # average over baseline cross-multiples
    std_2d = n.std(pspecs, axis=0) # get distribution as a function of time
    wgt_2d = 1. / std_2d**2 # inverse variance weighting
    avg_1d = n.sum(avg_2d * wgt_2d, axis=1) / n.sum(wgt_2d, axis=1)

    if PLOT:
        import capo as C
        p.subplot(131)
        C.arp.waterfall(avg_2d, mode='log', drng=3);p.colorbar(shrink=.5)
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

    outfile = 'pspec_boot%04d.npz'%(boot)
    if not opts.output == '':
        outfile =opts.output+'/'+outfile
    print "Writing", outfile
    n.savez(outfile, kpl=kpl, pk=avg_1d, err=std_1d, scalar=scalar, times=n.array(times),freq=fq,
        pk_vs_t=avg_2d, err_vs_t=std_2d, temp_noise_var=temp_noise_var, nocov_vs_t=n.average(dspecs,axis=0),
        cmd=' '.join(sys.argv))
if PLOT: p.show()

