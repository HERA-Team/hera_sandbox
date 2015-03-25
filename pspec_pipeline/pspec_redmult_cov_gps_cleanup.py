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
o.add_option('--ngps', type='int', default=4,
    help='Number of groups. Default is 4.')
o.add_option('--boot_number', type='int', 
    help='Bootstrap number to do. Used with qsub')
o.add_option('--noise', action='store_true',
    help='use noise uv files.')
o.add_option('--write2uv', action='store_true',
    help='Instead of forming power spectra, write to uv files after \
          applying covariance matrix.')
o.add_option('--applycov', action='store', 
    help='Apply covariance matrices given in these boot npz files. \
          Should have name cov_matrix.')
o.add_option('--savecov', action='store_true',
    help='Save covariance matrices')
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

class CoV:
    '''Covariance class.
        input : Data matrix and bls in set. (X, bls)

        bls   = bls
        X     = data matrix
        nprms = number of channels/kmodes/prms
        C     = covariance of data matrix.
    '''
    def __init__(self, X, bls, times):
        self.bls = bls
        self.X = X
        self.nprms = X.shape[0] / len(bls)
        self.C = cov(X)
        self.times = n.array(times)
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


def grid2ij(GRID):
    '''
        bl_str = given sep, returns bls in string format.
        bl_conj = given a baseline (miriad bl) gives separation.
        bl2sep_str = given baseline (miriad) return its separation.    
    '''
    bls, conj = {}, {}
    for ri in range(GRID.shape[0]):
        for ci in range(GRID.shape[1]):
            for rj in range(GRID.shape[0]):
                for cj in range(GRID.shape[1]):
                    if ci > cj: continue
#                    if ri > rj and ci == cj: continue
#                    if ci > cj and ri == rj: continue
                    sep = (rj-ri, cj-ci)
                    sep = '%d,%d'%sep
                    i,j = GRID[ri, ci], GRID[rj,cj]
                    bls[sep] = bls.get(sep,[]) + [(i,j)]
    for sep in bls.keys():
        if sep == '0,0' or len(bls[sep]) < 2 or (sep[-1] == '0' and sep[0] == '-'): del(bls[sep
])
    for sep in bls:
        conj[sep] = [i>j for i,j in bls[sep]]

    bl_str,bl_conj,bl2sep_str = {}, {}, {}
    for sep in bls:
        bl_str[sep],bl_list = [], []
        for (i,j),c in zip(bls[sep],conj[sep]):
            if c: i,j = j,i
            bl_list.append(a.miriad.ij2bl(i,j))
            bl_str[sep].append('%d_%d'%(i,j))
            bl2sep_str[a.miriad.ij2bl(i,j)] = bl2sep_str.get(a.miriad.ij2bl(i,j),'') + sep
            bl_conj[a.miriad.ij2bl(i,j)] = c
        bls[sep] = bl_list
        bl_str[sep] = ','.join(bl_str[sep])
    return bl_str,bl_conj,bl2sep_str



if not opts.ant:
    print 'NO BASELINES. PLEASE INPUT BASELINES ON COMMAND LINE.'
    exit()
else:
    print 'Using the following baselines'
    print opts.ant.split(',')
    print 'There are %d of them'%len(opts.ant.split(','))

sep2bl, conjbl, bl2sep = grid2ij(ANTPOS)


#checking that our grid indexing is working

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


def read_noise_uv(files, ants, pol, chs, conj, win):
    print 'Reading noise files'
    print files
    NRMS = {}
    window = a.dsp.gen_window(len(chs), win)
    for fname in files:
        uvi = a.miriad.UV(fname)
        a.scripting.uv_selector(uvi, ants, pol)
        for (crd,t,(i,j)),d,f in uvi.all(raw=True):
            bl = a.miriad.ij2bl(i,j)
            if conj[bl]:
                d = n.conj(d)
            d,f = d.take(chans), f.take(chans)
            w = n.logical_not(f).astype(n.float)            
            #window /= np.sum(window)
            Nfft = n.fft.ifft(window * d* w)
            Nffts = n.fft.fftshift(Nfft)
            NRMS[bl] = NRMS.get(bl, []) + [Nffts]
    return NRMS

def write_2_uv(files, CC, afreqs, channels=(110,150), ending='V'):
    '''Give input uv files. CC = covariance class that'''
    window = a.dsp.gen_window(len(afreqs), WINDOW)

    def mfunc(uv, p, d, f):
        uvw, t, (i,j) = p
        print (i,j)
        bl = a.miriad.ij2bl(i,j)
        print t
        ti = n.where(CC.times == t); print ti
#        import IPython
#        IPython.embed()
        data = n.zeros(len(d), n.complex64).flatten()
#        print data.shape
#        print CC.get_x(bl).shape
        try:
            data[channels[0]:channels[1]] = n.fft.fft(n.fft.ifftshift(CC.get_x(bl)[:,ti].flatten()))/(window*capo.pspec.jy2T(afreqs))
            
        except:
            data = None
            flags = None
        flags = n.ones(len(f), dtype=f.dtype)
        flags[channels[0]:channels[1]] = 0
        return p, data, flags

    for ff in files:
        outfile = ff + ending
        print 'writing to %s'%outfile
        uvi = a.miriad.UV(ff)
        uvo = a.miriad.UV(outfile, status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc, append2hist='covariance applied\n', raw=True)
          
        
    

            

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
            #newrandnoise1 = n.random.normal()*n.exp(2*n.pi*1j*n.random.uniform())
            #newrandnoise2 = n.random.normal()*n.exp(2*n.pi*1j*n.random.uniform())
            #noise_tm_fq = n.random.normal(size=chans.size) * n.exp(2j*n.pi*n.random.uniform(size=chans.size))
            if len(times) % 44 == 0:
                #For every 8th integration make random normal noise that is our eor model. Note for every block of 8, this noise is the same. 222
                eor_mdl[t] = n.random.normal(size=chans.size) * n.exp(2j*n.pi*n.random.uniform(size=chans.size))
            else: eor_mdl[t] = eor_mdl[times[-1]]
            times.append(t)

            #For 32 array inside 64 array. skip bls not in the subarray, but may be in the data set.
        if not (( i in ANTPOS ) and ( j in ANTPOS )) : continue
        bl = a.miriad.ij2bl(i,j)
        sys.stdout.flush()
        if conjbl[bl]:
#            print 'Conjugating baseline (%d,%d) or %d'%(i,j,bl)
            d= n.conj(d)
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
    
        #XXX
        #_Trms[20] += newrandnoise1*.01
        #_Trms[21] += (newrandnoise1*.5 + newrandnoise2*.5)*.01

        if False: # swap in a simulated signal post delay transform
            _Trms = n.random.normal(size=_Trms.size) * n.exp(2j*n.pi*n.random.uniform(size=_Trms.size))
            mask = n.ones(_Trms.size); mask[15:25] = 0
            _Trms += .3*eor_mdl[times[-1]] * mask
        #Makes list of visibilities for each baseline for all times. number of integrations by number of channels.
        T[bl] = T.get(bl, []) + [_Trms]
        N[bl] = N.get(bl, []) + [_Nrms]
        W[bl] = W.get(bl, []) + [_Wrms]


if opts.noise:
    #overwrite the noise above with that from the noise files.
    print 'Over writing noise from uv noise files'
    noise_files = [f+'_noiseL' for f in args]
    N = read_noise_uv(noise_files, opts.ant, opts.pol, chans, conjbl, opts.window)

#444
n_k = chans.size / NTAPS
bls = T.keys()
for bl in bls:
    print 'bl shape (nints, nchan):'
    T[bl],N[bl] = n.array(T[bl]),n.array(N[bl])
    print '\t',T[bl].shape
#if True:
if False:
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
    while leftover > 0:
        rc = random.choice(range(ngps))
        gps[rc].append(bls_[-1*leftover])
        leftover-=1
    for gpi in xrange(ngps):
        gps[gpi] = random.sample(gps[gpi],3) + [random.choice(gps[gpi]) for bl in gps[gpi][:len(gps[gpi])-3]]
    return gps

def same_group(bli, blj, gps):
    for gp in gps:
        if bli in gp and blj in gp:
            return gp
    return None

############################################

if opts.applycov:
    import glob
    covfiles = n.sort(glob.glob(opts.applycov))
    print covfiles

for boot in xrange(NBOOT):
    #continue with this itertion
    if not opts.boot_number: pass
    elif opts.boot_number != boot:
        continue 
    #if using previous covariance matrix set variables and skip 
    #diagonalization process.
    if opts.applycov:
        print "Skipping Covariance Diagonalization and over writing bls_, _Cxtot, _Cntot, gps"
        cov_file = n.load(covfiles[boot])
        print "Reading from %s" %covfiles[boot]
        bls_ = list(cov_file['bls']) #need list to work well with CoV
        _Cxtot = cov_file['cov_matrix']
        _Cntot = cov_file['cov_matrix_noise']
        gps = cov_file['gps']
#    if True: # pick a sample of baselines with replacement
    else: # pick a sample of baselines with replacement
        bls_ = random.sample(bls, len(bls))
        nbls = len(bls)
        gps = make_groups(bls_, ngps=opts.ngps)
        bls_=[]
        for gp in gps : bls_+=gp
        
        print 'Bootstrap sample %d:' % boot,
        for gp in gps: print '(%s)' % (','.join(['%d_%d' % a.miriad.bl2ij(bl) for bl in gp])),
        print
        #again Ts is the number of bls*channels X number of integrations. Note this rearragnes the order of bls in Ts to be that of bls_. May be repititions.
        Ts = n.concatenate([T[bl] for bl in bls_], axis=1).T
        Ns = n.concatenate([N[bl] for bl in bls_], axis=1).T # this noise copy processed as if it were the data
        L = len(bls_)
        #temp_noise_var = n.var(n.array([T[bl] for bl in bls_]), axis=0).T
        temp_noise_var = n.average(n.array([T[bl] for bl in bls_]), axis=0).T
        print Ts.shape, temp_noise_var.shape
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
            MCx,MCn = CoV(Ts, bls_, times), CoV(Ns, bls_, times)
            #Cx,Cn = cov(Ts), cov(Ns)
            Cx,Cn = MCx.C, MCn.C
            if PLOT:
                #if cnt%10==0:
                #    p.figure(7)
                #    capo.arp.waterfall(Cx*scalar, mode='log', mx=8,drng=4); p.colorbar(shrink=.5)
                #capo.arp.waterfall(cov(Ts), mode='log', mx=-1,  drng=4); p.colorbar(shrink=.5)
                #p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(Cx*scalar, mode='log', mx=8,  drng=4); p.colorbar(shrink=.5)
                p.figure(10)
                #p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(Cn*scalar, mode='log', mx=8,  drng=4); p.colorbar(shrink=.5)
                p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(Cx*scalar, mode='log', mx=8,  drng=4); p.colorbar(shrink=.5)
                #p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(cov(Ns), mode='log', mx=0, drng=2)
                print "max(cov(Ts))",n.max(Cx)
                print "max(cov(Ns))",n.max(Cn)
                sys.stdout.flush()
            print "max(cov(Ts))",n.max(Cx)
            print "max(cov(Ns))",n.max(Cn)
            #999
            #for c in [Cx,Cn]: # Normalize covariance matrices
            dx = n.copy(n.diag(Cx)); dx.shape = (1,SZ); Cx /= dx
            dn = n.copy(n.diag(Cn)); dn.shape = (1,SZ); Cn /= dn
            #g = .3 # for 1*7 baselines
            g = opts.gain # for 4*7 baselines
            print 'gain factor = ', g
            # begin with off-diagonal covariances to subtract off
            # (psuedo-inv for limit of small off-diagonal component)
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
            #_Cx, avg_Cx = subtract_average(_Cx, bls_, n_k, [gp1,gp2,gp3,gp4])
            #_Cn, avg_Cn = subtract_average(_Cn, bls_, n_k, [gp1,gp2,gp3,gp4])
            _Cx, avg_Cx = subtract_average(_Cx, bls_, n_k, gps)
            _Cn, avg_Cn = subtract_average(_Cn, bls_, n_k, gps)
    ######This is the subtract the average part which has been moved to a function above.######



    #        if True: # estimate and remove signal covariance from diagonalization process
    #            # do this twice: once for the signal (Cx) and once for the noise (Cn)
    #            # using the statistics of the signal and noise, respectively
    #            #for _C in [_Cx,_Cn]:
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
                pass
    #            if cnt%10==0:
    #                p.figure(100)
                p.figure(11)
                #p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(avg_Cn*scalar, mode='log', mx=8,  drng=4); p.colorbar(shrink=.5)
                p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(avg_Cx*scalar, mode='log', mx=8,  drng=4); p.colorbar(shrink=.5)
    #            #correct for diagonal, scalar, and gain factor
    #                capo.arp.waterfall(avg_Cx*scalar*dx/g, mode='log', mx=8, drng=4); p.colorbar(shrink=.5)
    #                p.figure(99)
    #                p.plot(dx.flatten())
                #p.subplot(131);capo.arp.waterfall(sub_C, mode='log',mx=-1, drng=4); p.colorbar(shrink=.5)
                #p.subplot(132);capo.arp.waterfall(_Cx, mode='log',mx=-1, drng=4); p.colorbar(shrink=.5)
                #p.subplot(133);capo.arp.waterfall(-g*Cx, mode='log',mx=-1, drng=4); p.colorbar(shrink=.5)
    #                p.show()
                #p.figure(1)

            if True:
                # divide bls into two independent groups to avoid cross-contamination of noise
                # this is done by setting mask=0 for all panels pairing bls between different groups
                # this masks covariances between groups, so no covariance from a bl in one group is subtracted from
                # a bl in another group
               mask = n.ones_like(Cx)
               for i, bli in enumerate(bls_):
                    for j, blj in enumerate(bls_):
                        if not same_group(bli, blj, gps) is None: continue
                        mask[i*n_k:(i+1)*n_k,j*n_k:(j+1)*n_k] = 0
               _Cx *= mask; _Cn *= mask
            #make diagonal 1 after applying mask.
            _Cx[ind,ind] = _Cn[ind,ind] = 1
            Ts,Ns = n.dot(_Cx,Ts), n.dot(_Cn,Ns)
            # These are a running tally of all the diagonalization steps applied
            _Cxtot,_Cntot = n.dot(_Cx,_Cxtot), n.dot(_Cn,_Cntot)
            p.figure(12)
            p.subplot(PLT1,PLT2,cnt+1); capo.arp.waterfall(_Cxtot*scalar, mode='log', mx=8, drng=4);p.colorbar(shrink=.5)
    #        p.figure(cnt)
    #        p.subplot(111); capo.arp.waterfall(Ts, mode='log', drng=3);p.colorbar(shrink=.5)
        if opts.savecov:
            print 'Saving covariance matrix'
            cov_outfile = 'pspec_cov%04d.npz'%boot
            print 'Writing to %s'%cov_outfile
            n.savez(cov_outfile, kpl=kpl, scalar=scalar, times=n.array(times),freq=fq, cov_matrix=_Cxtot, cov_matrix_noise=_Cntot, cov_matrix_orig_data=cov(Ts), bls=bls_, gps=gps, cmd=' '.join(sys.argv))
            #NOTE: Ts gets over written below and so is the original covariance matrix. 
            #The cov(Ts) above is the fully diagonalized covariance matrix from before.
            continue
            
        if PLOT:
            #capo.arp.waterfall(cov(Ts), mode='log', mx=-1, drng=4);p.colorbar(shrink=.5)
            #p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ts)*scalar, mode='log', mx=8, drng=4);p.colorbar(shrink=.5)
            p.figure(10)
            #p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ns)*scalar, mode='log', mx=8, drng=4);p.colorbar(shrink=.5)
            p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ts)*scalar, mode='log', mx=8, drng=4);p.colorbar(shrink=.5)
            p.figure(11)
            #p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(avg_Cn*scalar, mode='log', mx=8,  drng=4); p.colorbar(shrink=.5)
            p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(avg_Cx*scalar, mode='log', mx=8,  drng=4); p.colorbar(shrink=.5)
            p.figure(12)
            p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(_Cxtot*scalar, mode='log', mx=8, drng=4);p.colorbar(shrink=.5)
            #p.subplot(PLT1,PLT2,cnt+2); capo.arp.waterfall(cov(Ns), mode='log', mx=0, drng=3)
            p.show()

#    p.show()
#    import IPython
#    IPython.embed()
#    exit()
    print bls_
    Ts = n.concatenate([T[bl] for bl in bls_], axis=1).T
    Ns = n.concatenate([N[bl] for bl in bls_], axis=1).T # this noise copy processed as if it were the data
    #need this var in the boot strapping
    temp_noise_var = n.average(n.array([T[bl] for bl in bls_]), axis=0).T
    pspecs,dspecs = [], []
    nspecs,n1specs,n2specs = [], [], []
    Cx,Cn = CoV(Ts, bls_, times), CoV(Ns, bls_, times)
    Cx_ = CoV(n.dot(_Cxtot,Ts), bls_, times)
    if PLOT:
        p.figure(1)
        capo.arp.waterfall(_Cxtot*scalar, mode='log', mx=8, drng=4);p.colorbar(shrink=.5)
        p.title('C : Junk Matrix.')
        p.figure(2)
        capo.arp.waterfall(Cx_.C*scalar, mode='log', mx=8, drng=4);p.colorbar(shrink=.5)
        p.title('C : after removing off diagonal covariances')
        p.figure(3)
        capo.arp.waterfall(Cx.C*scalar, mode='log', mx=8, drng=4);p.colorbar(shrink=.5)
        p.title('C : Original covariance matrix of the data.')
        p.show()
    
    # Cn1 is the noise diagonalized as if it were the signal, Cn2 is the noise with the signal diagonalization applied
    Cn1_,Cn2_ = CoV(n.dot(_Cntot,Ns), bls_, times), CoV(n.dot(_Cxtot,Ns), bls_, times)
    if opts.write2uv:
        write_2_uv( args, Cx_, afreqs, ending='boot%d'%boot)
        continue
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
            if not same_group(bli, blj, gps) is None: continue
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
    n.savez(outfile, kpl=kpl, pk=avg_1d, err=std_1d, scalar=scalar, times=n.array(times),freq=fq, pk_vs_t=avg_2d, err_vs_t=std_2d, temp_noise_var=temp_noise_var, nocov_vs_t=n.average(dspecs,axis=0), cov_matrix=_Cxtot, bls=bls_, gps=gps, cmd=' '.join(sys.argv))
    #n.savez(outfile, kpl=kpl, pk=avg_1d, err=std_1d, scalar=scalar, times=n.array(times),freq=fq, pk_vs_t=avg_2d, err_vs_t=std_2d, nocov_vs_t=n.average(dspecs,axis=0), cov_matrix=_Cxtot, bls=bls_, gps=gps, cmd=' '.join(sys.argv))
    #kpl = k parallels
    #pk = 
    #err = 
    #scalara = cosmological scalar
    #times = times used in estimating covariances
    #freq = middle frequency
    #pk_vs_t = 
    #err_vs_t = 
    #temp_noise_var = 
    #nocov_vs_t = 
    #cov_matrix = Covariance matrix of data. 
    #bls = baselines used in this bootstrap
    #gps = gps of baselines.
    #cmd = command run
if PLOT: p.show()

