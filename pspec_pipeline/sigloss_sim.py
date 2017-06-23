#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo
import capo.frf_conv as fringe
<<<<<<< HEAD
import glob, optparse, sys, random
=======
import glob
import sys
import optparse
import random
>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
import capo.zsa as zsa
from IPython import embed
import capo.oqe as oqe

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
o.add_option('--output', type='string', default='',
<<<<<<< HEAD
    help='output directory for pspec_boot files (default "")')
o.add_option('--noise_only',action='store_true',
    help='Replace data with noise.')
opts,args = o.parse_args(sys.argv[1:])

#Basic Parameters
=======
             help='output directory for pspec_boot files (default "")')
o.add_option('--noise_only', action='store_true',
             help='Replace data with noise.')
o.add_option('--rmbls', dest='rmbls', type='string',
             help=('List of baselines (ex:1_4,2_33) '
                   'to remove from the power spectrum analysis.'))
opts, args = o.parse_args(sys.argv[1:])

# Basic Parameters

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
random.seed(0)
#n.random.seed(1235813)
POL = opts.pol #'I'
LST_STATS = False
DELAY = False
NGPS = 5
INJECT_SIG = opts.inject_sig
PLOT = opts.plot
<<<<<<< HEAD

### FUNCTIONS ###

#function uses aa, ij, afreqs, inttime, POL
def frf(shape,loc=0,scale=1):
    shape = shape[1]*2,shape[0] #(2*times,freqs)
    dij = noise(shape,loc=loc,scale=scale)
    dij = oqe.noise(shape)
    #bins = fringe.gen_frbins(inttime)
    #frp, bins = fringe.aa_to_fr_profile(aa, ij, len(afreqs)/2, bins=bins)
    #timebins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(), fq0=aa.get_freqs()[len(afreqs)/2])
    #_,blconj,_ = zsa.grid2ij(aa.ant_layout)
    #if blconj[a.miriad.ij2bl(ij[0],ij[1])]: fir = {(ij[0],ij[1],POL):n.conj(firs)}
    #else: fir = {(ij[0],ij[1],POL):firs}
    wij = n.ones(shape,dtype=bool) #XXX flags are all true (times,freqs)
    #dij and wij are (times,freqs)
    _d,_w,_,_ = fringe.apply_frf(aa,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)
    _d = n.transpose(_d)
    _d = _d[:,shape[0]/4:shape[0]/2+shape[0]/4]
    return _d

def get_data(filenames, antstr, polstr, verbose=False):
=======
cov_reg_level = 0

try:
    rmbls = []
    rmbls_list = opts.rmbls.split(',')
    for bl in rmbls_list:
        i, j = bl.split('_')
        rmbls.append(a.miriad.ij2bl(int(i), int(j)))
    print 'Removing baselines:', rmbls
    # rmbls = map(int, opts.rmbls.split(','))
except:
    rmbls = []

#  FUNCTIONS #

# function uses aa, ij, afreqs, inttime, POL


def frf(shape, loc=0, scale=1):
    """Generate White Noise and apply FRF."""
    shape = shape[1]*4, shape[0]  # (4*times,freqs)
    dij = noise(shape, loc=loc, scale=scale)
    # dij = n.ones(shape)
    # bins = fringe.gen_frbins(inttime)
    # frp, bins = fringe.aa_to_fr_profile(aa, ij, len(afreqs)/2, bins=bins)
    # timebins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(),
    #                                      fq0=aa.get_freqs()[len(afreqs)/2])
    # _,blconj,_ = zsa.grid2ij(aa.ant_layout)
    # if blconj[a.miriad.ij2bl(ij[0],ij[1])]:
    #    fir = {(ij[0],ij[1],POL):n.conj(firs)}
    # else:
    #    fir = {(ij[0],ij[1],POL):firs}
    wij = n.ones(shape, dtype=bool)  # XXX flags are all true (times,freqs)
    # dij and wij are (times,freqs)
    _d, _w, _, _ = fringe.apply_frf(aa, dij, wij, ij[0], ij[1],
                                    pol=POL, bins=bins, firs=fir)
    # _d,_w,_,_ = fringe.apply_frf(aa, _d, wij, ij[0] ,ij[1] ,
    #                               pol=POL,bins=bins,firs=fir)
    _d = n.transpose(_d)
    _d = _d[:, shape[0]/2. - shape[0]/4./2.: shape[0]/2. + shape[0]/4./2.]
    return _d


def get_data(filenames, antstr, polstr, rmbls, verbose=False):
    """Return data dictionary from files."""

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
    # XXX could have this only pull channels of interest to save memory
    lsts, dat, flg = [], {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
<<<<<<< HEAD
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            lst = uv['lst']
            if len(lsts) == 0 or lst != lsts[-1]: lsts.append(lst)
            bl = a.miriad.ij2bl(i,j)
            if not dat.has_key(bl): dat[bl],flg[bl] = [],[]
=======
        for (crd, t, (i, j)), d, f in uv.all(raw=True):
            bl = a.miriad.ij2bl(i, j)
            if bl in rmbls:
                continue
            lst = uv['lst']
            if len(lsts) == 0 or lst != lsts[-1]:
                lsts.append(lst)
            # bl = a.miriad.ij2bl(i,j)
            if bl not in dat.keys():
                dat[bl], flg[bl] = [], []
>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
            dat[bl].append(d)
            flg[bl].append(f)
            #if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            #pol = a.miriad.pol2str[uv['pol']]
            #if not dat[bl].has_key(pol):
            #    dat[bl][pol],flg[bl][pol] = [],[]
            #dat[bl][pol].append(d)
            #flg[bl][pol].append(f)
    order_lst = n.argsort(lsts) #sometimes data is not in LST order!
    lsts = n.array(lsts)[order_lst]
    for bl in dat:
        dat[bl] = n.array(dat[bl])[order_lst]
        flg[bl] = n.array(flg[bl])[order_lst]
    return lsts, dat, flg

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (n.dot(X, X.T.conj()) / fact).squeeze()

<<<<<<< HEAD
def noise(size,loc=0,scale=1): #loc is mean, scale is stdev (sqrt(var))
    #sig = 1./n.sqrt(2)
    #return n.random.normal(scale=sig, size=size) + 1j*n.random.normal(scale=sig, size=size)
    return (n.random.normal(size=size,scale=scale) * n.exp(1j*n.random.uniform(0,2*n.pi,size=size))) + loc
=======

def noise(size, loc=0, scale=1):  # loc is mean, scale is stdev (sqrt(var))
    """Generate Complex Gaussian Noise."""
    sig = 1./n.sqrt(2)
    return (n.random.normal(scale=sig, size=size) +
            1j*n.random.normal(scale=sig,  size=size))
    # return ((n.random.normal(size=size,scale=scale) *
    #          n.exp(1j*n.random.uniform(0,2*n.pi,size=size))) + loc)


>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a

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

#SEP = 'sep0,2'

#Read even&odd datasets
dsets = {
    'even': [x for x in args if 'even' in x],
    'odd': [x for x in args if 'odd' in x]
    #'even': glob.glob('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/lstbin_manybls/even/sep0,2/*I.uvGAL'),
    #'odd': glob.glob('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/lstbin_manybls/odd/sep0,2/*I.uvGAL')
    #'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even_xtalk_removed/'+SEP+'/newfrf/*242.[3456]*uvGL'),
    #'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd_xtalk_removed/'+SEP+'/newfrf/*243.[3456]*uvGL'),
}
#for i in xrange(10): dsets[i] = glob.glob('lstbinX%d/%s/lst.24562[45]*.[3456]*.uvAL'%(i,SEP))

#Get uv file info
WINDOW = opts.window
uv = a.miriad.UV(dsets.values()[0][0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
#inttime = uv['inttime'] #* 4 # XXX hack for *E files that have inttime set incorrectly
(uvw,t1,(i,j)),d = uv.read()
(uvw,t2,(i,j)),d = uv.read()
while t1 == t2: (uvw,t2,(i,j)),d = uv.read()
inttime = (t2-t1)* (3600*24)
del(uv)

afreqs = freqs.take(chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

#aa = a.cal.get_aa(opts.cal, n.array([.150]))
aa = a.cal.get_aa(opts.cal, afreqs)
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
jy2T = capo.pspec.jy2T(afreqs)
window = a.dsp.gen_window(nchan, WINDOW)
if not WINDOW == 'none': window.shape=(nchan,1)
sep2ij, blconj, bl2sep = zsa.grid2ij(aa.ant_layout)

#B = sdf * afreqs.size / capo.pfb.NOISE_EQUIV_BW[WINDOW] # this is wrong if we aren't inverting
# the window post delay transform (or at least dividing out by the gain of the window)
# For windowed data, the FFT divides out by the full bandwidth, B, which is
# then squared.  Proper normalization is to multiply by B**2 / (B / NoiseEqBand) = B * NoiseEqBand
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

#Acquire data
#antstr = '41_49,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_47,3_25,1_58,35_61,42_63,2_33'
antstr = 'cross'
lsts,data,flgs = {},{},{}
days = dsets.keys()
for k in days:
<<<<<<< HEAD
    lsts[k],data[k],flgs[k] = get_data(dsets[k], antstr=antstr, polstr=POL, verbose=True)
    #data has keys 'even' and 'odd'
    #inside that are baseline keys
    #inside that has shape (#lsts, #freqs)
=======
    lsts[k], data[k], flgs[k] = get_data(dsets[k], antstr=antstr, polstr=POL,
                                         rmbls=rmbls, verbose=True)
# data has keys 'even' and 'odd'
# inside that are baseline keys
# inside that has shape (#lsts, #freqs)

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a

#Get some statistics
if LST_STATS:
    #collect some metadata from the lst binning process
    cnt, var = {}, {}
    times = []
    for filename in dsets.values()[0]:
        print 'Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, '64_49', POL)
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            if len(times) == 0 or times[-1] != uv['lst']: times.append(uv['lst'])
            bl = '%d,%d,%d' % (i,j,uv['pol'])
            cnt[bl] = cnt.get(bl, []) + [uv['cnt']]
            var[bl] = var.get(bl, []) + [uv['var']]
    cnt = n.array(cnt.values()[0]) #all baselines should be the same
    var = n.array(var.values()[0]) #all baselines should be the same
    times = n.array(times)
else: cnt,var = n.ones_like(lsts.values()[0]), n.ones_like(lsts.values()[0])


#Align data in LST (even/odd data might have a different number of LSTs)
lstr, order = {}, {}
lstres = 0.001
for k in lsts: #orders LSTs to find overlap
    order[k] = n.argsort(lsts[k])
    lstr[k] = n.around(lsts[k][order[k]] / lstres) * lstres
lsts_final = None
for i,k1 in enumerate(lstr.keys()):
    for k2 in lstr.keys()[i:]:
        if lsts_final is None: lsts_final = n.intersect1d(lstr[k1],lstr[k2]) #XXX LSTs much match exactly
        else: lsts_final = n.intersect1d(lsts_final,lstr[k2])
inds = {}
for k in lstr: #selects correct LSTs from data
    inds[k] = order[k].take(lstr[k].searchsorted(lsts_final))
lsts = lsts[lsts.keys()[0]][inds[lsts.keys()[0]]]
for k in days:
    for bl in data[k]:
        data[k][bl],flgs[k][bl] = data[k][bl][inds[k]],flgs[k][bl][inds[k]]

"""# XXX found a bug in this original code (lsts['even'] and lsts['odd'] are different!)
lstmax = max([lsts[k][0] for k in days]) #the larger of the initial lsts
for k in days:
    #print k
    for i in xrange(len(lsts[k])):
        #allow for small numerical differences (which shouldn't exist!)
        if lsts[k][i] >= lstmax - .001: break
    lsts[k] = lsts[k][i:]
    for bl in data[k]:
        data[k][bl],flgs[k][bl] = data[k][bl][i:],flgs[k][bl][i:]
j = min([len(lsts[k]) for k in days])
for k in days:
    lsts[k] = lsts[k][:j]
    for bl in data[k]:
        data[k][bl],flgs[k][bl] = n.array(data[k][bl][:j]),n.array(flgs[k][bl][:j])
lsts = lsts.values()[0] #same set of LST values for both even/odd data
"""
daykey = data.keys()[0]
blkey = data[daykey].keys()[0]
ij = a.miriad.bl2ij(blkey)
<<<<<<< HEAD
if blconj[blkey]: ij = (ij[1], ij[0])
#ij = (64, 49) 
=======
if blconj[blkey]:
    ij = (ij[1], ij[0])

sep = bl2sep[blkey]
c = 0
while c != -1:
    ij = map(int, sep2ij[sep].split(',')[c].split('_'))
    blkey = a.miriad.ij2bl(*ij)
    if blconj[blkey]:
        c += 1
    else:
        break
# ij = (64, 49)

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a

#Prep FRF Stuff
bins = fringe.gen_frbins(inttime)
frp, bins = fringe.aa_to_fr_profile(aa, ij, len(afreqs)/2, bins=bins)
<<<<<<< HEAD
timebins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(), fq0=aa.get_freqs()[len(afreqs)/2-4])
_,blconj,_ = zsa.grid2ij(aa.ant_layout)
#if blconj[a.miriad.ij2bl(ij[0],ij[1])]: fir = {(ij[0],ij[1],POL):n.conj(firs)}
#else: fir = {(ij[0],ij[1],POL):firs}
fir = {(ij[0],ij[1],POL):firs}

#Extract frequency range of data 
=======
# , alietal=True)
timebins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(),
                                    fq0=aa.get_freqs()[len(afreqs)/2])
#                                   , alietal=True)#, maxfr=1.3e-3,
#                                   frwidth=2.0)

_, blconj, _ = zsa.grid2ij(aa.ant_layout)
# if blconj[a.miriad.ij2bl(ij[0],ij[1])]:
#    fir = {(ij[0],ij[1],POL):n.conj(firs)}
# else:
#    fir = {(ij[0],ij[1],POL):firs}
fir = {(ij[0], ij[1], POL): firs}

# Extract frequency range of data

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
xi = {}
f = {}
#NOISE = frf((len(chans),len(lsts)),loc=0,scale=1) #same noise on each bl
for k in days:
    xi[k] = {}
    f[k] = {}
    for bl in data[k]:
        d = data[k][bl][:,chans] * jy2T
        flg = flgs[k][bl][:,chans]
        if conj[a.miriad.bl2ij(bl)]: d = n.conj(d) #conjugate if necessary
        shape = d.shape #(times,freqs)
        if opts.noise_only:
            xi[k][bl] =  frf((len(chans),len(lsts)),loc=0,scale=1) #diff noise for each bl
        else:
             xi[k][bl] = n.transpose(d, [1,0]) #swap time and freq axes
        f[k][bl] = n.transpose(flg, [1,0])
bls_master = xi.values()[0].keys()
#bls_master = bls_master[:5]
nbls = len(bls_master)
print 'Baselines:', nbls
<<<<<<< HEAD
=======
print 'N times:', n.shape(xi[k][bl])[-1]

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a

#Bootstrapping
for boot in xrange(opts.nboot):

<<<<<<< HEAD
    print '%d / %d' % (boot+1,opts.nboot)   
=======
    print '%d / %d' % (boot+1, opts.nboot)

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a

    ### Calculate pC just based on the data/simulation noise (no eor injection) ###
    print 'Getting pCv'
    x = {}
    for k in days:
        x[k] = {}
        for bl in xi[k]: x[k][bl] = xi[k][bl].copy()
    #XXX x is JUST data

    Q = [get_Q(i,nchan) for i in xrange(nchan)] #get Q matrix (does FT from freq to delay)
    I,_I,_Ix = {},{},{}
    C,_C,_Cx = {},{},{}
    for k in days:
        I[k],_I[k],_Ix[k] = {},{},{}
        C[k],_C[k],_Cx[k] = {},{},{}
        for bl in bls_master:
<<<<<<< HEAD
            C[k][bl] = cov(x[k][bl]) 
=======
            C[k][bl] = cov(x[k][bl])
            C[k][bl] += cov_reg_level * n.identity(C[k][bl].shape[0])
>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
            I[k][bl] = n.identity(C[k][bl].shape[0])
            U,S,V = n.linalg.svd(C[k][bl].conj()) #singular value decomposition
            _C[k][bl] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
            _I[k][bl] = n.identity(_C[k][bl].shape[0])
            _Cx[k][bl] = n.dot(_C[k][bl], x[k][bl]) # XXX
            _Ix[k][bl] = x[k][bl] # XXX
            if PLOT and True:
                print a.miriad.bl2ij(bl), k
                p.subplot(311); capo.arp.waterfall(x[k][bl], mode='real')
                p.colorbar()
                p.title('Data x')
                p.subplot(323); capo.arp.waterfall(C[k][bl])
                p.colorbar()
                p.title('C')
                p.subplot(324); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
                p.subplot(313); capo.arp.waterfall(_Cx[k][bl], mode='real')
                p.colorbar()
                p.title('C^-1 x')
                p.suptitle('%d_%d'%a.miriad.bl2ij(bl)+' '+k)
                p.tight_layout()
                p.show()
<<<<<<< HEAD
        
    #Make boots
=======

    # Make boots
>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
    bls = bls_master[:]
    if True: #shuffle and group baselines for bootstrapping
        random.shuffle(bls)
        gps = [bls[i::NGPS] for i in range(NGPS)]
        gps = [[random.choice(gp) for bl in gp] for gp in gps]
    else: gps = [bls[i::NGPS] for i in range(NGPS)]
    bls = [bl for gp in gps for bl in gp]
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
            for ch in xrange(nchan): # XXX this loop makes computation go as nchan^3
                _IsumQ[k][i][ch] = n.dot(_Isum[k][i], Q[ch])
                _CsumQ[k][i][ch] = n.dot(_Csum[k][i], Q[ch])
        if PLOT:
            NGPS = len(gps)
            _Csumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            _Isumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            for i in xrange(len(gps)): _Isumk[i,:,i,:] = _Isum[k][i]
            _Isumk.shape = (NGPS*nchan, NGPS*nchan)
            for i in xrange(len(gps)): _Csumk[i,:,i,:] = _Csum[k][i]
            _Csumk.shape = (NGPS*nchan, NGPS*nchan)
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
                        if not Q_Iz[k2].has_key(bl2): Q_Iz[k2][bl2] = [n.dot(Q[i], _Iz[k2][bl2]) for i in xrange(nchan)]
                        if not Q_Cz[k2].has_key(bl2): Q_Cz[k2][bl2] = [n.dot(Q[i], _Cz[k2][bl2]) for i in xrange(nchan)]
                        _qI = n.array([_Iz[k1][bl1].conj() * Q_Iz[k2][bl2][i] for i in xrange(nchan)])
                        qI += n.sum(_qI, axis=1)
                        _qC = n.array([_Cz[k1][bl1].conj() * Q_Cz[k2][bl2][i] for i in xrange(nchan)]) #C^-1 Q C^-1
                        qC += n.sum(_qC, axis=1)
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

    print "   Getting M"
<<<<<<< HEAD
    order = n.array([10,11,9,12,8,20,0,13,7,14,6,15,5,16,4,17,3,18,2,19,1])
    iorder = n.argsort(order)
    FC_o = n.take(n.take(FC,order, axis=0), order, axis=1)
    L_o = n.linalg.cholesky(FC_o)
    U,S,V = n.linalg.svd(L_o.conj())
    MC_o = n.dot(n.transpose(V), n.dot(n.diag(1./S), n.transpose(U)))
    MC = n.take(n.take(MC_o,iorder, axis=0), iorder, axis=1)
    MI  = n.identity(nchan, dtype=n.complex128)
=======
    # order = n.array([10, 11, 9, 12, 8, 20, 0, 13, 7,
    #                  14, 6, 15, 5, 16, 4, 17, 3, 18, 2, 19, 1])
    # iorder = n.argsort(order)
    # FC_o = n.take(n.take(FC, order, axis=0), order, axis=1)
    # L_o = n.linalg.cholesky(FC_o)
    # U, S, V = n.linalg.svd(L_o.conj())
    # MC_o = n.dot(n.transpose(V), n.dot(n.diag(1./S), n.transpose(U)))
    # MC = n.take(n.take(MC_o, iorder, axis=0), iorder, axis=1)
    MC = n.identity(nchan, dtype=n.complex128)
    MI = n.identity(nchan, dtype=n.complex128)

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a

    print "   Getting W"
    WI = n.dot(MI, FI)
    norm  = WI.sum(axis=-1); norm.shape += (1,)
    MI /= norm; WI = n.dot(MI, FI)
    WC = n.dot(MC, FC)
    norm  = WC.sum(axis=-1); norm.shape += (1,)
    MC /= norm; WC = n.dot(MC, FC)

    print '   Generating ps'
    #if opts.noise_only: scalar = 1
    pC = n.dot(MC, qC) * scalar
    pI = n.dot(MI, qI) * scalar 

    #XXX Overwriting to new variables
    pCv = pC.copy()
<<<<<<< HEAD
    pIv = pI

    
    ### Loop to calculate pC of (data/noise+eor) and pI of eor ###
=======
    pIv = pI.copy()

    # Loop to calculate pC of (data/noise+eor) and pI of eor #

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
    print 'Getting pCr and pIe'

    if INJECT_SIG > 0.: #Create a fake EoR signal to inject
        print 'INJECTING SIMULATED SIGNAL'
<<<<<<< HEAD
        eor = frf((shape[1],shape[0]),loc=0,scale=1) * INJECT_SIG #create FRF-ered noise
=======
        eor = frf((shape[1], shape[0]), loc=0, scale=1) * INJECT_SIG
        # create FRF-ered noise
        # eor = frf((shape[1],shape[0]), loc=0, scale=1) * (boot+1)
        # create FRF-ered noise
        # eor  = noise(size=(shape[1],shape[0])) * INJECT_SIG

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
        x = {}
        for k in days:
            x[k] = {}
            for bl in xi[k]: x[k][bl] = xi[k][bl].copy() + eor.copy() #add injected signal to data
        if False and PLOT:
            p.subplot(211); capo.arp.waterfall(eor1, mode='real'); p.colorbar()
            p.subplot(212); capo.arp.waterfall(eor2, mode='real'); p.colorbar(); p.show()
    Q = [get_Q(i,nchan) for i in xrange(nchan)] #get Q matrix (does FT from freq to delay)
    I,_I,_Ix = {},{},{}
    C,_C,_Cx = {},{},{}
    for k in days:
        I[k],_I[k],_Ix[k] = {},{},{}
        C[k],_C[k],_Cx[k] = {},{},{}
        for bl in bls_master:
<<<<<<< HEAD
            C[k][bl] = cov(x[k][bl]) 
=======
            C[k][bl] = cov(x[k][bl])
            C[k][bl] += cov_reg_level * n.identity(C[k][bl].shape[0])

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
            I[k][bl] = n.identity(C[k][bl].shape[0])
            U,S,V = n.linalg.svd(C[k][bl].conj()) #singular value decomposition
            _C[k][bl] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
            _I[k][bl] = n.identity(_C[k][bl].shape[0])
            #_C[k][bl] = n.identity(_C[k][bl].shape[0]) #XXX overwriting C with I
            _Cx[k][bl] = n.dot(_C[k][bl], x[k][bl]) # XXX x = v+e
            _Ix[k][bl] = eor.copy() # XXX only e!
            if PLOT and True:
                print a.miriad.bl2ij(bl), k
                p.subplot(311); capo.arp.waterfall(x[k][bl], mode='real')
                p.colorbar()
                p.title('Data x')
                p.subplot(323); capo.arp.waterfall(C[k][bl])
                p.colorbar()
                p.title('C')
                p.subplot(324); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
                p.subplot(313); capo.arp.waterfall(_Cx[k][bl], mode='real')
                p.colorbar()
                p.title('C^-1 x')
                p.suptitle('%d_%d'%a.miriad.bl2ij(bl)+' '+k)
                p.tight_layout()
                p.show()
<<<<<<< HEAD
        
    #Make boots
=======

    # Make boots

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
    bls = bls_master[:]
    #if True: #XXX GPS already defined the first time
    #shuffle and group baselines for bootstrapping
    #    random.shuffle(bls)
    #    gps = [bls[i::NGPS] for i in range(NGPS)]
    #    gps = [[random.choice(gp) for bl in gp] for gp in gps]
    #else: gps = [bls[i::NGPS] for i in range(NGPS)]
    bls = [bl for gp in gps for bl in gp]
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
            for ch in xrange(nchan): # XXX this loop makes computation go as nchan^3
                _IsumQ[k][i][ch] = n.dot(_Isum[k][i], Q[ch])
                _CsumQ[k][i][ch] = n.dot(_Csum[k][i], Q[ch])
        if PLOT:
            NGPS = len(gps)
            _Csumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            _Isumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            for i in xrange(len(gps)): _Isumk[i,:,i,:] = _Isum[k][i]
            _Isumk.shape = (NGPS*nchan, NGPS*nchan)
            for i in xrange(len(gps)): _Csumk[i,:,i,:] = _Csum[k][i]
            _Csumk.shape = (NGPS*nchan, NGPS*nchan)
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
                        if not Q_Iz[k2].has_key(bl2): Q_Iz[k2][bl2] = [n.dot(Q[i], _Iz[k2][bl2]) for i in xrange(nchan)]
                        if not Q_Cz[k2].has_key(bl2): Q_Cz[k2][bl2] = [n.dot(Q[i], _Cz[k2][bl2]) for i in xrange(nchan)]
                        _qI = n.array([_Iz[k1][bl1].conj() * Q_Iz[k2][bl2][i] for i in xrange(nchan)])
                        qI += n.sum(_qI, axis=1)
                        _qC = n.array([_Cz[k1][bl1].conj() * Q_Cz[k2][bl2][i] for i in xrange(nchan)]) #C^-1 Q C^-1
                        qC += n.sum(_qC, axis=1)
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

    print "   Getting M"
<<<<<<< HEAD
    #Cholesky decomposition
    order = n.array([10,11,9,12,8,20,0,13,7,14,6,15,5,16,4,17,3,18,2,19,1])
    iorder = n.argsort(order)
    FC_o = n.take(n.take(FC,order, axis=0), order, axis=1)
    L_o = n.linalg.cholesky(FC_o)
    U,S,V = n.linalg.svd(L_o.conj())
    MC_o = n.dot(n.transpose(V), n.dot(n.diag(1./S), n.transpose(U)))
    MC = n.take(n.take(MC_o,iorder, axis=0), iorder, axis=1)
    MI  = n.identity(nchan, dtype=n.complex128)
    
=======
    # Cholesky decomposition
    # order = n.array([10, 11, 9, 12, 8, 20, 0, 13, 7, 14, 6,
    #                  15, 5, 16, 4, 17, 3, 18, 2, 19, 1])
    # iorder = n.argsort(order)
    # FC_o = n.take(n.take(FC, order, axis=0), order, axis=1)
    # L_o = n.linalg.cholesky(FC_o)
    # U, S, V = n.linalg.svd(L_o.conj())
    # MC_o = n.dot(n.transpose(V), n.dot(n.diag(1./S), n.transpose(U)))
    # MC = n.take(n.take(MC_o, iorder, axis=0), iorder, axis=1)
    MC = n.identity(nchan, dtype=n.complex128)
    MI = n.identity(nchan, dtype=n.complex128)


>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
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
    #if opts.noise_only: scalar = 1
    pC = n.dot(MC, qC) * scalar
    pI = n.dot(MI, qI) * scalar 

    if PLOT:
        p.subplot(411); capo.arp.waterfall(qC, mode='real'); p.colorbar(shrink=.5)
        p.subplot(412); capo.arp.waterfall(pC, mode='real'); p.colorbar(shrink=.5)
        p.subplot(413); capo.arp.waterfall(qI, mode='real'); p.colorbar(shrink=.5)
        p.subplot(414); capo.arp.waterfall(pI, mode='real'); p.colorbar(shrink=.5)
        p.show()

    #XXX Overwriting to new variables
    pCr = pC
    pIe = pI
    #XXX Final variables
    pI = pIe
<<<<<<< HEAD
    pC = pCr - pCv
    
    print 'pI=', n.average(pI.real), 'pC=', n.average(pC.real), 'pI/pC=', n.average(pI.real)/n.average(pC.real)
   
=======
    pC = pCr  # - pCv

    print 'pI=', n.average(pI.real), 'pC=', n.average(pC.real),
    print 'pI/pC=', n.average(pI.real)/n.average(pC.real)
    # embed()

>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a
    if PLOT:
        p.plot(kpl, n.average(pC.real, axis=1), 'b.-')
        p.plot(kpl, n.average(pI.real, axis=1), 'k.-')
        p.show()

    print '   Writing pspec_bootsigloss%04d.npz' % boot

    if len(opts.output) > 0: outpath = opts.output+'/pspec_bootsigloss%04d.npz' % boot
    else: outpath = 'pspec_bootsigloss%04d.npz' % boot
    n.savez(outpath, kpl=kpl, scalar=scalar, times=n.array(lsts),
<<<<<<< HEAD
        pk_vs_t=pC, pCv = pCv, err_vs_t=1./cnt, temp_noise_var=var, nocov_vs_t=pI, freq=fq,
        cmd=' '.join(sys.argv))

=======
            pk_vs_t=n.real(pC), pCv=n.real(pCv), err_vs_t=1./cnt,
            temp_noise_var=var, nocov_vs_t=n.real(pI),
            freq=fq, pIv=n.real(pIv), cmd=' '.join(sys.argv))
>>>>>>> 37b76786b503ec7ac46bc3a33ba1824e8074254a

