#! /usr/bin/env python

from matplotlib import pylab as p
import aipy as a, numpy as n, pylab as p, capo, capo.frf_conv as fringe
import glob, optparse, sys, random
from scipy.linalg import fractional_matrix_power
import capo.zsa as zsa
import capo.oqe as oqe

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
o.add_option('-b', '--nboot', type='int', default=20,
    help='Number of bootstraps.  Default is 20.')
o.add_option('--plot', action='store_true',
    help='Generate plots')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
o.add_option('--sep', default='sep0,1', action='store',
    help='Which separation directory to use for signal loss data.')
o.add_option('--loss', action='store', 
    help='Uses signal loss mode to measure the signal loss. Uses default data in my path. Give it the path to the simulated signal data.')
o.add_option('--level', type='float', default=-1.0,
    help='Scalar by which to multiply the default signal level for simulation runs.')
o.add_option('--rmbls', dest='rmbls',type='string', 
    help='List of baselines (ex:1_4,2_33) to remove from the power spectrum analysis.')
#o.add_option('--noise', type='float', default=0.0,
#    help='amplitude of noise to inject into data')
#o.add_option('--filter_noise',action='store_true',
#    help='Apply fringe rate filter to injected noise')
o.add_option('--noise_only',action='store_true',
    help='Instead of injecting noise, Replace data with noise')
o.add_option('--output', type='string', default='',
    help='Output directory for pspec_boot files (default "")')

opts,args = o.parse_args(sys.argv[1:])

#Basic parameters
random.seed(0)
POL = opts.pol #'xx'
LST_STATS = False
DELAY = False
NGPS = 5 #number of groups to break the random sampled bls into
INJECT_SIG = 0.
SAMPLE_WITH_REPLACEMENT = False
#NOISE = opts.noise
PLOT = opts.plot

#Remove baselines if specified
try:
    rmbls = []
    rmbls_list = opts.rmbls.split(',')
    for bl in rmbls_list:
        i,j = bl.split('_')
        rmbls.append(a.miriad.ij2bl(int(i),int(j)))   
    print 'Removing baselines:',rmbls
    #rmbls = map(int, opts.rmbls.split(','))
except:
    rmbls = []

#Set up signal loss mode if specified
if opts.loss:
    if opts.level >= 0.0:
        INJECT_SIG = opts.level
        print 'Running in signal loss mode, with an injection signal of %s*default level'%(opts.level)
    else:
        print 'Exiting. If in signal loss mode, need a signal level to input.'
        exit()


### FUNCTIONS ###

def frf(shape,loc=0,scale=1): #FRF NOISE
    shape = shape[1]*2,shape[0] #(2*times,freqs)
    dij = noise(shape,loc=loc,scale=scale)
    #dij = oqe.noise(size=shape)
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

def frf_data(data): #FRF FULL LENGTH DATA
    aa2 = a.cal.get_aa(opts.cal, freqs) #All freqs for data instead of subset
    bins = fringe.gen_frbins(inttime)
    frp, bins = fringe.aa_to_fr_profile(aa2, ij, len(freqs)/2, bins=bins)
    timebins, firs = fringe.frp_to_firs(frp, bins, aa2.get_freqs(), fq0=aa2.get_freqs()[len(freqs)/2])
    _,blconj,_ = zsa.grid2ij(aa2.ant_layout)
    if blconj[a.miriad.ij2bl(ij[0],ij[1])]: fir = {(ij[0],ij[1],POL):n.conj(firs)}
    else: fir = {(ij[0],ij[1],POL):firs}
    dij = n.transpose(data)
    wij = n.ones(dij.shape,dtype=bool) #XXX flags are all true (times,freqs)
    #dij and wij are (times,freqs)
    _d,_w,_,_ = fringe.apply_frf(aa2,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)
    _d = n.transpose(_d)
    return _d

def get_data(filenames, antstr, polstr, rmbls, verbose=False):
    # XXX could have this only pull channels of interest to save memory
    lsts, dat, flg = [], {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            bl = a.miriad.ij2bl(i,j)
            #print bl,uv['pol']
            if bl in rmbls: continue
            lst = uv['lst']
            if len(lsts) == 0 or lst != lsts[-1]: lsts.append(lst)
            if not dat.has_key(bl): dat[bl],flg[bl] = [],[]
            dat[bl].append(d)
            flg[bl].append(f)
            #if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            #pol = a.miriad.pol2str[uv['pol']]
            #if not dat[bl].has_key(pol):
            #    dat[bl][pol],flg[bl][pol] = [],[]
            #dat[bl][pol].append(d)
            #flg[bl][pol].append(f)
    #print len(dat[bl]),len(lsts)
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
    fact = float(N - 1) #normalization
    return (n.dot(X, X.T.conj()) / fact).squeeze()

def noise(size,loc=0,scale=1): #loc is mean, scale is stdev (sqrt(var))
    #return n.random.normal(size=size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=size))
    return (n.random.normal(size=size,scale=scale) * n.exp(1j*n.random.uniform(0,2*n.pi,size=size))) + loc

def get_Q(mode, n_k): #encodes the fourier transform from freq to delay
    if not DELAY:
        _m = n.zeros((n_k,), dtype=n.complex)
        _m[mode] = 1. #delta function at specific delay mode
        m = n.fft.fft(n.fft.ifftshift(_m)) * a.dsp.gen_window(nchan, WINDOW) #FFT it to go to freq
        Q = n.einsum('i,j', m, m.conj()) #dot it with its conjugate
        return Q
    else:
        # XXX need to have this depend on window
        Q = n.zeros_like(C)
        Q[mode,mode] = 1
        return Q


SEP = opts.sep #XXX used only for signal loss paths?

#Read even&odd data
if 'even' in args[0] or 'odd' in args[0]:
    dsets = {
    'even': [x for x in args if 'even' in x],
    'odd' : [x for x in args if 'odd' in x]
    }
else:
    jds = []
    for arg in args:
        jds.append(arg.split('zen')[1].split('.')[1])
    lenjds = len(n.unique(jds)) #length of unique JDs
    if lenjds == 1: dsets = {'even': args, 'odd': args}
    elif lenjds == 2:
        dsets = {
        'even': [x for x in args if n.unique(jds)[0] in x],
        'odd': [x for x in args if n.unique(jds)[1] in x]
        }
print dsets

"""
dsets = {
    #'only': glob.glob('sep0,1/*242.[3456]*uvL'),
    #'vissim1' : glob.glob('/home/cacheng/capo/ctc/tables/203files/pspec_Jy.uv'), 
    #'vissim2' : glob.glob('/home/cacheng/capo/ctc/tables/203files/pspec_Jy.uv'),
    #vissim3' : glob.glob('/home/cacheng/capo/ctc/tables/203files/gsm_K_2.uv'),
    #'vissim4' : glob.glob('/home/cacheng/capo/ctc/tables/203files/gsm_K_2.uv'),
    'even': [x for x in args if 'even' in x],
    'odd' : [x for x in args if 'odd' in x]
}
"""
print dsets

#Read signal loss data
if opts.loss:
    dsets = {
    'even': glob.glob('/data2/ali_et_al_2015_apj_data/even/*242.[3456]*uvALG'),
    'odd' : glob.glob('/data2/ali_et_al_2015_apj_data/odd/*243.[3456]*uvALG'),
}

#Get uv file info
WINDOW = opts.window
uv = a.miriad.UV(dsets.values()[0][0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])

#manually calculate the inttime based off time difference in data
#(uvw,t1,(i,j)),d = uv.read()
#(uvw,t2,(i,j)),d = uv.read()
#while t1 == t2: (uvw,t2,(i,j)),d = uv.read()
#inttime = (t2-t1)* (3600*24)
#del(uv)

inttime = uv['inttime']

afreqs = freqs.take(chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

aa = a.cal.get_aa(opts.cal, afreqs)
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
sep2ij, blconj, bl2sep = capo.zsa.grid2ij(aa.ant_layout)
jy2T = capo.pspec.jy2T(afreqs)
window = a.dsp.gen_window(nchan, WINDOW)
if not WINDOW == 'none': window.shape=(nchan,1)

#B = sdf * afreqs.size / capo.pfb.NOISE_EQUIV_BW[WINDOW] #this is wrong if we aren't inverting
# the window post delay transform (or at least dividing out by the gain of the window)
# For windowed data, the FFT divides out by the full bandwidth, B, which is
# then squared.  Proper normalization is to multiply by B**2 / (B / NoiseEqBand) = B * NoiseEqBand
# XXX NEED TO FIGURE OUT BW NORMALIZATION
B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW] #proper normalization
etas = n.fft.fftshift(capo.pspec.f2eta(afreqs)) #create etas (fourier dual to frequency)
#etas = capo.pspec.f2eta(afreqs) #create etas (fourier dual to frequency)
kpl = etas * capo.pspec.dk_deta(z) 
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
    lsts[k],data[k],flgs[k] = get_data(dsets[k], antstr=antstr, polstr=POL, rmbls=rmbls, verbose=True) 
    #data has keys 'even' and 'odd'
    #inside that are baseline keys
    #inside that has shape (#lsts, #freqs)

#Get some statistics
if LST_STATS:
    #collect some metadata from the lst binning process
    cnt, var = {}, {}
    for filename in dsets.values()[0]:
        print 'Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, '64_49', POL)
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            bl = '%d,%d,%d' % (i,j,uv['pol'])
            cnt[bl] = cnt.get(bl, []) + [uv['cnt']]
            var[bl] = var.get(bl, []) + [uv['var']]
    cnt = n.array(cnt.values()[0]) #all baselines should be the same
    var = n.array(var.values()[0]) #all baselines should be the same
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
if True:
#if False:
    #print [lsts[k][0] for k in days] #initial lst for even/odd data
    lstmax = max([lsts[k][0] for k in days]) #the larger of the  initial lsts
    for k in days:
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
else:
    for k in days:
        for bl in data[k]:
            data[k][bl], flgs[k][bl] = n.array(data[k][bl][:]), n.array(flgs[k][bl][:])
lsts = lsts.values()[0] #same set of LST values for both even/odd data
"""
daykey = data.keys()[0]
blkey = data[daykey].keys()[0]
ij = a.miriad.bl2ij(blkey)

#Prep FRF Stuff
bins = fringe.gen_frbins(inttime)
frp, bins = fringe.aa_to_fr_profile(aa, ij, len(afreqs)/2, bins=bins)
timebins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(), fq0=aa.get_freqs()[len(afreqs)/2])
_,blconj,_ = zsa.grid2ij(aa.ant_layout)
if blconj[a.miriad.ij2bl(ij[0],ij[1])]: fir = {(ij[0],ij[1],POL):n.conj(firs)}
else: fir = {(ij[0],ij[1],POL):firs}

"""
#FRF data #XXX
for k in days:
    for bl in data[k]:
        d = frf_data(n.transpose(data[k][bl]))
        d = n.transpose(d)
        data[k][bl] = d
"""

#Extract frequency range of data
x = {}
#if opts.noise_only: NOISE = frf((len(chans),len(lsts)),loc=0,scale=1) #same noise for all bl
for k in days:
    x = {}
    f = {}
    for k in days:
        x[k] = {}
        f[k] = {}
        for bl in data[k]:
            d = data[k][bl][:,chans] * jy2T
            flg = flgs[k][bl][:,chans]
            if conj[a.miriad.bl2ij(bl)]: d = n.conj(d) #conjugate if necessary
            if opts.noise_only:
                x[k][bl] = frf((len(chans),len(lsts)),loc=0,scale=1) #diff noise for each bl
            else:
                d = n.transpose(d) #now (freqs,times)
                x[k][bl] = d 
            f[k][bl] = n.transpose(flg, [1,0])
bls_master = x.values()[0].keys()
#bls_master = bls_master[:5] #XXX
nbls = len(bls_master)
print 'Baselines:', nbls

"""
#Inject fake EoR signal
if INJECT_SIG > 0.: 
    print 'INJECTING SIMULATED SIGNAL'
    eor_sets = {
    #    'only': glob.glob('sep0,1/*242.[3456]*uvL'),
#        'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even/'+SEP+'/*242.[3456]*uvALG_signalL'),
#        'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd/'+SEP+'/*243.[3456]*uvALG_signalL'),
        'even': glob.glob(opts.loss+'/even/*242.[3456]*uvALG_signalL'),
        'odd' : glob.glob(opts.loss+'/odd/*242.[3456]*uvALG_signalL'),
    }
    eorlsts,eordata,eorflgs = {},{},{}
    for k in days:
        eorlsts[k],eordata[k],eorflgs[k] = get_data(eor_sets[k], antstr=antstr, polstr=POL, verbose=True)
        #cut with same lst cut
        for bl in eordata[k]:
            eordata[k][bl], eorflgs[k][bl] = n.array(eordata[k][bl][:j]), n.array(eorflgs[k][bl][:j])
    eor = {}
    for k in days:
        eor[k] = {}
        for bl in eordata[k]:
            ed = eordata[k][bl][:,chans] * jy2T * INJECT_SIG
            if conj[bl]: ed = n.conj(ed)
            eor[k][bl] = n.transpose(ed, [1,0]) 

    for k in days:
        for bl in x[k]:
            x[k][bl] += eor[k][bl] 
    
    if PLOT:
        capo.arp.waterfall(x[k][bl], mode='real'); p.colorbar(); p.show()
""" 
"""
if NOISE > 0.: # Create a fake EoR signal to inject
    print 'INJECTING WHITE NOSIE at {0}Jy ...'.format(NOISE) ,
    sys.stdout.flush()
    noise_array = {}
    if opts.filter_noise:
        sep = bl2sep[bls_master[0]]
        c=0
        while c != -1:
            ij = map(int, sep2ij[sep].split(',')[c].split('_'))
            bl1 = a.miriad.ij2bl(*ij)
            if blconj[bl1]: c+=1
            else: break
        bl_ij= a.miriad.bl2ij(bl1)
        frbins = capo.frf_conv.gen_frbins(inttime)
        #frbins = n.arange( -.5/inttime+5e-5/2, .5/inttime,5e-5)
        #use channel 101 like in the frf_filter script
        frp, bins = capo.frf_conv.aa_to_fr_profile(aa, bl_ij,101, bins= frbins)
        timebins, firs = capo.frf_conv.frp_to_firs(frp, bins, aa.get_afreqs(), fq0=aa.get_afreqs()[101])

    for k in days:
        noise_array[k] = {}
        for bl in x[k]:

            ns = noise(x[k][bls_master[0]].shape) * NOISE
            wij = n.transpose(f[days[0]][bls_master[0]], [1,0]) #flags (time,freq)

            if opts.filter_noise:
                if blconj[a.miriad.ij2bl(ij[0],ij[1])]: fir = {(ij[0],ij[1],POL):n.conj(firs)} #conjugate fir if needed
                else: fir = {(ij[0],ij[1],POL):firs}
                dij,wij = n.transpose(ns, [1,0]),n.logical_not(wij)
                ns,_w,_,_ = fringe.apply_frf(aa,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)
                ns = n.transpose(ns, [1,0])


            if conj[bl]: ns = n.conj(ns)
            noise_array[k][bl] = n.transpose( ns.T * jy2T, [1,0])

    for k in days:
        for bl in x[k]:
            if opts.noise_only: x[k][bl] = n.copy( noise_array[k][bl])
            else: x[k][bl] += noise_array[k][bl]
            if PLOT and True:
               fig,axes = p.subplots(nrows=2,ncols=1)
               p.subplot(211); capo.arp.waterfall(x[k][bl], mode='real',mx=16,drng=32); #p.colorbar();
               p.ylabel('Data + Noise')
               p.subplot(212); im=  capo.arp.waterfall(noise_array[k][bl], mode='real',mx=16,drng=32); #p.colorbar(); p.show()
               p.ylabel('Noise')
               cbar_ax =fig.add_axes([.85,0.15,.05,.7])
               fig.subplots_adjust(right=.8)
               cbar=fig.colorbar(im,cax=cbar_ax)
               cbar.set_label('K')
               p.suptitle('%d_%d'%a.miriad.bl2ij(bl))
               p.show()
    print '[Done]'
sys.stdout.flush()
"""

#Power spectrum stuff
Q = [get_Q(i,nchan) for i in xrange(nchan)] #get Q matrix (does FT from freq to delay)

#Compute baseline auto-covariances and apply inverse to data
I,_I,_Ix = {},{},{}
C,_C,_Cx = {},{},{}
for k in days:
    I[k],_I[k],_Ix[k] = {},{},{}
    C[k],_C[k],_Cx[k] = {},{},{}
    for bl in bls_master:
        C[k][bl] = cov(x[k][bl])
        I[k][bl] = n.identity(C[k][bl].shape[0])
        #C[k][bl] = C[k][bl] + 1*I[k][bl] #C+NI noise
        U,S,V = n.linalg.svd(C[k][bl].conj()) #singular value decomposition
        _C[k][bl] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
        _I[k][bl] = n.identity(_C[k][bl].shape[0])
        _Cx[k][bl] = n.dot(_C[k][bl], x[k][bl])
        _Ix[k][bl] = x[k][bl].copy()
        if PLOT and True:
            #p.plot(S); p.show()
            p.subplot(311); capo.arp.waterfall(x[k][bl], mode='real')#,drng=0.5,mx=0)
            p.colorbar()
            p.title('Data x')
            p.subplot(323); capo.arp.waterfall(C[k][bl])
            p.colorbar()
            p.title('C')
            p.subplot(324); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
            p.subplot(313); capo.arp.waterfall(_Cx[k][bl], mode='real')#,drng=6000,mx=3000)
            p.colorbar()
            p.title('C^-1 x')
            p.suptitle('%d_%d'%a.miriad.bl2ij(bl)+' '+k)
            #p.figure(2); p.plot(n.diag(S))
            p.tight_layout()
            p.show()

#Make boots        
for boot in xrange(opts.nboot):
    print '%d / %d' % (boot+1,opts.nboot)
    bls = bls_master[:]
    if True: #shuffle and group baselines for bootstrapping
        if not SAMPLE_WITH_REPLACEMENT:
            random.shuffle(bls)
            #bls = bls[:-5] # XXX
        else: #sample with replacement (XXX could be biased)
            bls = [random.choice(bls) for bl in bls]
        gps = [bls[i::NGPS] for i in range(NGPS)]
        gps = [[random.choice(gp) for bl in gp] for gp in gps]
    else: #assign each baseline its own group
        gps = [bls[i::NGPS] for i in range(NGPS)]
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
                #XXX need to take fft of _Csum, _Isum here
            for ch in xrange(nchan): #XXX this loop makes computation go as nchan^3
                _IsumQ[k][i][ch] = n.dot(_Isum[k][i], Q[ch])
                _CsumQ[k][i][ch] = n.dot(_Csum[k][i], Q[ch]) #C^-1 Q
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
        for k2 in days[cnt1:]: #loop over even with even, even with odd, etc.
            if not Q_Iz.has_key(k2): Q_Iz[k2] = {}
            if not Q_Cz.has_key(k2): Q_Cz[k2] = {}
            for bl1 in _Cz[k1]:
                for bl2 in _Cz[k2]:
                    #if k1 == k2 and bl1 == bl2: continue #this results in a significant bias
                    if k1 == k2 or bl1 == bl2: continue 
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

    #other choices for M
    #U,S,V = n.linalg.svd(FC.conj())
    #_S = n.sqrt(1./S)
    # _S = 1./S
    # _S = n.ones_like(S)
    #MC = n.dot(n.transpose(V), n.dot(n.diag(_S), n.transpose(U)))
    #order = n.array([10,11,9,12,8,13,7,14,6,15,5,16,4,17,3,18,2,19,1,20,0])
    
    print "   Getting M"
    #Cholesky decomposition to get M
    order = n.array([10,11,9,12,8,20,0,13,7,14,6,15,5,16,4,17,3,18,2,19,1])
    iorder = n.argsort(order)
    FC_o = n.take(n.take(FC,order, axis=0), order, axis=1)
    L_o = n.linalg.cholesky(FC_o)
    U,S,V = n.linalg.svd(L_o.conj())
    MC_o = n.dot(n.transpose(V), n.dot(n.diag(1./S), n.transpose(U)))
    MC = n.take(n.take(MC_o,iorder, axis=0), iorder, axis=1)
    MI  = n.identity(nchan, dtype=n.complex128)
    #MC = n.linalg.inv(FC) #different choice for MC
    #MI = MC #use this if always want MC instead of MI (can choose C=I later still though)

    print "   Getting W"
    #print 'Normalizing M/W'
    WI = n.dot(MI, FI)
    norm  = WI.sum(axis=-1); norm.shape += (1,)
    #norm  = WI.max(axis=-1); norm.shape += (1,) # XXX
    MI /= norm; WI = n.dot(MI, FI)
    if PLOT:
        capo.arp.waterfall(MI,mode='real');p.colorbar(shrink=.5)
        p.show()
    WC = n.dot(MC, FC)
    norm  = WC.sum(axis=-1); norm.shape += (1,)
    #norm  = WC.max(axis=-1); norm.shape += (1,) # XXX
    MC /= norm; WC = n.dot(MC, FC)

    print '   Generating ps'
    #if opts.noise_only: scalar = 1
    pC = n.dot(MC, qC) * scalar
    #pC[m] *= 1.81 # signal loss, high-SNR XXX
    #pC[m] *= 1.25 # signal loss, low-SNR XXX
    #MI = fractional_matrix_power(FI,-0.5)
    pI = n.dot(MI, qI) * scalar 

    print 'pI ~ ', n.median(pI)    
    print 'pC ~ ', n.median(pC)
 
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

    if len(opts.output) > 0: outpath = opts.output+'/pspec_boot%04d.npz' % boot
    else: outpath = 'pspec_boot%04d.npz' % boot
    print '   Writing '+outpath
    n.savez(outpath, kpl=kpl, scalar=scalar, times=n.array(lsts),
        pk_vs_t=pC, err_vs_t=1./cnt, temp_noise_var=var, nocov_vs_t=pI,
        afreqs=afreqs,chans=chans,cmd=' '.join(sys.argv))


