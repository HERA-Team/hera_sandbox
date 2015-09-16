#! /usr/bin/env python
import matplotlib
#matplotlib.use('Agg')
from IPython import embed
from scipy import interpolate
from matplotlib.widgets import Slider
import aipy as a, numpy as n, pylab as p, capo
import glob, optparse, sys, random
import ipdb

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
o.add_option('-b', '--nboot', type='int', default=20,
    help='Number of bootstraps.  Default is 20')
o.add_option('--plot', action='store_true',
    help='Generate plots')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
o.add_option('--sep', default='sep0,1', action='store',
    help='Which separation type?')
o.add_option('--loss', action='store', 
    help='In signal loss mode to measure the signal loss. Uses default data in my path. Give it the path to the simulated signal data. Assumes ends in ')
o.add_option('--level', type='float', default=-1.0,
    help='Scalar to multiply the default signal level for simulation runs.')
o.add_option('--rmbls', action='store', 
    help='List of baselines, in miriad format, to remove from the power spectrum analysis.')
o.add_option('--output', type='string', default='',
    help='output directory for pspec_boot files (default "")')
o.add_option('--band', default='80_150', action='store',
    help='Channels from which to Calculate full band Covariance')
o.add_option('--auto',  action='store_true',
    help='Auto-scale covariance matrix')
o.add_option('--noise',  type="float", default=0,
    help='Creates White Noise spectrum. Input is rms in mk')
o.add_option('--diff',  action='store_true',
    help='Differences Even and Odd data to create Noise Estimator')

opts,args = o.parse_args(sys.argv[1:])

random.seed(0)
POL = 'I'
LST_STATS = False
DELAY = False
NGPS = 5
INJECT_SIG = 0.
SAMPLE_WITH_REPLACEMENT = True
#NOISE = .0
PLOT = opts.plot
FRF_WIDTH=401
NOISE= opts.noise
try:
    rmbls = map(int, opts.rmbls.split(','))
except:
    rmbls = []

if opts.loss:
    if opts.level >= 0.0:
        INJECT_SIG = opts.level
        print 'Running in signal loss mode, with an injection signal of %s*default level'%(opts.level)
    else:
        print 'Exiting. If in signal loss mode, need a signal level to input.'
        exit()

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

SEP = opts.sep
#dsets = {
#
#    'only': glob.glob('sep0,1/*242.[3456]*uvL'),
#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even/'+SEP+'/*242.[3456]*uvAL'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd/'+SEP+'/*243.[3456]*uvAL'),
#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even_nomni/'+SEP+'/*242.[3456]*uvALG'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd_nomni/'+SEP+'/*243.[3456]*uvALG'),
#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even/'+SEP+'/*242.[3456]*uvAFG'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd/'+SEP+'/*243.[3456]*uvAFG'),
#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even/'+SEP+'/*242.[3456]*uvALG'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd/'+SEP+'/*243.[3456]*uvALG'),

#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even_xtalk_removed/'+SEP+'/*242.[3456]*uvGL'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd_xtalk_removed/'+SEP+'/*243.[3456]*uvGL'),

#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even_nomni_xtalk/'+SEP+'/*242.[3456]*uvGL'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd_nomni_xtalk/'+SEP+'/*243.[3456]*uvGL'),
#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even_xtalk_removed/'+SEP+'/*242.[3456]*uvGF'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd_xtalk_removed/'+SEP+'/*243.[3456]*uvGF'),
#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even_xtalk_removed/'+SEP+'/*242.[3456]*uvGL'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd_xtalk_removed/'+SEP+'/*243.[3456]*uvGL'),
#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even_xtalk_removed_optimal/'+SEP+'/*242.[3456]*uvGL'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd_xtalk_removed_optimal/'+SEP+'/*243.[3456]*uvGL'),


#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even_fg/*242.[3456]*uvA'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd_fg/*243.[3456]*uvA'),
#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/signal_loss/data/even/*242.[3456]*uvALG'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/signal_loss/data/odd/*243.[3456]*uvALG'),
#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/signal_loss/signal/even/*242.[3456]*uv_perf'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/signal_loss/signal/odd/*243.[3456]*uv_perf'),
#    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_even/'+SEP+'/*242.[3456]*uvALG_signalL'),
#    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/lstbin_odd/'+SEP+'/*243.[3456]*uvALG_signalL'),
#}
#for i in xrange(10): dsets[i] = glob.glob('lstbinX%d/%s/lst.24562[45]*.[3456]*.uvAL'%(i,SEP))
#dsets = {
#    #'only': glob.glob('sep0,1/*242.[3456]*uvL'),
#    'even': glob.glob('/home/mkolopanis/psa64/lstbin_even_noxtalk/'+SEP+'/*242.[3456]*uvGL'),
#    'odd' : glob.glob('/home/mkolopanis/psa64/lstbin_odd_noxtalk/'+SEP+'/*243.[3456]*uvGL'),
#}
#dsets['even'].sort()
#dsets['odd'].sort()
dsets = {
    'even': [x for x in args if 'even' in x],
    'odd' : [x for x in args if 'odd' in x]
}

print 'Number of even data sets: {0:d}'.format(len(dsets['even']))
print 'Number of odd data sets: {0:d}'.format(len(dsets['odd']))
for dset_count in xrange(len(dsets['even'])):
        print dsets['even'][dset_count].split('/')[-1], dsets['odd'][dset_count].split('/')[-1]
#sys.exit()
if opts.loss:
    dsets = {
    'even': glob.glob('/home/mkolopanis/psa64/lstbin_even_noxtalk/sep0,1/*242.[3456]*uvGL'),
    'odd' : glob.glob('/home/mkolopanis/psa64/lstbin_odd_noxtalk/sep0,1/*243.[3456]*uvGL'),
}

WINDOW = opts.window
uv = a.miriad.UV(dsets.values()[0][0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
band_chans = a.scripting.parse_chans(opts.band, uv['nchan'])
inttime = uv['inttime'] * 4 # XXX hack for *E files that have inttime set incorrectly
#inttime=8*60.
del(uv)

afreqs = freqs.take(chans)
allfreqs = freqs.take(band_chans)
nchan_band = len(band_chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

aa = a.cal.get_aa(opts.cal, allfreqs)
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
jy2T = capo.pspec.jy2T(allfreqs)
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
antstr = 'cross'
lsts,data,flgs = {},{},{}
days = dsets.keys()
for k in days:
    lsts[k],data[k],flgs[k] = get_data(dsets[k], antstr=antstr, polstr=POL, rmbls=rmbls, verbose=True)
    print data[k].keys()

if LST_STATS:
    # collect some metadata from the lst binning process
    cnt, var = {}, {}
    for filename in dsets.values()[0]:
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

if True:
#if False:
    # Align data sets in LST
    print [lsts[k][0] for k in days]
    lstmax = max([lsts[k][0] for k in days])
    for k in days:
        print k
        for i in xrange(len(lsts[k])):
            # allow for small numerical differences (which shouldn't exist!)
            if lsts[k][i] >= lstmax - .001: break
        lsts[k] = lsts[k][i:]
        for bl in data[k]:
            data[k][bl],flgs[k][bl] = data[k][bl][i:],flgs[k][bl][i:]
    print [len(lsts[k]) for k in days]
    j = min([len(lsts[k]) for k in days])
    for k in days:
        lsts[k] = lsts[k][:j]
        for bl in data[k]:
            data[k][bl],flgs[k][bl] = n.array(data[k][bl][:j]),n.array(flgs[k][bl][:j])
else:
    for k in days:
        for bl in data[k]:
            data[k][bl], flgs[k][bl] = n.array(data[k][bl][:]), n.array(flgs[k][bl][:])
lsts = lsts.values()[0]
index = map(lambda x: n.where( x == band_chans)[0][0], chans)
x = {}
x1 = {}
print len(data[k][bl])
print type(chans)
for k in days:
    x[k] = {}
    x1[k] = {}
    for bl in data[k]:
        print k, bl
        d = data[k][bl][:,band_chans] *jy2T
        d1 = n.copy(d)
        d[:,index] = 0. + 0j
        if conj[bl]:
             d = n.conj(d)
             d1 = n.conj(d1)
        x[k][bl] = n.transpose(d, [1,0]) # swap time and freq axes
        x1[k][bl] = n.transpose(d1, [1,0]) # swap time and freq axes

bls_master = x.values()[0].keys()
nbls = len(bls_master)
print 'Baselines:', nbls

wnx={}
if NOISE > 0:
        print 'Creating white noise at {0} mk rms'.format(NOISE)
        tmp_nx = {}
        bl1 = a.miriad.bl2ij(bls_master[0])
        beam_w_fr = capo.frf_conv.get_beam_w_fr(aa,bl1,ref_chan=0)
        t, firs, frbins, frspace = capo.frf_conv.get_fringe_rate_kernels(beam_w_fr, inttime, FRF_WIDTH)
        for k in days:
            wnx[k]={}
            tmp_nx[k] = {}
            for bl in x[k]:
                noise1 = noise(x[days[0]][bls_master[0]].shape)*6.6408 * NOISE ### mk 
                for cnt, ch in enumerate(band_chans):
                    noise1[cnt] = n.convolve(noise1[cnt], firs[cnt], mode='same')

                tmp_nx[k][bl] = noise1.copy() 
                wnx[k][bl]=noise1.copy()

        #if set(['even', 'odd']) == set(days):
        #    k1,k2 = x.keys()
        #    wnx[k1], wnx[k2] = {},{}
        #    for bl in bls_master:



if INJECT_SIG > 0.: # Create a fake EoR signal to inject
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
        #cut with same lst cut.
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
#            p.figure(1)
#            p.plot(x[k][bl])
#            p.figure(2)
#            p.plot(eor[k][bl])
#            p.show()
            x[k][bl] += eor[k][bl] 
    
    if PLOT:
        capo.arp.waterfall(x[k][bl], mode='real'); p.colorbar(); p.show()
       
    
    
    
    
#    eor = noise(x.values()[bls_master[0]].shape) * INJECT_SIG
#    fringe_filter = n.ones((44,))
#    # Maintain amplitude of original noise
#    fringe_filter /= n.sqrt(n.sum(fringe_filter))
#    for ch in xrange(eor.shape[0]):
#        eor[ch] = n.convolve(eor[ch], fringe_filter, mode='same')
#    _eor = n.fft.ifft(eor, axis=0)
#    #wgt = n.exp(-n.fft.ifftshift(kpl)**2/(2*.3**2))
#    wgt = n.zeros(_eor.shape[0]); wgt[0] = 1
#    wgt.shape = wgt.shape + (1,)
#    #_eor *= wgt
#    #_eor = n.fft.ifft(eor, axis=0); _eor[4:-3] = 0
#    #eor = n.fft.fft(_eor, axis=0)
#    eor *= wgt

#Q = {} # Create the Q's that extract power spectrum modes
#for i in xrange(nchan):
#    Q[i] = get_Q(i, nchan)
Q = [get_Q(i,nchan) for i in xrange(nchan)]

if False: ##Looking for signal below noise level
    print 'Averaging over Baselines'
    x_blavg, _Cx_blavg = {},{}
    C_blavg ,_C_blavg = {},{}
    for k in days:
        x_blavg[k]=n.mean( [x1[k][bl] for bl in bls_master],axis=0)
        C_blavg[k] = cov(x_blavg[k])
        #U,S,V= n.linalg.svd(C_blavg[k].conj())
        #_C_blavg[k] =  n.einsum('ij,j,jk', V.T, 1./S,U.T)
        #norm  = _C_blavg[k].sum(axis=-1); norm.shape += (1,)
        #_C_blavg[k] /= norm
        #_Cx_blavg[k] = n.dot(_C_blavg[k],x_blavg[k])[index]
    if PLOT and True:
        for k in days:
            p.subplot(211); capo.arp.waterfall(x_blavg[k], mode='real'); p.colorbar()
            p.ylabel('X')
            p.subplot(212); capo.arp.waterfall(C_blavg[k], mode='real'); p.colorbar()
            p.ylabel('Cov')
            #p.xticks([])
            #p.subplot(324); capo.arp.waterfall(_C_blavg[k], mode='real',mx=1.2,drng=2.4);# p.colorbar()
            #p.xticks([])
            #p.subplot(313); capo.arp.waterfall(_Cx_blavg[k], mode='real')
            p.subplots_adjust(right=.9,left=.1)
            p.suptitle('Baseline Averaged '+k)
            p.show()


##MAKE Average covariance
I,_I,_Ix = {},{},{}
C,Cav,Cst,_C,_Cav,_Cavx,_Cx = {},{},{},{},{},{},{}
x_blavg, x1_blavg= {},{}
W = {}
for k in days:
    I[k],_I[k],_Ix[k] = {},{},{}
    C[k],Cav[k],Cst[k],_Cav[k],_C[k],_Cavx[k],_Cx[k] = {},{},{},{},{},{},{}
    x_blavg[k]=n.mean( [x[k][bl] for bl in bls_master],axis=0)
    x1_blavg[k]=n.mean( [x1[k][bl] for bl in bls_master],axis=0)
    W[k]= {}
    for bl in x[k]:
        C[k][bl] = cov(x_blavg[k])
        I[k][bl] = n.identity(C[k][bl].shape[0])
        U,S,V = n.linalg.svd(C[k][bl].conj())
        _C[k][bl] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
        _I[k][bl] = n.identity(_C[k][bl].shape[0])
        _Cx[k][bl] = n.dot(_C[k][bl], x1[k][bl])
#/n.dot(_C[k][bl],n.ones_like(x1[k][bl]))
        _Ix[k][bl] = x[k][bl].copy()

        #Compute Cav by taking diags, averaging and reforming the matrix
        diags =[C[k][bl].diagonal(count) for count in xrange(nchan_band-1, -nchan_band,-1)]
        Cav[k][bl]=n.zeros_like(C[k][bl])
        Cst[k][bl]=n.zeros_like(C[k][bl])

        for count,count_chan in enumerate(range(nchan_band-1,-nchan_band,-1)):
                Cav[k][bl] += n.diagflat( n.mean(diags[count]).repeat(len(diags[count])), count_chan)
                Cst[k][bl] += n.diagflat( n.sqrt( n.mean( (diags[count]-n.mean(diags[count]) ) *n.conj(diags[count] - n.mean(diags[count]) )  )).repeat(len(diags[count])), count_chan)

        new_c = cov(x1_blavg[k])
        U2,S2,V2=n.linalg.svd(new_c.conj())
        _new_c=n.einsum('ij,j,jk',V2.T,1./S2,U2.T)
        #max_diff=  abs( Cav[k][bl] - new_c)/new_c 
        #print n.shape(max_diff),n.max(max_diff)
        #if n.max(max_diff) >= .15:
        #    print 'Setting Diagonal to Auto-Covariance'
        #    opts.auto = True
        #else:
        #    print 'Not Normalizing by Auto-Covariance'
        #    opts.auto = False

        if opts.auto:
            #Need full covariance for auto-covariance
            #points=C[k][bl].nonzero()

            #scale rows and columns by sqrt of auto-covariance
            for count in xrange(nchan_band):
                tmp=n.copy(Cav[k][bl])
                Cav[k][bl][count,:] *= n.sqrt(new_c[count,count]/tmp[count,count])
                Cav[k][bl][:,count] *= n.sqrt(new_c[count,count]/tmp[count,count])
            #Ensure symmetric matrix, some interpolation alogrithms 
            #have artefacts
            #Cav[k][bl] /= n.mean(Cav[k][bl].diagonal())
            #Cav[k][bl] = n.array(Cav[k][bl]+Cav[k][bl].T)/2.
        U1,S1,V1 = n.linalg.svd(Cav[k][bl].conj())
        S1[35:]=n.Inf
        _Cav[k][bl] = n.einsum('ij,j,jk', V1.T, 1./S1, U1.T)

        #norm  = n.sqrt(n.sum(_Cav[k][bl]**2,axis=-1)); norm.shape += (1,)
        #_Cav[k][bl] /= norm
        W[k][bl] = _Cav[k][bl].copy()
        W[k][bl] /= n.sqrt(n.dot(W[k][bl].T,W[k][bl]).sum(axis=-1))
        
        eig_freqs=n.fft.fftshift(n.fft.fft(n.einsum('ij,jk', n.diag(S1),V1).T.real,axis=0))


        if not n.allclose( Cav[k][bl], Cav[k][bl].T.conj()):
            good_values=n.isclose(Cav[k][bl].Cav[k][bl].T.conj())
            bad_values=~good_values
            print('There are {0:d} elements which do not match'.format(n.sum(bad_values)))
            sys.exit(0)

        _Cavx[k][bl] = n.dot(W[k][bl], x1[k][bl]) 
        #print('Cov:',n.max(_new_c),n.min(_new_c))
        #print('Cav:', n.max(_Cav[k][bl]),n.min(_Cav[k][bl]))
#/ n.dot(_Cav[k][bl], n.ones(n.shape(x1[k][bl])))
        #_Cavx[k][bl] = n.zeros_like(_Cx[k][bl])
        #_Cavx[k][bl][index] += temp[index]

        x_power= n.sum(x1[k][bl].conj() * x1[k][bl],axis=0)
        c_power= n.sum(_Cavx[k][bl].conj() * _Cavx[k][bl],axis=0)


        if PLOT and True:
            #p.plot(S); p.show()
                f1=p.figure(1)
                p.subplot(511); capo.arp.waterfall(x1[k][bl], mode='real',mx=6,drng=12);
                p.subplot(534); capo.arp.waterfall(C[k][bl],mx=1.2,drng=2.4); p.ylabel('C')
                p.subplot(535); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
                p.subplot(536); capo.arp.waterfall(_C[k][bl])
                p.subplot(537); capo.arp.waterfall(Cav[k][bl],mx=1.2,drng=2.4); p.ylabel('Cav')
                #p.subplot(537); capo.arp.waterfall(C[k][bl] - Cav[k][bl],mx=1.2,drng=2.4); p.ylabel('Cav')
                p.subplot(538); p.plot(n.einsum('ij,jk',n.diag(S1),V1).T.real)
                p.subplot(539); capo.arp.waterfall(_Cav[k][bl])
                p.subplot(514); capo.arp.waterfall(_Cx[k][bl], mode='real',mx=6,drng=12); p.ylabel('_Cx')
                p.subplot(515); capo.arp.waterfall(_Cavx[k][bl], mode='real',mx=6,drng=12); p.ylabel('_Cavx')
                p.suptitle('%d_%d'%a.miriad.bl2ij(bl))
                #f1.savefig('Cov_Baseline_%d_%d'%a.miriad.bl2ij(bl))
                #f1.clf()
                f2=p.figure(2);
                p.plot(S2,label='C')
                p.plot(S1,label='Cav')
                p.yscale('log')
                p.legend(loc='best')
                #f2.savefig('Cov_Eigenvalues_%d_%d'%a.miriad.bl2ij(bl))
                #f2.clf()
                f3=p.figure(3)
                p.subplot(231); capo.arp.waterfall(new_c,mode='real', mx=2,drng=4)
                p.title('C_xblavg')
                p.subplot(232); capo.arp.waterfall(Cav[k][bl],mode='real',mx=2,drng=4)
                p.title('Cav_xblavg')
                p.subplot(233); capo.arp.waterfall(new_c-Cav[k][bl],mode='real',mx=2,drng=4); p.colorbar()
                p.title('resid')

                p.subplot(234); capo.arp.waterfall(_new_c,mode='real')#, mx=12,drng=24)
                p.subplot(235); capo.arp.waterfall(_Cav[k][bl])#,mx=12,drng=24)
                p.subplot(236); capo.arp.waterfall(_new_c-_Cav[k][bl])#,mode='real',mx=12,drng=24); p.colorbar()
                #f3.savefig('Cov_Residual_%d_%d'%a.miriad.bl2ij(bl))
                #f3.clf()
                f4=p.figure(4)
                p.plot(eig_freqs)
                f5=p.figure(5)
                p.subplot(211)
                p.xticks([])
                p.plot(x_power)
                p.plot(c_power)
                p.subplot(212)
                p.subplots_adjust(hspace=0)
                p.plot(c_power/x_power)
                f6=p.figure(6)
                p.subplot(131)
                capo.arp.waterfall(W[k][bl],mode='real',mx=2,drng=4)
                p.subplot(132)
                capo.arp.waterfall(n.dot(W[k][bl].T,W[k][bl]),mode='real',mx=2,drng=4)
                p.subplot(133)
                capo.arp.waterfall(n.dot(W[k][bl].T,W[k][bl]) - n.ones_like(W[k][bl]),mode='real',mx=2,drng=4)
                p.colorbar()
                if False: 
                    fig=p.figure(4); 
                    S2_=n.zeros_like(S2)
                    S2_[-1:]+=S2[-1:]
                    C1 = n.einsum('ij,j,jk', U2, S2_, V2)
                    S1_=n.zeros_like(S1)
                    S1_[-1:]+=S1[-1:]
                    Cav1= n.einsum('ij,j,jk', U1, S1_, V1)

                    p.subplot(231); capo.arp.waterfall(new_c, mode='real',mx=7,drng=14); p.ylabel('C')
                    p.subplot(232);l= p.plot(n.einsum('ij,jk',n.diag(S2_),V2).T.real) 
                    p.subplot(233); capo.arp.waterfall(C1,mode='real',mx=7,drng=14) 
                    p.subplot(234); capo.arp.waterfall(Cav[k][bl], mode='real',mx=7,drng=14);  p.ylabel('Cav')
                    
                    p.subplot(235);l1=p.plot(n.einsum('ij,jk',n.diag(S1_),V1).T.real) 
                    p.subplot(236); capo.arp.waterfall(Cav1,mode='real',mx=7,drng=14) 
                    axeig=p.axes([0.15,.01,.65,.03])
                    seig = Slider(axeig, 'Eig val',1,71,valinit=1, valfmt='%0.0f')
                    #embed()
                    def update(val):
                        S2_=n.zeros_like(S2)
                        S1_=n.zeros_like(S1)
                        s= int(seig.val)
                        S2_[-s:]+=S2[-s:]
                        S1_[-s:]+=S1[-s:]
                        C1 = n.einsum('ij,j,jk', U2, S2_, V2)
                        Cav1= n.einsum('ij,j,jk', U1, S1_, V1)
                        vals1 = n.einsum('ij,jk',n.diag(S2_),V2).T.real
                        valsav1=n.einsum('ij,jk', n.diag(S1_),V1).T.real
                        for j in xrange(len(l)):
                            l[j].set_ydata(vals1[:,j])
                            l1[j].set_ydata(valsav1[:,j])
                        fig.axes[1].relim()
                        fig.axes[1].autoscale_view()
                        fig.axes[4].relim()
                        fig.axes[4].autoscale_view()
                        p.subplot(233); capo.arp.waterfall(C1,mode='real',mx=7,drng=14) 
                        p.subplot(236); capo.arp.waterfall(Cav1, mode='real',mx=7,drng=14); 
                        fig.canvas.draw_idle()
                    seig.on_changed(update)
                p.show()
                p.close



print 'Differencing Even and Odd Days'
nx={}
if set(['even', 'odd']) == set(days):
    k1,k2 = x.keys()
    nx[k1], nx[k2] = {},{}
    for bl in bls_master:
        #Divide by Sqrt(2) because differnce gives 2*noise^2 in power
        nx[k1][bl] = n.copy(_Cavx[k1][bl] - _Cavx[k2][bl])/n.sqrt(2)
        nx[k2][bl] = n.copy(_Cavx[k1][bl] - _Cavx[k2][bl])/n.sqrt(2)

if opts.diff:
    if set(['even','odd']) == set(days):
        k1,k2=x.keys()
        wnx[k1],wnx[k2] = {},{}
        for bl in bls_master:
            wnx[k1][bl] = n.copy(tmp_nx[k1][bl] - tmp_nx[k2][bl])/n.sqrt(2)
            wnx[k2][bl] = n.copy(tmp_nx[k1][bl] - tmp_nx[k2][bl])/n.sqrt(2)
# Compute baseline auto-covariances and apply inverse to data
x = {}
grid = n.meshgrid(index,index)
I,_I,_Ix = {},{},{}
_Inx, _Iwnx, _Icavx = {},{},{}
C,_C,_Cx = {},{},{}
if True:
    for k in days:
        x[k] = {}
        for bl in data[k]:
            print k, bl
            x[k][bl] = x1[k][bl][index].copy()
            nx[k][bl] = nx[k][bl][index].copy()   
            wnx[k][bl] = wnx[k][bl][index].copy()   
    for k in days:
        I[k],_I[k],_Ix[k] = {},{},{}
        _Inx[k], _Iwnx[k], _Icavx[k] = {},{},{}
        C[k],_C[k],_Cx[k] = {},{},{}
        for bl in x[k]:
            C[k][bl] = cov(x[k][bl])
            I[k][bl] = n.identity(C[k][bl].shape[0])
            U,S,V = n.linalg.svd(C[k][bl].conj())
            _C[k][bl] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
            _I[k][bl] = n.identity(_C[k][bl].shape[0])
            _Cx[k][bl] = n.dot(_C[k][bl],x[k][bl])
            _Ix[k][bl] = x[k][bl].copy()
            _Icavx[k][bl] = _Cavx[k][bl][index].copy()
            _Inx[k][bl] = nx[k][bl].copy()
            _Iwnx[k][bl] = wnx[k][bl].copy()
else:
    for k in days:
        x[k] = {}
        for bl in data[k]:
            print k, bl
            x[k][bl] = _Cavx[k][bl][index].copy()
            nx[k][bl] = nx[k][bl][index].copy()
            wnx[k][bl] = wnx[k][bl][index].copy()   
    for k in days:
        I[k],_I[k],_Ix[k] = {},{},{}
        _Inx[k], _Iwnx[k], _Icavx[k] = {},{},{}
        C[k],_C[k],_Cx[k] = {},{},{}
        for bl in x[k]:
            C[k][bl] = cov(x[k][bl])
            I[k][bl] = n.identity(C[k][bl].shape[0])
            U,S,V = n.linalg.svd(C[k][bl].conj())
            _C[k][bl] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
            _I[k][bl] = n.identity(_C[k][bl].shape[0])
            _Cx[k][bl] = n.dot(_C[k][bl], x[k][bl])
            _Ix[k][bl] = x[k][bl].copy()
            _Icavx[k][bl] = x[k][bl].copy()
            _Inx[k][bl] = nx[k][bl].copy()
            _Iwnx[k][bl] = wnx[k][bl].copy()

Nt=13
nlst=n.shape(lsts)[0]
pk_array=[]
rms_array=[]

bls_master.sort()
#for bl in bls_master:
#    dlst= n.ceil(nx['even'][bl][index].shape[1]/13)
#    dlst= n.ceil(nlst/13.)
#    pk_array.append( n.mean(nx['even'][bl][index,::dlst]*nx['even'][bl][index,::dlst].conj())*scalar )
#    rms_array.append( n.sqrt(n.mean(nx['even'][bl][index[0]]*nx['even'][bl][index[0]].conj()))  ) 
#    print '', bl, ' Pk [mk^2/Mpc^3] ', n.sqrt(n.mean(nx['even'][bl][index[0]].conj()*nx['even'][bl][index[0]])) * scalar, 'rms [mk]  ', n.sqrt(n.mean(nx['even'][bl][index[0]].conj()*nx['even'][bl][index[0]]))
#
#print('dlst: {0}'.format(dlst))
#diff_blavg= n.mean( [nx['even'][bl] for bl in bls_master],axis=0)
#Trms_blavg= n.sqrt( n.mean(diff_blavg[index,::dlst]*diff_blavg[index,::dlst].conj()))
#print( 'BL Average Pk [mk^2/Mpc^3]: {0:3e}'.format(Trms_blavg.real**2*scalar/n.sqrt(nbls*13.)))
#print('Average rms [mk]: {0:f}'.format(n.mean(rms_array).real))
#print('BL Average rms [mk]: {0:f}'.format(Trms_blavg.real))
#
#
#print('\n\n')
#print('White Noise RMS')
#for bl in bls_master:
#    dlst= n.ceil(wnx['even'][bl][index].shape[1]/13)
#    dlst= n.ceil(nlst/13.)
#    pk_array.append( n.mean(wnx['even'][bl][index,::dlst]*wnx['even'][bl][index,::dlst].conj())*scalar )
#    rms_array.append( n.sqrt(n.mean(wnx['even'][bl][index[0]]*wnx['even'][bl][index[0]].conj()))  ) 
#    print '', bl, ' Pk [mk^2/Mpc^3] ', n.sqrt(n.mean(wnx['even'][bl][index[0]].conj()*wnx['even'][bl][index[0]])) * scalar, 'rms [mk]  ', n.sqrt(n.mean(wnx['even'][bl][index[0]].conj()*wnx['even'][bl][index[0]]))
#
#print('dlst: {0}'.format(dlst))
#diff_blavg= n.mean( [wnx['even'][bl] for bl in bls_master],axis=0)
#Trms_blavg= n.sqrt( n.mean(diff_blavg[index,::dlst]*diff_blavg[index,::dlst].conj()))
#print( 'BL Average Pk [mk^2/Mpc^3]: {0:3e}'.format(Trms_blavg.real**2*scalar/n.sqrt(nbls*13.)))
#print('Average rms [mk]: {0:f}'.format(n.mean(rms_array).real))
#print('BL Average rms [mk]: {0:f}'.format(Trms_blavg.real))
#
#if PLOT and False:
#    capo.arp.waterfall(nx['even'][1555][index], mode='real')
#    p.show()


for boot in xrange(opts.nboot):
    print '%d / %d' % (boot+1,opts.nboot)
    bls = bls_master[:]
    if True: # shuffle and group baselines for bootstrapping
        if not SAMPLE_WITH_REPLACEMENT:
            random.shuffle(bls)
            bls = bls[:-1] # XXX
        else: # sample with replacement
            bls = [random.choice(bls) for bl in bls]
        gps = [bls[i::NGPS] for i in range(NGPS)]
        gps = [[random.choice(gp) for bl in gp] for gp in gps]
    else: # assign each baseline its own group
        #gps = [[bl] for bl in bls]
        #gps = [bls]
        gps = [bls[i::NGPS] for i in range(NGPS)]
    #gps = [[bl for bl in gp] for gp in gps]
    bls = [bl for gp in gps for bl in gp]
    print '\n'.join([','.join(['%d_%d'%a.miriad.bl2ij(bl) for bl in gp]) for gp in gps])
    _Iz,_Isum,_IsumQ = {},{},{}
    _Inz, _Iwnz, _Icavz = {},{},{}
    _Cz,_Csum,_CsumQ = {},{},{}
    Csum = {}
    for k in days:
        _Iz[k],_Isum[k],_IsumQ[k] = {},{},{}
        _Inz[k], _Iwnz[k], _Icavz[k] = {},{},{}
        _Cz[k],_Csum[k],_CsumQ[k] = {},{},{}
        Csum[k]={}
        for i,gp in enumerate(gps):
            _Iz[k][i] = sum([_Ix[k][bl] for bl in gp])
            _Icavz[k][i] = sum([_Icavx[k][bl] for bl in gp])
            _Inz[k][i] = sum([_Inx[k][bl] for bl in gp])
            _Iwnz[k][i] = sum([_Iwnx[k][bl] for bl in gp])
            _Cz[k][i] = sum([_Cx[k][bl] for bl in gp])
            _Isum[k][i] = sum([_I[k][bl] for bl in gp])
            _Csum[k][i] = sum([_C[k][bl] for bl in gp])
            Csum[k][i] = sum([C[k][bl] for bl in gp])
            _IsumQ[k][i] = {}
            _CsumQ[k][i] = {}
            if DELAY: # this is much faster
                _Iz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Iz[k][i], axis=0), axes=0)
                _Icavz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Iz[k][i], axis=0), axes=0)
                _Inz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Inz[k][i], axis=0), axes=0)
                _Iwnz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Iwnz[k][i], axis=0), axes=0)
                _Cz[k][i] = n.fft.fftshift(n.fft.ifft(window*_Cz[k][i], axis=0), axes=0)
                # XXX need to take fft of _Csum, _Isum here
            for ch in xrange(nchan): # XXX this loop makes computation go as nchan^3
                _IsumQ[k][i][ch] = n.dot(_Isum[k][i], Q[ch])
                _CsumQ[k][i][ch] = n.dot(_Csum[k][i], Q[ch])
        if PLOT and  False:
            NGPS = len(gps)
            _Csumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            Csumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            _Isumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            for i in xrange(len(gps)): _Isumk[i,:,i,:] = _Isum[k][i]
            _Isumk.shape = (NGPS*nchan, NGPS*nchan)
            #_Isum[k] = _Isumk
            for i in xrange(len(gps)):
                     _Csumk[i,:,i,:] = _Csum[k][i]
                     Csumk[i,:,i,:] = Csum[k][i] 
            _Csumk.shape = (NGPS*nchan, NGPS*nchan)
            Csumk.shape = (NGPS*nchan, NGPS*nchan)
            #_Csum[k] = _Csumk
            _Czk = n.array([_Cz[k][i] for i in _Cz[k]])
            _Izk = n.array([_Iz[k][i] for i in _Iz[k]])
            _Icavzk = n.array([_Icavz[k][i] for i in _Icavz[k]])
            _Inzk = n.array([_Inz[k][i] for i in _Inz[k]])
            _Iwnzk = n.array([_Iwnz[k][i] for i in _Iwnz[k]])
            _Czk = n.reshape(_Czk, (_Czk.shape[0]*_Czk.shape[1], _Czk.shape[2]))
            _Izk = n.reshape(_Izk, (_Izk.shape[0]*_Izk.shape[1], _Izk.shape[2]))
            _Icavzk = n.reshape(_Icavzk, (_Icavzk.shape[0]*_Icavzk.shape[1], _Icavzk.shape[2]))
            _Inzk = n.reshape(_Inzk, (_Inzk.shape[0]*_Inzk.shape[1], _Inzk.shape[2]))
            _Iwnzk = n.reshape(_Iwnzk, (_Iwnzk.shape[0]*_Iwnzk.shape[1], _Iwnzk.shape[2]))
            C_I=cov(_Izk)
            I_U,I_S,I_V=n.linalg.svd(C_I)
            _C_I=n.einsum('ij,j,jk',I_V.T,1./I_S,I_U.T)
            p.subplot(411); capo.arp.waterfall(_Izk, mode='real')
            p.subplot(423); capo.arp.waterfall(Csumk)
            p.subplot(424); capo.arp.waterfall(_Csumk)
            p.subplot(425); capo.arp.waterfall(C_I)
            p.subplot(426); capo.arp.waterfall(_C_I)
            #p.subplot(426); capo.arp.waterfall(cov(_Czk))
            p.subplot(414); capo.arp.waterfall(_Czk, mode='real')
            fig_file='Data_Covariance_boot{0:0>4d}'.format(boot) 
            if not opts.output == '':
                fig_file= opts.output + '/' + fig_file
            p.savefig(fig_file)
            p.close()

    FI = n.zeros((nchan,nchan), dtype=n.complex)
    FC = n.zeros((nchan,nchan), dtype=n.complex)
    qI = n.zeros((nchan,_Iz.values()[0].values()[0].shape[1]), dtype=n.complex)
    qCav = n.zeros((nchan,_Icavz.values()[0].values()[0].shape[1]), dtype=n.complex)
    qN = n.zeros((nchan,_Inz.values()[0].values()[0].shape[1]), dtype=n.complex)
    qW = n.zeros((nchan,_Iwnz.values()[0].values()[0].shape[1]), dtype=n.complex)
    qC = n.zeros((nchan,_Cz.values()[0].values()[0].shape[1]), dtype=n.complex)
    Q_Iz = {}
    Q_Icavz = {}
    Q_Inz = {}
    Q_Iwnz = {}
    Q_Cz = {}
    for cnt1,k1 in enumerate(days):
        for k2 in days[cnt1:]:
            if not Q_Iz.has_key(k2): Q_Iz[k2] = {}
            if not Q_Icavz.has_key(k2): Q_Icavz[k2] = {}
            if not Q_Inz.has_key(k2): Q_Inz[k2] = {}
            if not Q_Iwnz.has_key(k2): Q_Iwnz[k2] = {}
            if not Q_Cz.has_key(k2): Q_Cz[k2] = {}
            for bl1 in _Cz[k1]:
                for bl2 in _Cz[k2]:
                    #if k1 == k2 and bl1 == bl2: continue # this results in a significant bias
                    if k1 == k2 or bl1 == bl2: continue
                    #if k1 == k2: continue
                    #if bl1 == bl2: continue # also a significant noise bias
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
                    if False: # use ffts to do q estimation fast
                        qI += n.conj(_Iz[k1][bl1]) * _Iz[k2][bl2]
                        qI += n.conj(_Icavz[k1][bl1]) * _Icavz[k2][bl2]
                        qN += n.conj(_Inz[k1][bl1]) * _Inz[k2][bl2]
                        qW += n.conj(_Iwnz[k1][bl1]) * _Iwnz[k2][bl2]
                        qC += n.conj(_Cz[k1][bl1]) * _Cz[k2][bl2]
                    else: # brute force with Q to ensure normalization
                        #_qI = n.array([_Iz[k1][bl1].conj() * n.dot(Q[i], _Iz[k2][bl2]) for i in xrange(nchan)])
                        #_qC = n.array([_Cz[k1][bl1].conj() * n.dot(Q[i], _Cz[k2][bl2]) for i in xrange(nchan)])
                        if not Q_Iz[k2].has_key(bl2): Q_Iz[k2][bl2] = [n.dot(Q[i], _Iz[k2][bl2]) for i in xrange(nchan)]
                        if not Q_Icavz[k2].has_key(bl2): Q_Icavz[k2][bl2] = [n.dot(Q[i], _Icavz[k2][bl2]) for i in xrange(nchan)]
                        if not Q_Inz[k2].has_key(bl2): Q_Inz[k2][bl2] = [n.dot(Q[i], _Inz[k2][bl2]) for i in xrange(nchan)]
                        if not Q_Iwnz[k2].has_key(bl2): Q_Iwnz[k2][bl2] = [n.dot(Q[i], _Iwnz[k2][bl2]) for i in xrange(nchan)]
                        if not Q_Cz[k2].has_key(bl2): Q_Cz[k2][bl2] = [n.dot(Q[i], _Cz[k2][bl2]) for i in xrange(nchan)]
                        _qI = n.array([_Iz[k1][bl1].conj() * Q_Iz[k2][bl2][i] for i in xrange(nchan)])
                        qI += n.sum(_qI, axis=1)
                        
                        _qCav = n.array([_Icavz[k1][bl1].conj() * Q_Icavz[k2][bl2][i] for i in xrange(nchan)])
                        qCav += n.sum(_qCav, axis=1)
                        
                        _qN = n.array([_Inz[k1][bl1].conj() * Q_Inz[k2][bl2][i] for i in xrange(nchan)])
                        qN += n.sum(_qN, axis=1)
                        
                        _qW = n.array([_Iwnz[k1][bl1].conj() * Q_Iwnz[k2][bl2][i] for i in xrange(nchan)])
                        qW += n.sum(_qW, axis=1)
                        _qC = n.array([_Cz[k1][bl1].conj() * Q_Cz[k2][bl2][i] for i in xrange(nchan)])
                        qC += n.sum(_qC, axis=1)
                    if DELAY: # by taking FFT of CsumQ above, each channel is already i,j separated
                        FI += n.conj(_IsumQ[k1][bl1]) * _IsumQ[k2][bl2]
                        FC += n.conj(_CsumQ[k1][bl1]) * _CsumQ[k2][bl2]
                    else:
                        for i in xrange(nchan):
                            for j in xrange(nchan):
                                FI[i,j] += n.einsum('ij,ji', _IsumQ[k1][bl1][i], _IsumQ[k2][bl2][j])
                                FC[i,j] += n.einsum('ij,ji', _CsumQ[k1][bl1][i], _CsumQ[k2][bl2][j])
    if n.allclose(qN, qI):
        print('Noise equals data \n exiting')
        sys.exit(0)
    
    if n.allclose(qN, qW):
        print('Noise equals White Noise \n exiting')
        sys.exit(0)

    if PLOT and False:
        p.subplot(141); capo.arp.waterfall(FC, drng=4)
        p.subplot(142); capo.arp.waterfall(FI, drng=4)
        p.subplot(143); capo.arp.waterfall(qC, mode='real')
        p.subplot(144); capo.arp.waterfall(qI, mode='real')
        fig_file = 'FC_FI_qC_qI_boot{0:0>4d}'.format(boot)
        p.show()
        if not opts.output == '':
            fig_file= opts.output + '/' + fig_file
        p.savefig(fig_file)
        p.close()

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
    MC_o = n.dot(n.transpose(V), n.dot(n.diag(1./S), n.transpose(U)))
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
    pCav = n.dot(MI, qCav) * scalar
    pN = n.dot(MI, qN) * scalar
    pW = n.dot(MI, qW) * scalar

    if PLOT and False:
        p.subplot(411); capo.arp.waterfall(qC, mode='real'); p.colorbar(shrink=.5)
        p.subplot(412); capo.arp.waterfall(pC, mode='real'); p.colorbar(shrink=.5)
        p.subplot(413); capo.arp.waterfall(qI, mode='real'); p.colorbar(shrink=.5)
        p.subplot(414); capo.arp.waterfall(pI, mode='real'); p.colorbar(shrink=.5)
        p.show()

    if PLOT and False:
        p.plot(kpl, n.average(pC.real, axis=1), 'b.-')
        p.plot(kpl, n.average(pI.real, axis=1), 'k.-')
        p.show()

    outfile = 'pspec_boot%04d.npz'%(boot)
    if not opts.output == '':
        outfile =opts.output+'/'+outfile
    print "Writing", outfile
    n.savez(outfile, kpl=kpl, scalar=scalar, times=n.array(lsts),
        pk_vs_t=pC, err_vs_t=1./cnt, temp_noise_var=var, nocov_vs_t=pI, 
        afreqs=afreqs, chans=chans, nk_vs_t=pN, wnk_vs_t=pW, cav_vs_t=pCav,
        cmd=' '.join(sys.argv))


###    #n.savez(outfile, kpl=kpl, scalar=scalar, times=n.array(lsts),
###        pk_vs_t=pI, err_vs_t=1./cnt, temp_noise_var=var, nocov_vs_t=pI, 
###        afreqs=afreqs, chans=chans,
###        cmd=' '.join(sys.argv))
###

