#! /usr/bin/env python/
import matplotlib
#matplotlib.use('Agg')
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
o.add_option('--noise', type='float', default=0.0,
    help='amplitude of noise to inject into data')
o.add_option('--filter_noise',action='store_true',
    help='Apply fringe rate filter to injected noise')
o.add_option('--noise_only',action='store_true',
    help='Instead of injecting noise, Replace data with noise')
o.add_option('--frpad', type='float',default=1.0,
    help='frpad to use while filtering noise')
o.add_option('--rmbls', action='store', 
    help='List of baselines, in miriad format, to remove from the power spectrum analysis.')
o.add_option('--output', type='string', default='',
    help='output directory for pspec_boot files (default "")')

opts,args = o.parse_args(sys.argv[1:])

random.seed(0)
POL = 'I'
LST_STATS = False
DELAY = False
NGPS = 5
INJECT_SIG = 0.
SAMPLE_WITH_REPLACEMENT = True
NOISE = opts.noise
PLOT = opts.plot
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
    return n.random.random(size=size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=size))

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
for dset_count in xrange(len(min([dsets['even'],dsets['odd']]))):
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

(uvw,t1,(i,j)),d = uv.read()
(uvw,t2,(i,j)),d = uv.read()
while t1 == t2:
    (uvw,t2,(i,j)),d = uv.read()
inttime = (t2-t1)* (3600*24)
del(uv)
#inttime = 50

afreqs = freqs.take(chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

aa = a.cal.get_aa(opts.cal, freqs)
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
sep2ij, blconj, bl2sep = capo.zsa.grid2ij(aa.ant_layout)
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
antstr = 'cross'
lsts,data,flgs = {},{},{}
days = dsets.keys()
for k in days:
    if False:
        lsts[k],data[k],flgs[k] = get_data(dsets[k], antstr=antstr, polstr=POL, rmbls=rmbls, verbose=True)
    else:
        lsts[k],data[k],flgs[k] = capo.zsa.get_dict_of_uv_data(dsets[k], antstr=antstr, polstr=POL,  verbose=True,decimate=1,return_lsts=True,recast_as_array=False)
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
            data[k][bl][POL],flgs[k][bl][POL] = data[k][bl][POL][i:],flgs[k][bl][POL][i:]
    print [len(lsts[k]) for k in days]
    j = min([len(lsts[k]) for k in days])
    for k in days:
        lsts[k] = lsts[k][:j]
        for bl in data[k]:
            data[k][bl][POL],flgs[k][bl][POL] = n.array(data[k][bl][POL][:j]),n.array(flgs[k][bl][POL][:j])
else:
    for k in days:
        for bl in data[k]:
            data[k][bl][POL], flgs[k][bl][POL] = n.array(data[k][bl][POL][:]), n.array(flgs[k][bl][POL][:])
lsts = lsts.values()[0]

x = {}
xf = {}
print len(data[k][bl][POL])
print type(chans)
for k in days:
    x[k] = {}
    xf[k] = {}
    for bl in data[k]:
        print k, bl
        d = data[k][bl][POL][:,chans] * jy2T
        f = flgs[k][bl][POL][:,chans]
        if conj[bl]: d = n.conj(d)
        x[k][bl] = n.transpose(d, [1,0]) # swap time and freq axes
        xf[k][bl] = n.transpose(f,[1,0])

bls_master = x.values()[0].keys()
nbls = len(bls_master)
print 'Baselines:', nbls

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
        frp, bins = capo.frf_conv.aa_to_fr_profile(aa, bl_ij,101, bins= frbins, frpad=opts.frpad)
        timebins, firs = capo.frf_conv.frp_to_firs(frp, bins, aa.get_afreqs(), fq0=aa.get_afreqs()[101], startprms=(.001*opts.frpad,.0001))

    for k in days:
        noise_array[k] = {}
        for bl in x[k]:

            ed = noise(x[k][bls_master[0]].shape) * NOISE

            if opts.filter_noise:
                if conj[bl]: fir = n.conj(firs)
                else: fir = n.copy(firs)

                for cnt,ch in enumerate(chans):
                    fg = n.logical_not(xf[k][bl][cnt])
                    ed[cnt] = n.convolve(ed[cnt], fir[ch], mode='same')
                    w = n.convolve(fg, n.abs(n.conj(fir[ch])), mode='same')
                    ed[cnt] = n.where(fg>0, ed[cnt]/w, 1)


            if conj[bl]: ed = n.conj(ed)
            noise_array[k][bl] = n.transpose( ed.T * jy2T, [1,0])

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
#Q = {} # Create the Q's that extract power spectrum modes
#for i in xrange(nchan):
#    Q[i] = get_Q(i, nchan)
Q = [get_Q(i,nchan) for i in xrange(nchan)]

# Compute baseline auto-covariances and apply inverse to data
I,_I,_Ix = {},{},{}
C,_C,_Cx = {},{},{}
for k in days:
    I[k],_I[k],_Ix[k] = {},{},{}
    C[k],_C[k],_Cx[k] = {},{},{}
    for bl in x[k]:
        C[k][bl] = cov(x[k][bl])
        I[k][bl] = n.identity(C[k][bl].shape[0])
        U,S,V = n.linalg.svd(C[k][bl].conj())
        _C[k][bl] = n.einsum('ij,j,jk', V.T, 1./S, U.T)
        #_C[k][bl] = n.identity(_C[k][bl].shape[0])
        _I[k][bl] = n.identity(_C[k][bl].shape[0])
        _Cx[k][bl] = n.dot(_C[k][bl], x[k][bl])
        _Ix[k][bl] = x[k][bl].copy()
        if PLOT and True:
            #p.plot(S); p.show()
            p.subplot(311); capo.arp.waterfall(x[k][bl], mode='real')
            p.subplot(334); capo.arp.waterfall(C[k][bl])
            p.subplot(335); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
            p.subplot(336); capo.arp.waterfall(_C[k][bl])
            p.subplot(313); capo.arp.waterfall(_Cx[k][bl], mode='real')
            p.suptitle('%d_%d'%a.miriad.bl2ij(bl))
#            p.figure(2); p.plot(n.diag(S))
            p.show()
        


for boot in xrange(opts.nboot):
    print '%d / %d' % (boot+1,opts.nboot)
    bls = bls_master[:]
    if True: # shuffle and group baselines for bootstrapping
        if not SAMPLE_WITH_REPLACEMENT:
            random.shuffle(bls)
            bls = bls[:-5] # XXX
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
    _Cz,_Csum,_CsumQ = {},{},{}
    Csum = {}
    for k in days:
        _Iz[k],_Isum[k],_IsumQ[k] = {},{},{}
        _Cz[k],_Csum[k],_CsumQ[k] = {},{},{}
        Csum[k]={}
        for i,gp in enumerate(gps):
            _Iz[k][i] = sum([_Ix[k][bl] for bl in gp])
            _Cz[k][i] = sum([_Cx[k][bl] for bl in gp])
            _Isum[k][i] = sum([_I[k][bl] for bl in gp])
            _Csum[k][i] = sum([_C[k][bl] for bl in gp])
            Csum[k][i] = sum([C[k][bl] for bl in gp])
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
            _Czk = n.reshape(_Czk, (_Czk.shape[0]*_Czk.shape[1], _Czk.shape[2]))
            _Izk = n.reshape(_Izk, (_Izk.shape[0]*_Izk.shape[1], _Izk.shape[2]))
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
    qC = n.zeros((nchan,_Cz.values()[0].values()[0].shape[1]), dtype=n.complex)
    Q_Iz = {}
    Q_Cz = {}
    for cnt1,k1 in enumerate(days):
        for k2 in days[cnt1:]:
            if not Q_Iz.has_key(k2): Q_Iz[k2] = {}
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
                        qC += n.conj(_Cz[k1][bl1]) * _Cz[k2][bl2]
                    else: # brute force with Q to ensure normalization
                        #_qI = n.array([_Iz[k1][bl1].conj() * n.dot(Q[i], _Iz[k2][bl2]) for i in xrange(nchan)])
                        #_qC = n.array([_Cz[k1][bl1].conj() * n.dot(Q[i], _Cz[k2][bl2]) for i in xrange(nchan)])
                        if not Q_Iz[k2].has_key(bl2): Q_Iz[k2][bl2] = [n.dot(Q[i], _Iz[k2][bl2]) for i in xrange(nchan)]
                        if not Q_Cz[k2].has_key(bl2): Q_Cz[k2][bl2] = [n.dot(Q[i], _Cz[k2][bl2]) for i in xrange(nchan)]
                        _qI = n.array([_Iz[k1][bl1].conj() * Q_Iz[k2][bl2][i] for i in xrange(nchan)])
                        qI += n.sum(_qI, axis=1)
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

    if PLOT:
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

    outfile = 'pspec_boot%04d.npz'%(boot)
    if not opts.output == '':
        outfile =opts.output+'/'+outfile
    print "Writing", outfile
    n.savez(outfile, kpl=kpl, scalar=scalar, times=n.array(lsts),
        pk_vs_t=pC, err_vs_t=1./cnt, temp_noise_var=var, nocov_vs_t=pI, 
        afreqs=afreqs, chans=chans, frpad=opts.frpad,
        cmd=' '.join(sys.argv))

