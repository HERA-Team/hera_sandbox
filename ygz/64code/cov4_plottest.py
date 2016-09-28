#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo
import glob, optparse, sys, random, itertools

# sample run python cov4_plottest.py --window=none --sep=sep0,1 --sepd=sep0,1 --plot -b 1 -c 95_115 -p I -C psa6240_v003
o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
o.add_option('-b', '--nboot', type='int', default=20,
    help='Number of bootstraps.  Default is 20')
o.add_option('--plot', action='store_true',
    help='Generate plots')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
o.add_option('--sep', default='sep0,1', action='store',  #-1,1 is like 0_38
    help='Which separation type?')
o.add_option('--sepd', default='sep-1,1', action='store',  #-1,1 is like 0_38
    help='Which separation type?')
o.add_option('--loss', action='store',
    help='In signal loss mode to measure the signal loss. Uses default data in my path. Give it the path to the simulated signal data. Assumes ends in ')
o.add_option('--level', type='float', default=-1.0,
    help='Scalar to multiply the default signal level for simulation runs.')
o.add_option('--rmbls', action='store',
    help='List of baselines, in miriad format, to remove from the power spectrum analysis.')

opts,args = o.parse_args(sys.argv[1:])
 #args is calfile
 #################################################
SEP, SEPD = opts.sep, opts.sepd
DicT = {'sep0,1_sep-1,1':0.032557, 'sep0,1_sep0,1':0.,'sep-1,1_sep-1,1':0.,'sep1,1_sep1,1':0.}
DelT = DicT[SEP+'_'+SEPD]
print 'DelT=', DelT
###################################################
random.seed(0)
POL = 'I'
VERBOSE = True
LST_STATS = False
DELAY = False
NGPS = 5
INJECT_SIG = 0.
SAMPLE_WITH_REPLACEMENT = True
NOISE = .0
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
        #if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            bl = a.miriad.ij2bl(i,j)
            if bl in rmbls: continue
            lst = uv['lst']
            #print lst
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

def cov2(m1,m2):
    X1 = n.array(m1, ndmin=2, dtype=n.complex); X2 = n.array(m2, ndmin=2, dtype=n.complex)
    X1 -= X1.mean(axis=1)[(slice(None),n.newaxis)]; X2 -= X2.mean(axis=1)[(slice(None),n.newaxis)]
    N = X1.shape[1]
    fact = float(N - 1)
    return (n.dot(X1, X2.T.conj()) / fact).squeeze()

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

dsets1 = {
    'even': glob.glob('/Users/yunfanzhang/local/DATA64/ali_et_al_2015_apj_data/even/'+SEP+'/*242.[3456]*uvGL'),
    'odd' : glob.glob('/Users/yunfanzhang/local/DATA64/ali_et_al_2015_apj_data/odd/'+SEP+'/*243.[3456]*uvGL'),
}
dsets2 = {
    'even': glob.glob('/Users/yunfanzhang/local/DATA64/ali_et_al_2015_apj_data/even/'+SEPD+'/*242.[3456]*uvGL'),
    'odd' : glob.glob('/Users/yunfanzhang/local/DATA64/ali_et_al_2015_apj_data/odd/'+SEPD+'/*243.[3456]*uvGL'),
}
#for i in xrange(10): dsets[i] = glob.glob('lstbinX%d/%s/lst.24562[45]*.[3456]*.uvAL'%(i,SEP))
#print dsets
if opts.loss:
    dsets = {
    'even': glob.glob('/Users/sherlock/projects/paper/analysis/psa64/signal_loss/data/even/*242.[3456]*uvALG'),
    'odd' : glob.glob('/Users/sherlock/projects/paper/analysis/psa64/signal_loss/data/odd/*243.[3456]*uvALG'),
}

WINDOW = opts.window
uv = a.miriad.UV(dsets1.values()[0][0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
del(uv)

afreqs = freqs.take(chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

aa = a.cal.get_aa(opts.cal, n.array([.150]))
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
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
kpl = etas * capo.pspec.dk_deta(z) #111

if True:
    bm = n.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq) * 2.35 # correction for beam^2
    scalar = capo.pspec.X2Y(z) * bm * B
else: scalar = 1
if not DELAY:
    # XXX this is a hack
    if WINDOW == 'hamming': scalar /= 3.67
    elif WINDOW == 'blackman-harris': scalar /= 5.72

sys.stdout.flush()

# acquire the data
#antstr = '41_49,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_47,3_25,1_58,35_61,42_63,2_33'
antstr = 'cross'
lsts1,data1,flgs1 = {},{},{}; lsts2,data2,flgs2 = {},{},{}
days = dsets1.keys()
for k in days:
    print 'Reading', SEP,k
    lsts1[k],data1[k],flgs1[k] = get_data(dsets1[k], antstr=antstr, polstr=POL, rmbls=rmbls, verbose=False)
    print 'Reading', SEPD,k
    lsts2[k],data2[k],flgs2[k] = get_data(dsets2[k], antstr=antstr, polstr=POL, rmbls=rmbls, verbose=False)

## format  lsts[day]=array of lst
cnt,var = n.ones_like(lsts1.values()[0]), n.ones_like(lsts1.values()[0])

##################################################################################
#for k in days: lsts2[k][:] = lsts2[k][:] - DelT
#MAGIC!!!!!!!!!!!!!!!!!!!
##################################################################################

# Align data sets in LST to all start at the lastest start lstmax
lstmax = max([lsts1[k][0] for k in days]+[lsts2[k][0] for k in days])    #latest start from all days
for k in days:
    print 'aligning lst'
    #print len(lsts1[k]),len(lsts2[k])
    for i1 in xrange(len(lsts1[k])):
        if lsts1[k][i1] >= lstmax - .001: break
    # so now the ith time of day k is
    lsts1[k] = lsts1[k][i1:]
    for i2 in xrange(len(lsts2[k])):
        if lsts2[k][i2] >= lstmax - .001: break
    lsts2[k] = lsts2[k][i2:]
    for bl in data1[k]:
        data1[k][bl],flgs1[k][bl] = n.array(data1[k][bl][i1:]),n.array(flgs1[k][bl][i1:])
    for bl in data2[k]:
        data2[k][bl],flgs2[k][bl] = n.array(data2[k][bl][i2:]),n.array(flgs2[k][bl][i2:])
    #print len(lsts1[k]),len(lsts2[k])
j = min([len(lsts1[k]) for k in days]+[len(lsts2[k]) for k in days])
for k in days:
    lsts1[k] = lsts1[k][:j]; lsts2[k] = lsts2[k][:j]
    # print len(lsts1[k]),len(lsts2[k])
    # print lsts1[k].shape
    # print lsts2[k].shape
    for bl in data1[k]:
        data1[k][bl],flgs1[k][bl] = n.array(data1[k][bl][:j]),n.array(flgs1[k][bl][:j])
    for bl in data2[k]:
        data2[k][bl],flgs2[k][bl] = n.array(data2[k][bl][:j]),n.array(flgs2[k][bl][:j])
lsts = lsts1.values()[0]
#print lsts1
#print lsts2
###################################################################################
x1,x2 = {},{}
for k in days:
    x1[k],x2[k] = {},{}
    for bl in data1[k]:
        if VERBOSE: print k, bl
        d = data1[k][bl][:,chans] * jy2T
        if conj[bl]: d = n.conj(d)
        x1[k][bl] = n.transpose(d, [1,0]) # swap time and freq axes
        if VERBOSE: print x1[k][bl].shape
    for bl in data2[k]:
        if VERBOSE: print k, bl
        d = data2[k][bl][:,chans] * jy2T
        if conj[bl]: d = n.conj(d)
        x2[k][bl] = n.transpose(d, [1,0]) # swap time and freq axes
        if VERBOSE: print x2[k][bl].shape


bls_master1 = x1.values()[0].keys(); bls_master2 = x2.values()[0].keys()
nbls1 = len(bls_master1); nbls2 = len(bls_master2)
print 'number of baselines:', nbls1, nbls2

Q = [get_Q(i,nchan) for i in xrange(nchan)]


##########################################################################
# Compute baseline auto-covariances and apply inverse to data
I,_I,_Ix1,_Ix2 = {},{},{},{}
C,_C,_Cx1,_Cx2 = {},{},{},{}
C1,C2,_C1,_C2 = {},{},{},{}
I1,I2,_I1,_I2 = {},{},{},{}
for k in days:
    print len(x1[k]), len(x2[k])
    I[k],_I[k],_Ix1[k],_Ix2[k] = {},{},{},{}
    C[k],_C[k],_Cx1[k],_Cx2[k]  = {},{},{},{}
    C1[k],_C1[k],C2[k],_C2[k]  = {},{},{},{}
    I1[k],_I1[k],I2[k],_I2[k]  = {},{},{},{}
    for bl1 in x1[k]:
        C1[k][bl1] = cov(x1[k][bl1])
        I1[k][bl1] = n.identity(C1[k][bl1].shape[0])
        U1,S1,V1 = n.linalg.svd(C1[k][bl1].conj()) ; _C1[k][bl1] = n.einsum('ij,j,jk', V1.T, 1./S1, U1.T)
        _I1[k][bl1] = n.identity(_C1[k][bl1].shape[0])
        _Cx1[k][bl1] = n.dot(_C1[k][bl1], x1[k][bl1]); _Ix1[k][bl1] = x1[k][bl1].copy()
    for bl2 in x2[k]:
        C2[k][bl2] = cov(x2[k][bl2])
        I2[k][bl2] = n.identity(C2[k][bl2].shape[0])
        #print bl1,bl2
        #print n.linalg.eigvals(C[k][blt])
        U2,S2,V2 = n.linalg.svd(C2[k][bl2].conj()) ; _C2[k][bl2] = n.einsum('ij,j,jk', V2.T, 1./S2, U2.T)
        _I2[k][bl2] = n.identity(_C2[k][bl2].shape[0])
        _Cx2[k][bl2] = n.dot(_C2[k][bl2], x2[k][bl2]); _Ix2[k][bl2] = x2[k][bl2].copy()
        #mport IPython; IPython.embed()
        if PLOT and bl2 == a.miriad.ij2bl(0,26):
            #p.plot(S); p.show()
            p.subplot(311); capo.arp.waterfall(x2[k][bl2], mode='real')
            p.subplot(323); capo.arp.waterfall(C2[k][bl2])
            p.subplot(324); p.plot(n.einsum('ij,jk',n.diag(S2),V2).T.real)
            p.subplot(313); capo.arp.waterfall(_Cx2[k][bl2], mode='real')
            p.suptitle('%d_%d'%a.miriad.bl2ij(bl2))
#            p.figure(2); p.plot(n.diag(S))
            p.show()
##########################################################################
#C2,_C2,_Cx2,I2,_I2,_Ix2 = C1,_C1,_Cx1,I1,_I1,_Ix1# for test against v002
##########################################################################

#import IPython; IPython.embed

for boot in xrange(opts.nboot):
    print '%d / %d' % (boot+1,opts.nboot)
    bls1 = bls_master1[:]; bls2 = bls_master2[:]
    if True: # shuffle and group baselines for bootstrapping
        if not SAMPLE_WITH_REPLACEMENT:
            random.shuffle(bls1); random.shuffle(bls2)
            bls1 = bls1[:-5] ;bls2 = bls2[:-5]
        else: # sample with replacement
            bls1 = [random.choice(bls1) for bl in bls1]#; bls2 = [random.choice(bls2) for bl in bls2]
        gps1 = [bls1[i::NGPS] for i in range(NGPS)]#; gps2 = [bls2[i::NGPS] for i in range(NGPS)]
        gps1 = [[random.choice(gp) for bl in gp] for gp in gps1]#; gps2 = [[random.choice(gp) for bl in gp] for gp in gps2]
    ############################################################
    gps2 = gps1 #for test against v002
    #print gps1
    #print gps2
    ############################################################
    bls1 = [bl for gp in gps1 for bl in gp]; bls2 = [bl for gp in gps2 for bl in gp]
    #print '\n'.join([','.join(['%d_%d'%a.miriad.bl2ij(bl) for bl in gp]) for gp in gps1])
    #print '\n'.join([','.join(['%d_%d'%a.miriad.bl2ij(bl) for bl in gp]) for gp in gps2])
    #import IPython; IPython.embed()

    _Iz1,_Iz2,_Isum1,_IsumQ1,_Isum2,_IsumQ2 = {},{},{},{},{},{}
    _Cz1,_Cz2,_Csum1,_CsumQ1,_Csum2,_CsumQ2  = {},{},{},{},{},{}
    for k in days:
        _Iz1[k],_Iz2[k],_Isum1[k],_IsumQ1[k],_Isum2[k],_IsumQ2[k] = {},{},{},{},{},{}
        _Cz1[k],_Cz2[k],_Csum1[k],_CsumQ1[k],_Csum2[k],_CsumQ2[k] = {},{},{},{},{},{}
        #for i,gp in enumerate(gps1):
        #print gps1
        #print gps2
        #for i,gpt in enumerate(itertools.product(gps1, gps2)):   #must use index in dict
        for i, gp1 in enumerate(gps1):
            _Iz1[k][i] = sum([_Ix1[k][bl] for bl in gp1])
            _Cz1[k][i] = sum([_Cx1[k][bl] for bl in gp1])
            _Isum1[k][i] = sum([_I1[k][bl] for bl in gp1])
            _Csum1[k][i] = sum([_C1[k][bl] for bl in gp1])
            _IsumQ1[k][i] = {}
            _CsumQ1[k][i] = {}
            #if DELAY: # this is much faster
            #    _Iz1[k][(i,j)] = n.fft.fftshift(n.fft.ifft(window*_Iz1[k][gpt], axis=0), axes=0)
            #    _Cz1[k][(i,j)] = n.fft.fftshift(n.fft.ifft(window*_Cz1[k][gpt], axis=0), axes=0)
            for ch in xrange(nchan): # XXX this loop makes computation go as nchan^3
                _IsumQ1[k][i][ch] = n.dot(_Isum1[k][i], Q[ch])
                _CsumQ1[k][i][ch] = n.dot(_Csum1[k][i], Q[ch])
        for j, gp2 in enumerate(gps2):
            #gpt = (gp1,gp2)
            #print gpt.shape
            #for blt in itertools.product(*gpt): print blt
            #print [_Ix1[k][blt] for blt in itertools.product(*gpt)]
            _Iz2[k][j] = sum([_Ix2[k][bl] for bl in gp2])
            _Cz2[k][j] = sum([_Cx2[k][bl] for bl in gp2])
            _Isum2[k][j] = sum([_I2[k][bl] for bl in gp2])
            _Csum2[k][j] = sum([_C2[k][bl] for bl in gp2])
            _IsumQ2[k][j] = {}
            _CsumQ2[k][j] = {}

            #if DELAY: # this is much faster
            #    _Iz1[k][(i,j)] = n.fft.fftshift(n.fft.ifft(window*_Iz1[k][gpt], axis=0), axes=0)
            #    _Cz1[k][(i,j)] = n.fft.fftshift(n.fft.ifft(window*_Cz1[k][gpt], axis=0), axes=0)
            #    _Iz2[k][(i,j)] = n.fft.fftshift(n.fft.ifft(window*_Iz2[k][gpt], axis=0), axes=0)
            #    _Cz2[k][(i,j)] = n.fft.fftshift(n.fft.ifft(window*_Cz2[k][gpt], axis=0), axes=0)
                # XXX need to take fft of _Csum, _Isum here
            for ch in xrange(nchan): # XXX this loop makes computation go as nchan^3
                _IsumQ2[k][j][ch] = n.dot(_Isum2[k][j], Q[ch])
                _CsumQ2[k][j][ch] = n.dot(_Csum2[k][j], Q[ch])
        if PLOT:
            gps = gps1
            _Cz = _Cz1
            NGPS = len(gps)
            _Csumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            _Isumk = n.zeros((NGPS,nchan,NGPS,nchan), dtype=n.complex)
            for i in xrange(len(gps)): _Isumk[i,:,i,:] = _Isum1[k][i]
            _Isumk.shape = (NGPS*nchan, NGPS*nchan)
            #_Isum[k] = _Isumk
            for i in xrange(len(gps)): _Csumk[i,:,i,:] = _Csum1[k][i]
            _Csumk.shape = (NGPS*nchan, NGPS*nchan)
            #_Csum[k] = _Csumk
            _Czk = n.array([_Cz[k][i] for i in _Cz[k]])
            _Czk = n.reshape(_Czk, (_Czk.shape[0]*_Czk.shape[1], _Czk.shape[2]))
            p.subplot(211); capo.arp.waterfall(_Czk, mode='real')
            p.subplot(223); capo.arp.waterfall(_Csumk)
            p.subplot(224); capo.arp.waterfall(cov(_Czk))
            p.show()
################################################################################################
    #import IPython; IPython.embed()

    FI = n.zeros((nchan,nchan), dtype=n.complex)
    FC = n.zeros((nchan,nchan), dtype=n.complex)
    qI = n.zeros((nchan,_Iz1.values()[0].values()[0].shape[1]), dtype=n.complex)
    qC = n.zeros((nchan,_Cz1.values()[0].values()[0].shape[1]), dtype=n.complex)
    Q_Iz1,Q_Iz2 = {},{}
    Q_Cz1,Q_Cz2 = {},{}
    for cnt1,k1 in enumerate(days):
        for k2 in days[cnt1:]:
            #if not Q_Iz1.has_key(k2): Q_Iz1[k2] = {}
            #if not Q_Cz1.has_key(k2): Q_Cz1[k2] = {}
            if not Q_Iz2.has_key(k2): Q_Iz2[k2] = {}
            if not Q_Cz2.has_key(k2): Q_Cz2[k2] = {}
            for bl1,bl2 in itertools.product(_Cz1[k1],_Cz2[k2]):
                #for bl2 in _Cz2[k2]:
                #if k1 == k2 and bl1 == bl2: continue # this results in a significant bias
                #################################################################
                if k1 == k2 or bl1 == bl2: continue
                #if k1 == k2: continue
                #if bl1 == bl2: continue # also a significant noise bias
                #################################################################
                if VERBOSE: print k1, k2, bl1, bl2
                # if PLOT and False:
                #     p.subplot(231); capo.arp.waterfall(C[m], drng=3)
                #     p.subplot(232); capo.arp.waterfall(_C[m], drng=3)
                #     p.subplot(233); capo.arp.waterfall(n.dot(C[m],_C[m]), drng=3)
                #     p.subplot(234); p.semilogy(S)
                #     p.subplot(236); capo.arp.waterfall(V, drng=3)
                #     p.show()
                #     p.subplot(311); capo.arp.waterfall(x[m], mode='real', mx=5, drng=10); p.colorbar(shrink=.5)
                #     p.subplot(312); capo.arp.waterfall(_Cx, mode='real'); p.colorbar(shrink=.5)
                #     p.subplot(313); capo.arp.waterfall(_Ix, mode='real'); p.colorbar(shrink=.5)
                #     p.show()
                if False: # use ffts to do q estimation fast
                    print 'nothing'
                #     qI += n.conj(_Iz[k1][bl1]) * _Iz[k2][bl2]
                #     qC += n.conj(_Cz[k1][bl1]) * _Cz[k2][bl2]
                else: # brute force with Q to ensure normalization
                    #_qI = n.array([_Iz[k1][bl1].conj() * n.dot(Q[i], _Iz[k2][bl2]) for i in xrange(nchan)])
                    #_qC = n.array([_Cz[k1][bl1].conj() * n.dot(Q[i], _Cz[k2][bl2]) for i in xrange(nchan)])
                    blt = (bl1,bl2)
                    if not Q_Iz2[k2].has_key(bl2): Q_Iz2[k2][bl2] = [n.dot(Q[c], _Iz2[k2][bl2]) for c in xrange(nchan)]
                    if not Q_Cz2[k2].has_key(bl2): Q_Cz2[k2][bl2] = [n.dot(Q[c], _Cz2[k2][bl2]) for c in xrange(nchan)]
                    _qI = n.array([_Iz1[k1][bl1].conj() * Q_Iz2[k2][bl2][c] for c in xrange(nchan)])
                    qI += n.sum(_qI, axis=1)
                    _qC = n.array([_Cz1[k1][bl1].conj() * Q_Cz2[k2][bl2][c] for c in xrange(nchan)])
                    qC += n.sum(_qC, axis=1)
                if DELAY: # by taking FFT of CsumQ above, each channel is already i,j separated
                    FI += n.conj(_IsumQ1[k1][bl1]) * _IsumQ2[k2][bl2]
                    FC += n.conj(_CsumQ1[k1][bl1]) * _CsumQ2[k2][bl2]
                else:
                    for i in xrange(nchan):
                        for j in xrange(nchan):
                            #print 'i,j=',i,j
                            FI[i,j] += n.einsum('ij,ji', _IsumQ1[k1][bl1][i], _IsumQ2[k2][bl2][j])
                            FC[i,j] += n.einsum('ij,ji', _CsumQ1[k1][bl1][i], _CsumQ2[k2][bl2][j])
   # import IPython; IPython.embed()
    # if PLOT:
    #     p.subplot(121); capo.arp.waterfall(FC, drng=4)
    #     p.subplot(122); capo.arp.waterfall(FI, drng=4)
    #     p.show()
    #
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
    ####################################################################
    #print n.linalg.eigvals(FC_o)
    #import IPython; IPython.embed
    #FC_o = n.abs(FC_o)
    ####################################################################
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

    # if PLOT:
    #     p.subplot(411); capo.arp.waterfall(qC, mode='real'); p.colorbar(shrink=.5)
    #     p.subplot(412); capo.arp.waterfall(pC, mode='real'); p.colorbar(shrink=.5)
    #     p.subplot(413); capo.arp.waterfall(qI, mode='real'); p.colorbar(shrink=.5)
    #     p.subplot(414); capo.arp.waterfall(pI, mode='real'); p.colorbar(shrink=.5)
    #     p.show()
    #
    # if PLOT:
    #     p.plot(kpl, n.average(pC.real, axis=1), 'b.-')
    #     p.plot(kpl, n.average(pI.real, axis=1), 'k.-')
    #     p.show()

    print 'Writing pspec_boot%04d.npz' % boot
    n.savez('boot/'+str(SEP)+'_'+str(SEPD)+'_boot_%04d.npz'%boot, kpl=kpl, scalar=scalar, times=n.array(lsts),
        pk_vs_t=pC, err_vs_t=1./cnt, temp_noise_var=var, nocov_vs_t=pI,
        cmd=' '.join(sys.argv))


