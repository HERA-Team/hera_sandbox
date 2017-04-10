#! /usr/bin/env python

import aipy as a, numpy as n, pylab as p, capo, capo.frf_conv as fringe
import glob, optparse, sys, random, os
import capo.zsa as zsa
import oqe

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True)
o.add_option('-b', '--nboot', type='int', default=20,
    help='Number of bootstraps.  Default is 20.')
o.add_option('--plot', action='store_true',
    help='Generate plots')
# o.add_option('--ngp', action='store_true',
#     dest='NOGROUP')
# o.add_option('--ngp', action='store_true',
#     dest='NOGROUP')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
o.add_option('--sep', default='sep0,1', action='store',
    help='Which separation directory to use for signal loss data.')
o.add_option('--sepd', default='sep-1,1', action='store',  #-1,1 is like 0_38
    help='Second separation type')
o.add_option('--noise_only', action='store_true',
    help='Instead of injecting noise, Replace data with noise')
o.add_option('--same', action='store_true',
    help='Noise is the same for all baselines.')
o.add_option('--diff', action='store_true',
    help='Noise is different for all baseline.')
o.add_option('--frf', action='store_true',
    help='FRF noise.')
o.add_option('--lmode',type='int', default=None,
    help='Eigenvalue mode of C (in decreasing order) to be the minimum value used in C^-1')
o.add_option('--lbins',type='int', default=None,
    help='number of lst_bins to grid into, if None use lst_align')
o.add_option('--output', type='string', default='./boot',
    help='Output directory for pspec_boot files (default "")')

opts,args = o.parse_args(sys.argv[1:])

#Basic parameters
n.random.seed(0)
POL = opts.pol 
LST_STATS = False
DELAY = False
NGPS = 5 #number of groups to break the random sampled bls into
PLOT = opts.plot
LMODE = opts.lmode

 #################################################
SEP, SEPD = opts.sep, opts.sepd
DicT = {'sep0,1_sep-1,1':0.0535, 'sep0,1_sep0,1':0.,'sep-1,1_sep-1,1':0.,'sep1,1_sep1,1':0.}
DicS = {'sep0,1_sep-1,1':0.46000}
DelT = DicT[SEP+'_'+SEPD]*2*n.pi*(24*3600./a.const.sidereal_day)   #sidereal day?
print 'DelT=', DelT, 'radians'
EQUIV = (SEP == SEPD)
if EQUIV: s_factor = 1
else:
    s_factor = DicS[SEP+'_'+SEPD]
NOGROUP = True
if EQUIV: NOGROUP=False
print NOGROUP, SEP, SEPD
###################################################

### FUNCTIONS ###

def frf(shape): #FRF NOISE
    shape = shape[1]*2,shape[0] #(2*times,freqs)
    dij = oqe.noise(size=shape)
    wij = n.ones(shape,dtype=bool) #XXX flags are all true (times,freqs)
    #dij and wij are (times,freqs)
    _d,_w,_,_ = fringe.apply_frf(aa,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)
    _d = n.transpose(_d)
    _d = _d[:,shape[0]/4:shape[0]/2+shape[0]/4]
    return _d

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    X = n.array(m, ndmin=2, dtype=n.complex)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1) #normalization
    return (n.dot(X, X.T.conj()) / fact).squeeze()

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
#print dsets

#Get uv file info
WINDOW = opts.window
uv = a.miriad.UV(dsets.values()[0][0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sdf = uv['sdf']
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
inttime = uv['inttime']
afreqs = freqs.take(chans)
nchan = chans.size
fq = n.average(afreqs)
z = capo.pspec.f2z(fq)

aa = a.cal.get_aa(opts.cal, afreqs)
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
sep2ij, blconj, bl2sep = capo.zsa.grid2ij(aa.ant_layout)
#import IPython; IPython.embed()
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
kpl = etas * capo.pspec.dk_deta(z) 
if True:
    bm = n.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq) * 2.35 #correction for beam^2
    scalar = capo.pspec.X2Y(z) * bm * B / s_factor
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
data_dict = {}
flg_dict = {}
conj_dict = {}
antstr = 'cross' #?
#,blconj,_ = zsa.grid2ij(aa.ant_layout)
days = dsets.keys()
lsts,data,flgs = {},{},{}
def tup2usc(tup):
    return str(tup[0])+'_'+str(tup[1])
def prostr(sepstr):
    return sepstr[3:]

sepls = sep2ij[prostr(SEP)]; sepdls = sep2ij[prostr(SEPD)]
sep_avail, sepd_avail = [],[]

for k in days:
    lsts[k],data[k],flgs[k] = capo.miriad.read_files(dsets[k], antstr=antstr, polstr=POL, verbose=True)
    lsts[k] = n.array(lsts[k]['lsts'])
    for bl in data[k]:
        if tup2usc(bl) in sepls:# and k == days[0]:
            #sep = SEP
            sep_avail.append(bl)
            if EQUIV: sepd_avail.append(bl)
        elif tup2usc(bl) in sepdls:# and k == days[0]:
            #sep = SEPD
            sepd_avail.append(bl)
        else:
            #print 'skipping', bl, tup2usc(bl)
            continue
        d = n.array(data[k][bl][POL])[:,chans] * jy2T  #extract frequency range
        flg = n.array(flgs[k][bl][POL])[:,chans]
        key = (k,bl,POL)
        data_dict[key] = d
        flg_dict[key] = n.logical_not(flg)
        conj_dict[key[1]] = conj[bl]
keys = data_dict.keys()
bls_master = []
for key in keys: #populate list of baselines
    if key[0] == keys[0][0]: bls_master.append(key[1])
print 'Baselines:', len(bls_master)
#import IPython; IPython.embed()
#Align and create dataset
ds = oqe.DataSet(lmode=LMODE)
if opts.lbins:
    #import IPython; IPython.embed()
    lsts_dict = lsts
    for k in data_dict.keys():
        #print 'gridding lst for', k
        data_dict[k],flg_dict[k],lsts = oqe.lst_grid(lsts_dict[k[0]], data_dict[k], wgts=flg_dict[k], lstbins=opts.lbins)
else:
    #import IPython; IPython.embed()
    if not EQUIV: 
        lststemp = lsts[lsts.keys()[0]]
        dshift = int(DelT/(lststemp[1]-lststemp[0]))
        print 'shifting by ', dshift
        #import IPython; IPython.embed()
        data_dict,flg_dict, lsts = oqe.lst_shift(sep_avail, sepd_avail, data_dict, flg_dict, lsts, dshift)
    #import IPython; IPython.embed()
    inds = oqe.lst_align(lsts)
    data_dict,flg_dict,lsts = oqe.lst_align_data(inds,dsets=data_dict,wgts=flg_dict,lsts=lsts)
    #the lsts given is a dictionary with 'even','odd', etc., but the lsts returned is one array
    
#import IPython; IPython.embed()
#If data is replaced by noise
if opts.noise_only:
    if opts.same == None and opts.diff == None: 
        print 'Need to specify if noise is the same on all baselines (--same) or different (--diff)'
        sys.exit()
    #Prep FRF Stuff
    ij = bls_master[0] #ij = (1,4)
    if blconj[a.miriad.ij2bl(ij[0],ij[1])]: #makes sure FRP will be the same whether bl is a conjugated one or not
        if ij[0] < ij[1]: temp = (ij[1],ij[0]); ij=temp
    timelen = data_dict[keys[0]].shape[0] 
    bins = fringe.gen_frbins(inttime)
    frp, bins = fringe.aa_to_fr_profile(aa, ij, len(afreqs)/2, bins=bins)
    timebins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(), fq0=aa.get_freqs()[len(afreqs)/2])
    fir = {(ij[0],ij[1],POL):firs}
    if opts.same and opts.frf: NOISE = frf((len(chans),timelen)) #same noise on all bls
    if opts.same and opts.frf == None: NOISE = oqe.noise((len(chans),timelen))
    for key in data_dict:
        if opts.same: thing = NOISE.T
        if opts.diff and opts.frf: thing = frf((len(chans),timelen)).T
        if opts.diff and opts.frf == None: thing = oqe.noise((len(chans),timelen)).T
        if blconj[a.miriad.ij2bl(key[1][0],key[1][1])]: data_dict[key] = n.conj(thing)
        else: data_dict[key] = thing
        flg_dict[key] = n.ones_like(data_dict[key])

#Set data
ds.set_data(dsets=data_dict,conj=conj_dict,wgts=flg_dict)

#Get some statistics
if LST_STATS:
    #collect some metadata from the lst binning process
    cnt, var = {}, {}
    for filename in dsets.values()[0]:
        print 'Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, '64_49', POL) #XXX
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            bl = '%d,%d,%d' % (i,j,uv['pol'])
            cnt[bl] = cnt.get(bl, []) + [uv['cnt']]
            var[bl] = var.get(bl, []) + [uv['var']]
    cnt = n.array(cnt.values()[0]) #all baselines should be the same
    var = n.array(var.values()[0]) #all baselines should be the same
else: cnt,var = n.ones_like(lsts), n.ones_like(lsts)

if PLOT and True:
    for key in keys:
        p.subplot(311); capo.arp.waterfall(ds.x[key], mode='real')
        p.colorbar()
        p.title('Data x')
        p.subplot(323); capo.arp.waterfall(ds.C(key))
        p.colorbar()
        p.title('C')
        U,S,V = n.linalg.svd(ds.C(key).conj()) 
        if LMODE: S += S[LMODE-1]
        p.subplot(324); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
        p.subplot(313); capo.arp.waterfall(n.dot(ds.iC(key),ds.x[key]), mode='real')#,drng=6000,mx=3000)
        p.colorbar()
        p.title('C^-1 x')
        p.suptitle(key)
        p.tight_layout()
        p.show()
        p.semilogy(S);p.title('Eigenvalues');p.show()

#Bootstrapping        
for boot in xrange(opts.nboot):
    print 'Bootstrap %d / %d' % (boot+1,opts.nboot)
 
    if NOGROUP:
        #gp1 = np.random.choice(sep)
        gps = [sep_avail, sepd_avail]
        print 'NOGROUP'
        print len(gps), len(gps[0]), len(gps[1]), len(bls_master)
        newkeys,dsC = ds.group_data(keys,gps)
        newkeys,dsI = ds.group_data(keys,gps,use_cov=False)
        #import IPython; IPython.embed()
    elif True: 
        print 'shuffle and group baselines for bootstrapping'
        gps = ds.gen_gps(bls_master, ngps=NGPS)
        #import IPython; IPython.embed()
        newkeys,dsC = ds.group_data(keys,gps)
        newkeys,dsI = ds.group_data(keys,gps,use_cov=False)
    else: #no groups (slower)
        newkeys = [random.choice(keys) for key in keys] #sample w/replacement for bootstrapping
        dsI,dsC = ds,ds #identity and covariance case dataset is the same
    

    if PLOT and True:
        for newkey in newkeys:
            p.subplot(311); capo.arp.waterfall(dsC.x[newkey], mode='real')
            p.colorbar()
            p.title('Data x')
            p.subplot(323); capo.arp.waterfall(dsC.C(newkey))
            p.colorbar()
            p.title('C')
            U,S,V = n.linalg.svd(dsC.C(newkey).conj())
            p.subplot(324); p.plot(n.einsum('ij,jk',n.diag(S),V).T.real)
            p.subplot(313); capo.arp.waterfall(n.dot(dsC.iC(newkey),dsC.x[newkey]), mode='real')#,drng=6000,mx=3000)
            p.colorbar()
            p.title('C^-1 x')
            p.suptitle(newkey)
            p.tight_layout()
            p.show()
            p.semilogy(S);p.title('Eigenvalues');p.show()

    #OQE Stuff
    FI = n.zeros((nchan,nchan), dtype=n.complex)
    FC = n.zeros((nchan,nchan), dtype=n.complex) #{}
    qI = n.zeros((nchan,data_dict[key].shape[0]), dtype=n.complex)
    qC = n.zeros((nchan,data_dict[key].shape[0]), dtype=n.complex)
    for k,key1 in enumerate(newkeys):
        #print '   ',k+1,'/',len(keys)
        for key2 in newkeys[k:]:
            if key1[0] == key2[0] or key1[1] == key2[1]:
                continue #don't do even w/even or bl w/same bl
            else:
                print key1,key2 
                qC += dsC.q_hat(key1,key2,cov_flagging=False)
                qI += dsI.q_hat(key1,key2,use_cov=False,cov_flagging=False)
                FC += dsC.get_F(key1,key2,cov_flagging=False)
                FI += dsI.get_F(key1,key2,use_cov=False,cov_flagging=False)
                """
                qC += dsC.q_hat(key1,key2,cov_flagging=True)
                qI += dsI.q_hat(key1,key2,use_cov=False,cov_flagging=True) 
                FC_single = dsC.get_F(key1,key2,cov_flagging=True)
                FI += dsI.get_F(key1,key2,use_cov=False,cov_flagging=True)    
                for Fkey in FC_single: #XXX 
                    # keys of FC_single are (k1,m1,k2,m2), and there are different m's if using different times for covariance matrices
                    # need to add up F's for same m1 and m2's, but then what do I use for k1 and k2??
                    # p_hat excepts keys of form (k1,m1,k2,m2) and uses k's for weights, so I need to keep those separate?
                    try: FC[(('even', 0),Fkey[1],('odd', 1),Fkey[3])] += FC_single[Fkey]
                    except: FC[(('even', 0),Fkey[1],('odd', 1),Fkey[3])] = FC_single[Fkey]
                """
    #import IPython;IPython.embed()
    MC,WC = dsC.get_MW(FC,mode='L^-1') #Cholesky decomposition
    MI,WI = dsI.get_MW(FI,mode='I')
    pC = dsC.p_hat(MC,qC,scalar=scalar)
    pI = dsI.p_hat(MI,qI,scalar=scalar)
    print 'pI ~ ', n.median(pI)
    print 'pC ~ ', n.median(pC)
 
    if PLOT:
        p.subplot(121); capo.arp.waterfall(FC, drng=4)
        p.title('FC')
        p.subplot(122); capo.arp.waterfall(FI, drng=4)
        p.title('FI')
        p.show()
    
    if PLOT:
        p.subplot(411); capo.arp.waterfall(qC, mode='real'); p.colorbar(shrink=.5); p.title('qC')
        p.subplot(412); capo.arp.waterfall(pC, mode='real'); p.colorbar(shrink=.5); p.title('pC')
        p.subplot(413); capo.arp.waterfall(qI, mode='real'); p.colorbar(shrink=.5); p.title('qI')
        p.subplot(414); capo.arp.waterfall(pI, mode='real'); p.colorbar(shrink=.5); p.title('pI')
        p.show()

    if PLOT:
        p.plot(kpl, n.average(pC.real, axis=1), 'b.-', label='pC')
        p.plot(kpl, n.average(pI.real, axis=1), 'k.-', label='pI')
        p.legend()
        p.show()

    #Save Output
    out_dir = opts.output+SEP+SEPD
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    if len(opts.output) > 0: outpath = out_dir+'/pspec_boot%04d.npz' % boot
    else: outpath = 'pspec_boot%04d.npz' % boot
    print '   Writing '+outpath
    n.savez(outpath, kpl=kpl, scalar=scalar, times=n.array(lsts),
        pk_vs_t=pI, err_vs_t=1./cnt, temp_noise_var=var, nocov_vs_t=pI,
        afreqs=afreqs,chans=chans,cmd=' '.join(sys.argv))
