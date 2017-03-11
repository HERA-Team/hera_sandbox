#! /usr/bin/env python

import numpy as np, aipy, capo, pylab as plt, sys, glob

def dB(sig): return 10*np.log10(np.abs(np.average(sig.real, axis=1)))

def find_sep(aa, bls, drow=None, dcol=None):
    layout = aa.ant_layout
    rv = []
    for i,j in bls:
        irow,icol = np.where(layout == i)
        jrow,jcol = np.where(layout == j)
        if not drow is None and abs(irow - jrow) != drow: continue
        if not dcol is None and abs(icol - jcol) != dcol: continue
        rv.append((i,j))
    return rv

SEPS = [
    (0,103) , #  1
    #(1,4) ,   #  2
    (0,26), # 2 for psa128 v3
    (0,101) , #  3
    (0,62) ,  #  4
]
'''
    (0,100) , #  5
    (1,13) ,  #  6
    (1,70) ,  #  7
    (1,56) ,  #  8
    (1,71) ,  #  9
    (1,59) ,  # 10
    (0,97) ,  # 11
    (12,43) , # 12
    (9,71) ,  # 13
    (9,59) ,  # 14
    (57,64) , # 15
]'''

CONJ = [ # XXX why conjugating all of these?
    (0,103) , #  1
    #(1,4) ,   #  2
    (0,101) , #  3
    (0,62) ,  #  4
    (0,100) , #  5
    (0,97) ,  # 11
    (12,43) , # 12
    (57,64) ] # 15

#SEPS = [(0,103), (0,111), (0,95)]
#SEPS += [(2,105), (1,83)]
#SEPS += [(0,79), (0,78)]
#SEPS += [(0,70),(0,71)]
#SEPS += [(1,107),(0,51)]
#SEPS += [(3,105),(3,106)]
#CH0,NCHAN = 0, 203
bandpass = np.load('bandpass.npz')['bandpass']; bandpass.shape = (1,-1)
fqs = np.linspace(.1,.2,bandpass.size)
aa = aipy.cal.get_aa('psa6622_v003', fqs)
#CH0,NCHAN = 110,51
CH0,NCHAN = 110,1
#CH0,NCHAN = 0,203
aa.select_chans(np.arange(CH0,CH0+NCHAN))
afreqs = aa.get_afreqs()
WINDOW = 'blackman-harris'
win = aipy.dsp.gen_window(NCHAN, WINDOW); win.shape = (-1,1)
fq = np.average(afreqs)
z = capo.pspec.f2z(fq)
sdf = fqs[1]-fqs[0]
B = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW['none'] #proper normalization
etas = np.fft.fftshift(capo.pspec.f2eta(afreqs)) #create etas (fourier dual to frequency)
kpl = etas * capo.pspec.dk_deta(z)
bm = np.polyval(capo.pspec.DEFAULT_BEAM_POLY, fq) * 2.35 # correction for beam^2
scalar = capo.pspec.X2Y(z) * bm * B
B_win = sdf * afreqs.size * capo.pfb.NOISE_EQUIV_BW[WINDOW] #proper normalization
scalar_win = capo.pspec.X2Y(z) * bm * B_win

sets = {
    #'day0' : sys.argv[1:],
    'day0' : glob.glob('zen.2456680.*.xx.npz'),
    'day1' : glob.glob('zen.2456681.*.xx.npz'),
    'day2' : glob.glob('zen.2456682.*.xx.npz'),
    #'day3' : glob.glob('zen.2456683.*.xx.npz'),
    'day4' : glob.glob('zen.2456684.*.xx.npz'),
    'day5' : glob.glob('zen.2456685.*.xx.npz'),
}
data,wgts = {}, {}
lsts = {}
chisqs = {}
for s in sets:
    if not lsts.has_key(s):
        meta, gains, vismdl, xtalk = capo.omni.from_npz(sets[s], bls=SEPS, verbose=True)
        lsts[s] = meta['lsts']
    chisqs[s] = meta['chisq'][:,CH0:CH0+NCHAN]
    for pol in vismdl:
        #for bl in vismdl[pol]:
        for bl in SEPS:
            k = (s,pol,bl)
            data[k] = vismdl[pol][bl][:,CH0:CH0+NCHAN]
            if bl in CONJ: data[k] = data[k].conj()
            data[k] *= bandpass[:,CH0:CH0+NCHAN]
            wgts[k] = np.where(np.abs(data[k]) == 0, 0., 1)
            #wgts[k] = np.where(np.abs(data[k]) == 0, 0., 1./chisqs[s])
set1,set2 = sets.keys()[0], sets.keys()[-1]
lst_res = np.average(lsts[set1][1:] - lsts[set1][:-1])
'''
data_g,wgts_g = {}, {}
for k in data:
    print 'Gridding', k
    lsts_g,data_g[k],wgts_g[k] = capo.oqe.lst_grid(lsts[k[0]], data[k], wgts=wgts[k])

def k2eq(k):
    return 'g%s * bl%d_%d' % ((k[0],) + k[-1])
data_eqs, wgts_eqs = {}, {}
sol0 = {}
for k in data_g:
    if not k[-1] == (0,103): continue
    data_eqs[k2eq(k)], wgts_eqs[k2eq(k)] = data_g[k], wgts_g[k]
    sol0['g'+k[0]] = np.ones_like(data_g[k])
    sol0['bl%d_%d' % k[-1]] = data_g[k]
ls = capo.linsolve.LinProductSolver(data_eqs, wgts_eqs, sol0)
sol1 = ls.solve()
for cnt,i in enumerate([1,2,4,5]):
    plt.subplot(1,4,cnt+1)
    capo.plot.waterfall(sol1['gday%d' % i], mode='lin', mx=1.2, drng=.4)
plt.show()
import IPython; IPython.embed()
'''

inds = capo.oqe.lst_align(lsts, lstres=lst_res)
data,wgts = capo.oqe.lst_align_data(inds, dsets=data, wgts=wgts)
for s in sets: chisqs[s] = chisqs[s][inds[s]].T
ds = capo.oqe.DataSet(data, wgts)
    
#k1a,k1b,k1c = [(s,'xx',(0,103)) for s in sets]
#k2a,k2b,k2c = [(s,'xx',(0,111)) for s in sets]
#k3a,k3b,k3c = [(s,'xx',(0, 95)) for s in sets]
#ks = [k1a,k1b,k1c,k2a,k2b,k2c,k3a,k3b,k3c]
#k1a,k1b = [(s,'xx',(0,103)) for s in sets]
#k2a,k2b = [(s,'xx',(0,111)) for s in sets]
#k3a,k3b = [(s,'xx',(0, 95)) for s in sets]
#ks = [k1a,k1b,k2a,k2b,k3a,k3b]
ks = [(s,'xx',bl) for bl in SEPS for s in sets]
#k1a, = [(s,'xx',(0,103)) for s in sets]
#k2a, = [(s,'xx',(0,111)) for s in sets]
#k3a, = [(s,'xx',(0, 95)) for s in sets]
#ks = [k1a,k2a,k3a]
NK = len(ks)

def set_C(norm=3e-6, tstart=400, tend=None, skip=None):
    ds.clear_cache()
    Cs,iCs = {},{}
    for k in ks:
        #sumbls = [ki for ki in ks if (ki != k) and (ki[2] != k[2])]
        #sumbls = [ki for ki in ks if (ki != k) and (ki[2] == k[2])]
        sumbls = [ki for ki in ks if (ki == k) and (ki[2] == k[2])]
        #sumbls = [ki for ki in ks if (ki[2] == k[2])]
        print k, sumbls
        Cs[k] = sum([capo.oqe.cov(ds.x[ki][:,tstart:tend:skip],ds.w[ki][:,tstart:tend:skip])+norm*np.identity(NCHAN) for ki in sumbls])
        ds.set_C({k:Cs[k]})
        iCs[k] = ds.iC(k)

def get_p(k1,k2,mode):
    assert(mode in 'IWC')
    if mode == 'I':
        qI = ds.q_hat(k1,k2,use_cov=False)
        FI = ds.get_F(k1,k2,use_cov=False)
        MI,WI = ds.get_MW(FI, mode='I')
        pI = ds.p_hat(MI,qI,scalar=scalar)
        return pI
    elif mode == 'W':
        pW = 1.6*2*np.fft.fftshift(np.fft.ifft(win*data[k1].T, axis=0) * np.fft.ifft(win*data[k2].T, axis=0).conj(), axes=0)
        return pW * scalar_win
    elif mode == 'C':
        qC = ds.q_hat(k1,k2)
        FC = ds.get_F(k1,k2)
        MC,WC = ds.get_MW(FC, mode='F^-1/2')
        pC = ds.p_hat(MC,qC,scalar=scalar)
        return pC

#set_C(1e-6)
set_C(4e4)
#pI,pW,pC = get_p(ks[0],ks[1])

for cnt,k in enumerate(ks):
    plt.subplot(NK,1,cnt+1)
    capo.plot.waterfall(ds.x[k], drng=3)
    plt.colorbar()
plt.show()

for cnt,k in enumerate(ks):
    #plt.subplot(3,2,cnt+1)
    #capo.plot.waterfall(ds.iC(k), mode='phs')
    #capo.plot.waterfall(ds.iC(k), mx=0, drng=7)
    plt.plot(ds.iC(k).diagonal())
    #plt.colorbar()
plt.show()

import IPython; IPython.embed()

for cnt,k in enumerate(ks):
    #plt.subplot(NK,1,cnt+1)
    pC = get_p(k,k,'C')
    #pC = get_p(k,('day2',)+k[1:],'C')
    plt.title(k[0])
    capo.plot.waterfall(pC, mx=16, drng=7)
    #capo.plot.waterfall(pC.real, mx=1e9, drng=2e9)
    #capo.plot.waterfall(pC.real)
    #plt.plot(np.abs(np.median(pC[:,400:], axis=1).real))
    #plt.plot(np.average(pC[:,400:], axis=1).real)
    plt.colorbar()
    plt.show()

#k1 = ('day1','xx',(0,103))
#k2 = ('day2','xx',(0,103))
k1 = ('day1','xx',(0,111))
k2 = ('day2','xx',(0,111))
plt.subplot(131)
capo.plot.waterfall(ds.iC(k1), drng=3); plt.colorbar()
plt.subplot(132)
capo.plot.waterfall(ds.iC(k2), drng=3); plt.colorbar()
plt.subplot(133)
capo.plot.waterfall(ds.iC(k1)-ds.iC(k2), drng=3); plt.colorbar()
plt.show()

'''
pC1 = get_p(k1a,k1b,'C')
pC2 = get_p(k2a,k2b,'C')
pC3 = get_p(k3a,k3b,'C')
plt.plot(np.abs(np.median(pC1, axis=1).real), 'r', label='pC1',)
plt.plot(np.abs(np.median(pC2, axis=1).real), 'b', label='pC2',)
plt.plot(np.abs(np.median(pC3, axis=1).real), 'k', label='pC3',)
plt.legend()
plt.show()
'''

import IPython; IPython.embed()
