#! /usr/bin/env python
'''
Bootstrap antenna positions using only known source positions and measured visibilities.
1) Produce delay spectra for baselines and determine delays of peak flux (i.e. sources)
2) Using 3-antenna closures to co-identify peaks between baselines, so peaks are listed
   in a consistent order for each baseline.
3) Establish a (D = P A S) matrix equation, where:
    i) D is a (#bls,#dlys) matrix of measured delays in consistent source ordering
    ii) P is a (#bls,#ants) matrix pairing antennas for each baseline measurement
    iii) A is a (#ants,(x,y,z,tau)) matrix of coordinates for each antenna (including cable delay)
    iv) S is a ((x,y,z,1),#dlys) matrix of known source coordinates mapped to each delay bin.
4) Solve for A and permute source order in S to minimize z and tau in the antenna coordinates.
---------------------------------------------------------------------------------------
Currently this script works very well with simulated visibilities with several bright sources
up at one time.  However, it seems to fail pretty miserably on real data.  Things to improve:
1) Allow for only some sources to be up at a given time, and using the same source at different times.
2) Improve co-identification of delay spikes in real (noisy, many-source) data.
'''

import aipy as a, numpy as n, pylab as p
from aipy.miriad import bl2ij,ij2bl
import sys, optparse, itertools, random

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
aa.set_jultime(2455746.0)
cat = a.cal.get_catalog(opts.cal, srcs=['cen','for','pic','hyd','Sun'])
cat.compute(aa)
srcs = cat.keys()
freqs = aa.get_afreqs()
dlys_noshift = n.fft.fftfreq(uv['nchan'], uv['sdf'])
dlys = n.fft.fftshift(dlys_noshift)

#NDLY = 32
NDLY = 3*len(srcs)
NSRC = len(srcs)
NCHAN = 256
data = {}
closure = {}
ants = {}
times = []
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    NCHAN = uv['nchan']
    window = a.dsp.gen_window(uv['nchan'], window='blackman-harris')
    a.scripting.uv_selector(uv, ants=opts.ant)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or times[-1] != t:
            times.append(t)
        if len(times) > 10: continue
        ants[i] = ants[j] = None
        bl = ij2bl(i,j)
        if False: # Override inferred delays with computed delays
            taus = n.array([-aa.get_baseline(i,j,src=cat[s])[-1] for s in srcs] + [123, 456])
            taus = n.around(taus, -1) # Estimate resolution of delay bins
            #n.random.shuffle(taus) # Randomize source delay order
        else:
            _d = n.abs(n.fft.ifft(d * window))
            data[bl] = _d
            _d = n.fft.fftshift(_d)
            taubins = n.argsort(_d)
            taus = dlys[taubins[-NDLY:]][::-1]
            if False:
                true_taus = n.array([-aa.get_baseline(i,j,src=cat[s])[-1] for s in srcs])
                best_score,best_order = n.Inf, None
                for t in itertools.permutations(true_taus):
                    score = n.sqrt(n.average(n.abs(taus - t)**2))
                    if score < best_score: best_score,best_order = score, t
                if best_score > 1e3:
                    print (i,j), taus
                    print (i,j), n.around(n.array(best_order)), best_score
                    #p.plot(dlys[1:-1], _d_d2)
                    #p.plot(dlys, _d)
                    #p.show()
                    continue
        closure[bl] = taus
        #closure[bl] = n.array(best_order)

ants = ants.keys(); ants.sort()
NANT = len(ants)
bls = closure.keys()
N_BL = len(bls)

print 'Attempting to put delays in a consistent source order'
def get_bl(i,j):
    if i > j: return ij2bl(j,i)
    else: return ij2bl(i,j)
def do_conj(i,j,k, _dij, _djk):
    if i > j: _dij = n.conj(_dij)
    if j > k: _djk = n.conj(_djk)
    if i > k: _dij,_djk = n.conj(_dij),n.conj(_djk)
    return _dij, _djk
def gen_map(i,j,k, data, bins=None):
    '''Generate mapping of dly correspondence from ij to ik.'''
    blij,blik,bljk = get_bl(i,j), get_bl(i,k), get_bl(j,k)
    dij,dik,djk = data[blij],data[blik],data[bljk]
    _djk = n.fft.fft(djk)
    if bins is None: bins = n.arange(dij.size)
    else: bins = n.array(bins)
    ans = n.zeros((bins.size,dij.size), dtype=dij.dtype)
    for cnt, b in enumerate(bins):
        mask = n.zeros_like(dij); mask[b] = 1
        _dij = n.fft.fft(dij*mask)
        c_dij, c_djk = do_conj(i,j,k, _dij, _djk)
        ans[cnt] = n.abs(n.fft.ifft(c_dij * c_djk) * dik)**(1./3)
        #print b, bl2ij(blij), dlys_noshift[b], bl2ij(blik), n.sum(dlys_noshift * ans[b]) / n.sum(ans[b])
    return ans
def find_peaks(d):
    p = n.zeros(d.size)
    np = n.argsort(d)[-NDLY:]
    while n.any(p != np):
        p = np
        np = n.unique(n.where(d[p-1] > d[p], p-1, p))
        np = n.unique(n.where(d[(np+1)%d.size] > d[np], (np+1)%d.size, np))
    return np
    
    
bl01 = ij2bl(0,1)
bins = find_peaks(data[bl01])
ok_bls = {bl01: n.where(bins > NCHAN/2, bins-NCHAN, bins)}
NDLY = len(ok_bls[bl01])
if True:
    print bl2ij(bl01), bins, dlys_noshift[bins]
    p.plot(dlys_noshift, data[bl01])
    p.plot(dlys_noshift[bins], data[bl01][bins], '^')
    p.show()


#        inds = n.argsort(ans,axis=None)[-NDLY**2:]
#        best_inds = {}
#        for ind in inds:
#            x = ind / ans.shape[0]
#            y = ind % ans.shape[0]
#            _x = None
#            for cnt in range(100):
#                _y = n.argmax(ans[x])
#                _x = n.argmax(ans[:,_y])
#                if _x == x: break
#                else: x = _x
#            y = _y
#            ind = x * ans.shape[0] + y
#            best_inds[ind] = None
cnt, CNTMAX = {}, 10
prod = {}
done = {}
for k in range(2,NANT):
    for i,j in itertools.combinations(range(k),2):
        if k in (i,j): continue
        blij,blik,bljk = ij2bl(i,j),ij2bl(i,k),ij2bl(j,k)
        if not ok_bls.has_key(blij): continue
        cik = cnt.get(blik,[])
        cjk = cnt.get(bljk,[])
        if len(cik) < CNTMAX and not j in cik:
            prod[blik] = prod.get(blik,1) * gen_map(i,j,k,data,ok_bls[blij])
            cnt[blik] = cik + [j]
        if len(cjk) < CNTMAX and not i in cjk:
            prod[bljk] = prod.get(bljk,1) * gen_map(j,i,k,data,ok_bls[blij])
            cnt[bljk] = cjk + [i]
        if False:
            print i,j,k
            print ok_bls[bl01]
            p.imshow(n.log10(prod[blik]), interpolation='nearest', aspect='auto')
            p.colorbar()
            p.show()
    for bl in prod:
        if done.has_key(bl): continue
        bins = n.argmax(prod[bl], axis=-1)
        bins = n.where(bins > NCHAN/2, bins-NCHAN, bins)
        ok_bls[bl] = bins
        print bl2ij(bl), bins, len(cnt[bl])
        if len(cnt[bl]) >= CNTMAX: done[bl] = None
    # Redo internal baselines
    for _i,_k in itertools.combinations(range(k+1),2):
        blik = ij2bl(_i,_k)
        cik = cnt.get(blik,[])
        if len(cik) >= CNTMAX: continue
        for _j in range(k+1):
            if len(cik) >= CNTMAX: break
            if _j in [_i,_k]: continue
            prod[blik] = prod.get(blik,1) * gen_map(_i,_j,_k,data,ok_bls[ij2bl(_i,_j)])
            cik.append(_j)
        cnt[blik] = cik
        bins = n.argmax(prod[blik], axis=-1)
        bins = n.where(bins > NCHAN/2, bins-NCHAN, bins)
        if n.any(ok_bls[blik] != bins):
            print _i,_k
            print ok_bls[blik]
            print bins
        ok_bls[blik] = bins
    print k, '-'*30
        

#for i,j,k in itertools.combinations(range(NANT),3):
#    blij,blik,bljk = ij2bl(i,j),ij2bl(i,k),ij2bl(j,k)
#    if not ok_bls.has_key(blij): continue
#    if not ok_bls.has_key(blik):
#        ij2ik = gen_map(i,j,k,data, bins=ok_bls[blij])
#        bins = n.argmax(ij2ik, axis=-1)
#        bins = n.where(bins > NCHAN/2, bins-NCHAN, bins)
#        ok_bls[blik] = bins
#    if not ok_bls.has_key(bljk):
#        ij2jk = gen_map(j,i,k,data, bins=ok_bls[blij])
#        bins = n.argmax(ij2jk, axis=-1)
#        bins = n.where(bins > NCHAN/2, bins-NCHAN, bins)
#        ok_bls[bljk] = bins

#def unwrap(b): return n.where(b < 0, NCHAN+b,b)
#print (ok_bls[ij2bl(0,4)])
#p.subplot(131); p.imshow(n.log10(gen_map(0,1,4, data)), interpolation='nearest')
#p.subplot(132); p.imshow(n.log10(gen_map(0,2,4, data)), interpolation='nearest')
#p.subplot(133); p.imshow(n.log10(gen_map(0,3,4, data)), interpolation='nearest')
#for i,k in itertools.combinations(range(NANT),2):
#    blik = ij2bl(i,k)
#    prod = 1.
#    for j in range(NANT):
#        if j in [i,k]: continue
#        m = gen_map(i,j,k,data,ok_bls[ij2bl(i,j)])
#        prod *= m
#    bins = n.argmax(prod, axis=-1)
#    bins = n.where(bins > NCHAN/2, bins-NCHAN, bins)
#    print i,k
#    print ok_bls[blik]
#    print bins
#    ok_bls[blik] = bins
#p.show()


## Use more baselines to improve estimate of bin maps
## XXX I think that if one source is very strong, this process
## gradually reduces all bin entries to be to that source
## because of 2 srcs occasionally hitting the same delay bin on some baselines.
#redo_ants = {}
#for i in range(NANT): redo_ants[i] = None
#maps = {}
#for j in range(NANT):
#    new_redo = {}
#    new_ok_bls = {}
#    for i,k in itertools.combinations(range(NANT),2):
#        if not redo_ants.has_key(i) and not redo_ants.has_key(k): continue
#        jrange = range(cnt-NANT/4,cnt)
#        for j in jrange:
#            if j in [i,k]: continue
#            if i > j: blij = ij2bl(j,i)
#            else: blij = ij2bl(i,j)
#            blik = ij2bl(i,k)
#            maps[blik] = maps.get(blik,1) * gen_map(i,j,k,data, bins=ok_bls[blij])
#            bins = n.argmax(maps[blik], axis=-1)
#            bins = n.where(bins > NCHAN/2, bins-NCHAN, bins)
#            if n.average(n.abs(bins - ok_bls[blik])**2) > 1:
#                print i,j,k
#                print bins
#                print ok_bls[blik]
#                new_redo[i],new_redo[k] = None, None
#            new_ok_bls[blik] = bins
#    if len(new_redo) == 0: break
#    print cnt, new_redo.keys()
#    redo_ants = new_redo
#    ok_bls = new_ok_bls

## Cull repeat entries in ok_bls
#bins = {}
#print ok_bls[bl01]
#for i,b in enumerate(ok_bls[bl01]): bins[b] = i
#uniq = n.array(bins.values())
#print ok_bls[bl01][uniq]
#NDLY = len(uniq)
#print 'Culled to %d delay bins' % NDLY

#for bl in closure: closure[bl] = dlys_noshift[ok_bls[bl][uniq]]

for bl in closure: closure[bl] = dlys_noshift[ok_bls[bl]]

wgts = {}
print 'Checking closure'
for i,j,k in itertools.combinations(range(NANT),3):
    blij,blik,bljk = ij2bl(i,j),ij2bl(i,k),ij2bl(j,k)
    score = n.sqrt(n.average(n.abs(closure[blij] + closure[bljk] - closure[blik])**2))
    for bl in [blij,blik,bljk]:
        #wgts[bl] = min(wgts.get(bl,n.Inf), score)
        wgts[bl] = wgts.get(bl,0) + score
        if False:
            if wgts.has_key(bl): continue
            _i,_j = bl2ij(bl)
            actual = n.around(n.array([-aa.get_baseline(_i,_j,src=cat[s])[-1] for s in srcs]))
            #score = n.sqrt(n.average(n.abs(closure[bl][n.array([2,3,0,1])] - actual)))
            score = n.average(n.abs(closure[bl][n.array([2,3,0,1])] - actual))

print 'Done with source ordering on baseline level'

# Using well-ordered delays (modulo one overall permutation of sources), derive antenna positions
NCRD = 3
#NCRD = 2
P = n.zeros((N_BL+1,NANT), dtype=n.float)
V = n.zeros((N_BL+1,NDLY), dtype=n.float)

min_wgt = min([w for w in wgts.values() if w != 0])
bin_wgts = data[bl01][ok_bls[bl01]]
print bin_wgts
for nbl, bl in enumerate(bls):
    i,j = bl2ij(bl)
    #print i,j, wgts[bl] - min_wgt
    #if wgts[bl] != 0: continue
    #wgt = 1. / max(1, wgts[bl]**2)
    wgt = 1. / max(min_wgt**2/4, (wgts[bl]-min_wgt)**2)
    #wgt = 1.
    P[nbl,j], P[nbl,i] = wgt, -wgt
    V[nbl] = -n.array(closure[bl]) * wgt * bin_wgts
# Fix ant0 at origin
P[-1,0] = 1e6
V[-1,:] = 0

if True: # Write down what the answer is
    B = n.zeros((NANT,NCRD+1), dtype=n.float)
    for i in range(NANT):
        B[i,:NCRD] = aa.get_baseline(0,i)[:NCRD]
        B[i,NCRD] = aa[i]._phsoff[0] - aa[0]._phsoff[0]
    p.plot(B[:,0],B[:,1], 'ko')

print 'Deducing source order'

best_score, best_order = n.Inf, None
if True: # Override surmised source order with actual source order
    #best_order = ('cen','hyd','hyd','for','pic')
    best_order = ('cen','Sun','hyd','for','pic')
else:
    for src_order in itertools.product(*[srcs] * NDLY):
        S = n.zeros((NCRD+1,NDLY), dtype=n.float)
        for nsrc, src in enumerate(src_order):
            bwgt = bin_wgts[nsrc]
            S[:NCRD,nsrc] = cat[src].get_crds('top')[:NCRD] * bwgt
            S[NCRD,nsrc] = 1 * bwgt
        _S = n.linalg.pinv(S)
        B = n.linalg.lstsq(P, n.dot(V,_S))[0]
        score = n.sum(n.abs(B[:,2:])**2)
        if score < best_score: best_score, best_order = score, src_order
        print ' '.join(src_order), best_score, score
S = n.zeros((NCRD+1,NDLY), dtype=n.float)
for nsrc, src in enumerate(best_order):
    bwgt = bin_wgts[nsrc]
    S[:NCRD,nsrc] = cat[src].get_crds('top')[:NCRD] * bwgt
    S[NCRD,nsrc] = 1 * bwgt
_S = n.linalg.pinv(S)
    
print 'Running best permutation', best_order, best_score
print bin_wgts
for (i,j) in [(0,1), (0,2), (0,3), (0,4)]:
    bl = ij2bl(i,j)
    print (i,j)
    try:
        print closure[bl]
        print n.array([-aa.get_baseline(i,j,src=cat[s])[-1] for s in best_order])
    except(KeyError): pass
M = n.dot(V,_S)
Best = n.linalg.lstsq(P, M)[0]
p.plot(Best[:,0],Best[:,1], 'r.')

if False: # Plot 0,x baselines, which are proxy for ant pos
    for bl,crd in zip(bls,M):
        i,j = bl2ij(bl)
        print (i,j), n.around(crd)
        if i == 0: p.plot(crd[:1],crd[1:2], 'kx')

p.show()
