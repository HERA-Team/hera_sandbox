#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as P
import optparse
import sys, scipy

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

ij2bl,bl2ij = a.miriad.ij2bl,a.miriad.bl2ij

USE_CAL = True

A_ = [0,16,8,24,4,20,12,28]
B_ = [i+1 for i in A_]
C_ = [i+2 for i in A_]
D_ = [i+3 for i in A_]
ANTPOS = n.array([A_, B_, C_, D_])

bls = {}
conj = {}
for ri in range(ANTPOS.shape[0]):
    for ci in range(ANTPOS.shape[1]):
        for rj in range(ri,ANTPOS.shape[0]):
            for cj in range(ci,ANTPOS.shape[1]):
                sep = '%d,%d' % (rj-ri, cj-ci)
                bls[sep] = bls.get(sep,[]) + [(ANTPOS[ri,ci],ANTPOS[rj,cj])]
for sep in bls.keys():
    if sep == '0,0' or len(bls[sep]) < 2: del(bls[sep])
for sep in bls:
    conj[sep] = [i>j for i,j in bls[sep]]

strbls = {}
conj_bl = {}
for sep in bls:
    strbls[sep] = []
    bl_list = []
    for (i,j),c in zip(bls[sep],conj[sep]):
        if c: i,j = j,i
        bl_list.append(ij2bl(i,j))
        strbls[sep].append('%d_%d' % (i,j))
        conj_bl[ij2bl(i,j)] = c
    bls[sep] = bl_list
    strbls[sep] = ','.join(strbls[sep])
    #print sep, strbls[sep]

uv = a.miriad.UV(args[0])
fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
(uvw,t,(i,j)),d = uv.read()
del(uv)
if USE_CAL: 
    print 'Applying calibration information'
    aa = a.cal.get_aa(opts.cal, fqs)
    aa.set_jultime(t)

seps = ['0,1']
plot_bls = bls[seps[0]]
strbls = ','.join([strbls[sep] for sep in seps])
#print strbls

times, d, f = C.arp.get_dict_of_uv_data(args, strbls, 'xx', verbose=True)
for bl in d:
    i,j = bl2ij(bl)
    if USE_CAL:
        d[bl] = aa.phs2src(d[bl], 'z', i,j)
        d[bl] /= n.median(aa.passband(i,j))
        #if i == 8 or j == 8: d[bl] = -d[bl]  # XXX remove this line once correct script is run
    if conj_bl.has_key(bl) and conj_bl[bl]: d[bl] = n.conj(d[bl])

w = {}
for bl in f: w[bl] = n.logical_not(f[bl]).astype(n.float)

NANT = ANTPOS.size
ANTIND = dict(zip(ANTPOS.flatten(), n.arange(ANTPOS.size)))

calrow,calcol = 0,0
calant = ANTPOS[calrow,calcol]
cal_bl = {}
for sep in seps: cal_bl[sep] = [bl for bl in bls[sep] if bl2ij(bl)[0] == ANTPOS[calrow,calcol]][0]

#flags = n.load('../zen.2455907.00656.uvGR.npz')
flags = n.load(args[0]+'R.npz')
f = n.zeros_like(flags['corrflags'])
for thing in flags.files:
    f |= flags[thing]

dlys = []
amp = []
phs = []
P.figure(0)
for sep in seps:
    cbl = cal_bl[sep]
    d[cbl] = n.ma.array(d[cbl],mask=f)
    i0,j0 = bl2ij(cbl)
    if conj_bl.has_key(cbl) and conj_bl[cbl]: i0,j0 = j0,i0
    for bl in bls[sep]:
        d[bl] = n.ma.array(d[bl],mask=f)
        if bl == cbl or not bl in plot_bls: continue
        i,j = bl2ij(bl)
        if conj_bl.has_key(bl) and conj_bl[bl]: i,j = j,i
        if USE_CAL: 
            dlys.append(aa[j].get_params(['dly'])['dly'] - aa[i].get_params(['dly'])['dly'])
            print dlys[-1]
        g = n.ma.mean(d[bl]/d[cbl],axis=0) 
        amp.append(n.abs(g))
        phs.append(n.ma.arctan2(g.real,g.imag))
        P.subplot(211); P.semilogy(fqs,n.abs(g), label='%d,%d'%(i,j))
        P.subplot(212); P.plot(fqs,n.ma.arctan2(g.real,g.imag), label='%d,%d'%(i,j))
 
P.legend()
P.draw()
print 'FINAL DELAYS'
print dlys - n.mean(dlys)

P.figure(1)
m_amp = n.where(amp == 0, n.NaN,n.ma.mean(amp,axis=0)) 
s_amp = n.where(amp == 0, n.NaN, n.ma.std(amp,axis=0)) 
m_phs = n.where(amp == 0, n.NaN,n.ma.mean(phs,axis=0)) 
s_phs = n.where(amp == 0, n.NaN, n.ma.std(phs,axis=0)) 
P.subplot(211)
P.plot(fqs,m_amp,'b')
P.plot(fqs,m_amp+s_amp,'g')
P.plot(fqs,m_amp-s_amp,'g')
P.subplot(212)
P.plot(fqs,m_phs,'b')
P.plot(fqs,m_phs+s_phs,'g')
P.plot(fqs,m_phs-s_phs,'g')

P.draw()

P.show()


    #if True:
    #    print bl2ij(bl), 'Delay:', tau
    #    #P.subplot(211); P.plot(fqs, n.abs(g))
    #    #P.subplot(212); P.plot(fqs, n.angle(g))
    #    #P.plot(fqs, n.angle(n.sum(d[cal_bl]*n.conj(d[bl]), axis=0)))
    #else:
    #    taus = []
    #    tau = 0.
    #    for t in range(0,d[cal_bl].shape[0],S):
    #        g,tau = C.arp.redundant_bl_cal(d[cal_bl][t:t+S],w[cal_bl][t:t+S],d[bl][t:t+S],w[bl][t:t+S],fqs,use_offset=False, tau=tau)
    #        print 'Delay:', tau
    #        taus.append(tau)
    #    P.plot(taus)
    
#P.show()
