import numpy as np, pylab as plt
import capo, aipy
import capo.oqe as oqe, capo.frf_conv as fringe
import sys

def clip_array(d,size,axis=0):
    if size is None:
        print 'Must give valid size for clipping'
        print 'Returning origian larray'
        return d
    d_size = d.shape[axis]
    if d_size < size:
        print 'Array Must have large size than desired clip'
        print 'Clipping failed, returning original array'
        return d
    start= int( (d_size - size)/2.)
    end = int( (d_size + size)/2. )
    diff = size - ( end - start )
    if diff != 0: end += diff
    d = np.swapaxes(d,axis,0)
    _d = d[start:end]
    _d = np.swapaxes(_d,0, axis)
    return _d

CH0,NCHAN = 30, 21
NSAMP = 609
NTIMES = 3
NRUN = 10
SEP = '0,1'
POL = 'I'
inttime = 42.9499

chans = (np.arange(NCHAN) - int(NCHAN/2.)) + CH0
freqs = np.linspace(0.1,0.2,num=203)
afreqs = freqs.take(chans)
fq = np.average(afreqs)

aa = aipy.cal.get_aa('psa6240_v003', freqs)
sep2ij, blconj, bl2sep = capo.zsa.grid2ij(aa.ant_layout)
ijs = sep2ij[SEP].split(',')
all_bls= [ aipy.miriad.ij2bl(*map( int,x.split('_'))) for x in ijs]

if True: #this one is the exact one
    sep = bl2sep[all_bls[0]]
    ij_array =  sep2ij[sep].split(',')
    while True:
        ij = map( int, ij_array.pop().split('_') )
        bl = aipy.miriad.ij2bl(*ij )
        if not blconj[bl]: break
    print 'Using Baseline for FRP:',bl
    bins = fringe.gen_frbins(inttime)
    frp, bins = fringe.aa_to_fr_profile(aa, ij, len(afreqs)/2, bins=bins)

    timebins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(), fq0=aa.get_freqs()[len(afreqs)/2])

    if blconj[aipy.miriad.ij2bl(ij[0],ij[1])]: fir = {(ij[0],ij[1],POL):np.conj(firs)} #conjugate fir if needed
    else: fir = {(ij[0],ij[1],POL):firs}

v_scale=np.logspace(-3,5,10)
qs_e,qs_v,qs_r,qs_ev = [], [], [], []
ps_e,ps_v,ps_r,ps_ev = [], [], [], []
for sc in v_scale:
    tmp_qs_e,tmp_qs_v,tmp_qs_r,tmp_qs_ev = [], [], [], []
    tmp_ps_e,tmp_ps_v,tmp_ps_r,tmp_ps_ev = [], [], [], []
    for i in xrange(NRUN):

           v = oqe.noise(size=(NCHAN,NSAMP*NTIMES)) * sc
           wij_ = np.zeros_like(v.T)
           dij,wij = v.T , np.logical_not(wij_)
           v,_w,_,_ = fringe.apply_frf(aa,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)

           e = oqe.noise(size=(NCHAN,NSAMP*NTIMES)).T
           r = np.copy(e + v)

           wij_ = np.zeros_like(r)
           dij,wij = r , np.logical_not(wij_)
           r,_w,_,_ = fringe.apply_frf(aa,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)

           wij_ = np.zeros_like(e)
           dij,wij = e , np.logical_not(wij_)
           e,_w,_,_ = fringe.apply_frf(aa,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)

           e = clip_array(e.T,NSAMP,axis=1)
           v = clip_array(v.T,NSAMP,axis=1)
           r = clip_array(r.T,NSAMP,axis=1)

           k = ('even',(0,1),'I')
           k1 = ('e',(0,1),'I')
           k2 = ('v',(0,1),'I')

           ds = oqe.DataSet(dsets={k:r.T})
           #print i, np.linalg.cond(ds.C(k))
           iC_r = ds.iC(k)

           ds.set_data({k:e.T}); q_e = ds.q_hat(k,k); tmp_qs_e.append(q_e)
           ds.set_data({k:v.T}); q_v = ds.q_hat(k,k); tmp_qs_v.append(q_v)
           ds.set_data({k:r.T}); q_r = ds.q_hat(k,k); tmp_qs_r.append(q_r)
           F = ds.get_F(k,k)
           ds.set_iC({k:iC_r})
           ds.set_data({k1:e.T, k2:v.T})
           ds.set_iC({k1:iC_r, k2:iC_r})
           q_ev = ds.q_hat(k1,k2); tmp_qs_ev.append(q_ev)
           (M,W) = ds.get_MW(F,mode='F^-1')
           p_e = ds.p_hat(M,q_e); tmp_ps_e.append(p_e)
           p_v = ds.p_hat(M,q_v); tmp_ps_v.append(p_v)
           p_r = ds.p_hat(M,q_r); tmp_ps_r.append(p_r)
           p_ev = ds.p_hat(M,q_ev); tmp_ps_ev.append(p_ev)

    qs_e.append(tmp_qs_e)
    qs_v.append(tmp_qs_v)
    qs_r.append(tmp_qs_r)
    qs_ev.append(tmp_qs_ev)
    ps_e.append(tmp_ps_e)
    ps_v.append(tmp_ps_v)
    ps_r.append(tmp_ps_r)
    ps_ev.append(tmp_ps_ev)

qs_e, qs_v, qs_r, qs_ev = np.array(qs_e), np.array(qs_v), np.array(qs_r), np.array(qs_ev)
ps_e, ps_v, ps_r, ps_ev = np.array(ps_e), np.array(ps_v), np.array(ps_r), np.array(ps_ev)
ps_e, ps_v, ps_r, ps_ev = np.array(ps_e).reshape( ps_e.shape[0],ps_e.shape[1],NSAMP*NCHAN), np.array(ps_v).reshape(ps_v.shape[0],ps_v.shape[1],NSAMP*NCHAN), np.array(ps_r).reshape( ps_r.shape[0],ps_r.shape[1],NSAMP*NCHAN), np.array(ps_ev).reshape( ps_ev.shape[0],ps_ev.shape[1],NSAMP*NCHAN)
qs_e, qs_v, qs_r, qs_ev = np.array(qs_e).reshape( qs_e.shape[0],qs_e.shape[1],NSAMP*NCHAN), np.array(qs_v).reshape(qs_v.shape[0],qs_v.shape[1],NSAMP*NCHAN), np.array(qs_r).reshape( qs_r.shape[0],qs_r.shape[1],NSAMP*NCHAN), np.array(qs_ev).reshape( qs_ev.shape[0],qs_ev.shape[1],NSAMP*NCHAN)

#print '\nq:'
#print '\tq_e:',np.mean(qs_e,axis=-1)
#print '\tq_v:',np.mean(qs_v,axis=-1)
#print '\tq_r:',np.mean(qs_r,axis=-1)
#print '\tq_ev:',np.mean(qs_ev,axis=-1)
#print 'p:'
#print '\tp_e:',np.mean(ps_e,axis=-1)
#print '\tp_v:',np.mean(ps_v,axis=-1)
#print '\tp_r:',np.mean(ps_r,axis=-1)
#print '\tp_ev:',np.mean(ps_ev,axis=-1)


p_ins = np.abs(ps_v).mean(axis=-1)
p_outs= np.abs(ps_r).mean(axis=-1) #- np.abs(ps_e).mean(axis=-1) 

ks = np.arange(len(v_scale)) - int(len(v_scale)/2.)

p_in = np.mean(p_ins,axis=1)
p_out = np.mean(p_outs,axis=1)
p_out_err = p_outs.std(axis=1)

plt.errorbar(p_in,p_out,p_out_err)
plt.plot(p_in,p_in,'k--')
plt.yscale('log')
plt.xscale('log')
plt.grid()
plt.show()
