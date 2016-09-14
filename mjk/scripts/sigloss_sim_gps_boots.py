import matplotlib
matplotlib.use('Agg')
import numpy as np, pylab as plt
import capo, aipy
import capo.oqe as oqe, capo.frf_conv as fringe
import sys
from time import time
import sigloss_functions as sf
from progressbar import ProgressBar, ETA, Bar, Percentage
import ipdb
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

def errorbars(data,axis=1):
    mean = np.percentile(data,50,axis=axis)
    lower = mean - np.percentile(data, 2.5,axis=axis)
    upper = np.percentile(data,97.5,axis=axis) - mean
    return lower, upper

CH0,NCHAN = 95, 21
NSAMP = 609
NFOLDS = 2
NRUN = 100
SEP = '0,1'
POL = 'I'
inttime = 42.9499
normalize_mode='L^-1'
v_scale=np.logspace(-5,10,40)
#v_scale=[1e-3]
NSCALE=len(v_scale)
NBOOT=100
NGPS=5
PLOT=True

days=['even','odd']
toc = time()

chans = (np.arange(NCHAN) - int(NCHAN/2.)) + CH0
freqs = np.linspace(0.1,0.2,num=203)
afreqs = freqs.take(chans)
fq = np.average(afreqs)

aa = aipy.cal.get_aa('psa6240_v003', afreqs)
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
    if False: bl = 11072; ij =  aipy.miriad.bl2ij(bl);
    print 'Using Baseline for FRP:',bl
    bins = fringe.gen_frbins(inttime)
    frp, bins = fringe.aa_to_fr_profile(aa, ij, len(afreqs)/2, bins=bins)

    timebins, firs = fringe.frp_to_firs(frp, bins, aa.get_freqs(), fq0=aa.get_freqs()[len(afreqs)/2])#, fr_width_scale=1.3, maxfr=1.3e-3)

    if blconj[aipy.miriad.ij2bl(ij[0],ij[1])]: fir = {(ij[0],ij[1],POL):np.conj(firs)} #conjugate fir if needed
    else: fir = {(ij[0],ij[1],POL):firs}

qs_e,qs_v,qs_r,qs_ev = [], [], [], []
ps_e,ps_v,ps_r,ps_ev = [], [], [], []
c_nums = []
bar =ProgressBar(maxval=len(v_scale)*NRUN,widgets=['Performing MC:',Bar(),Percentage(),' ',ETA()]).start()
for cnt,sc in enumerate(v_scale):
    tmp_qs_e,tmp_qs_v,tmp_qs_r,tmp_qs_ev = [], [], [], []
    tmp_ps_e,tmp_ps_v,tmp_ps_r,tmp_ps_ev = [], [], [], []
    tmp_c=[]
    for run in xrange(NRUN):
        e,v,r = {} , {}, {}
        #Gen 1 eor for all bls
        v_ = oqe.noise(size=(NCHAN,NSAMP*NFOLDS)) * sc
        wij_ = np.zeros_like(v_.T)
        dij,wij = v_.T , np.logical_not(wij_)
        v_full,_w,_,_ = fringe.apply_frf(aa,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)

        for num,k in enumerate(days):
            e[k],v[k],r[k] = {}, {} ,{}
            for n1,bl in enumerate(all_bls):
                e[k][bl],v[k][bl],r[k][bl] = {}, {} ,{}


                e_ = oqe.noise(size=(NCHAN,NSAMP*NFOLDS)).T

                wij_ = np.zeros_like(e_)
                dij,wij = e_ , np.logical_not(wij_)
                e_,_w,_,_ = fringe.apply_frf(aa,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)

                r_ = np.copy(e_ + v_full)

                e_ = clip_array(e_.T,NSAMP,axis=1)
                v_1 = clip_array(v_full.T,NSAMP,axis=1)
                r_ = clip_array(r_.T,NSAMP,axis=1)

                r[k][bl][POL] = r_.T
                e[k][bl][POL] = e_.T
                v[k][bl][POL] = v_1.T
           #wij_ = np.zeros_like(r)
           #dij,wij = r , np.logical_not(wij_)
           #r,_w,_,_ = fringe.apply_frf(aa,dij,wij,ij[0],ij[1],pol=POL,bins=bins,firs=fir)

           #k = ('even',(0,1),'I')
           #k1 = ('e',(0,1),'I')
           #k2 = ('v',(0,1),'I')
        ds = oqe.DataSet(dsets=r)
        #print i, np.log10(np.linalg.cond(ds.C(k)))
        all_gps = ds.gen_bl_boots(NBOOT,ngps=NGPS)
        _qs_e, _qs_v, _qs_r =[],[],[]
        _ps_e, _ps_v, _ps_r =[],[],[]
        for nboot,gps in enumerate(all_gps):
            _qs_e.append([]); _qs_v.append([]); _qs_r.append([])
            _ps_e.append([]); _ps_v.append([]); _ps_r.append([])
            _tmp_c = []
            bls = [bl for gp in gps for bl in gp]
            _Cez, _Cvz, _Crz = {}, {},{}
            _Ce, _Cv, _Cr = {}, {},{}
            for k in days:
                _Cez[k], _Cvz[k], _Crz[k] = {}, {},{}
                _Ce[k], _Cv[k], _Cr[k] = {}, {},{}
                for i,gp in enumerate(gps):
                    _Crz[k][i],_Cez[k][i],_Cvz[k][i]= {},{},{}
                    _Cr[k][i],_Ce[k][i],_Cv[k][i]= {},{},{}
                    ds.set_data(r)
                    _Crz[k][i][POL] = sum( [np.dot(ds.iC((k,bl,POL)), ds.x[(k,bl,POL)]) for bl in gp]).T
                    _Cr[k][i][POL] = sum( [ds.iC((k,bl,POL)) for bl in gp])

                    ds.set_data(e)
                    _Cez[k][i][POL] = sum( [ds.x[(k,bl,POL)] for bl in gp]).T
                    _Ce[k][i][POL] = sum( [ np.identity(ds.iC((k,bl,POL)).shape[0]) for bl in gp])

                    ds.set_data(v)
                    _Cvz[k][i][POL] = sum( [ds.x[(k,bl,POL)] for bl in gp]).T
                    _Cv[k][i][POL] = sum( [ np.identity(ds.iC((k,bl,POL)).shape[0]) for bl in gp])
            for k_index,k1 in enumerate(days):
                for k2 in days[k_index:]:
                    for gp1 in _Crz[k1]:
                        for  gp2 in _Crz[k2]:
                            if gp1 == gp2 or k1 == k2: continue
                            key1 = (k1,gp1,POL)
                            key2 = (k2,gp2,POL)

                            ds.set_data(_Crz)
                            _tmp_c.append(np.log10(np.linalg.cond(ds.C(key1))))
                            ds.set_iC({key1:_Cr[k1][gp1][POL],key2:_Cr[k2][gp2][POL]})
                            _iC1x, _iC2x = np.fft.fft(_Crz[k1][gp1][POL].T.conj(), axis=0), np.fft.fft(_Crz[k2][gp2][POL].T.conj(), axis=0)
                            q_r =np.fft.fftshift(_iC1x,axes=0).conj() * np.fft.fftshift(_iC2x,axes=0)

                            _qs_r[nboot].append(q_r)
                            Fr = ds.get_F(key1,key2)
                            (Mr,Wr) = ds.get_MW(Fr,mode=normalize_mode)
                            p_r = ds.p_hat(Mr,q_r)
                            _ps_r[nboot].append(p_r)

                            ds.set_data(_Cez)
                            ds.set_iC({key1:_Ce[k1][gp1][POL],key2:_Ce[k2][gp2][POL]})
                            q_e = ds.q_hat(key1,key2,use_cov=False)
                            _qs_e[nboot].append(q_e)
                            Fe = ds.get_F(key1,key2)
                            (Me,We) = ds.get_MW(Fe,mode=normalize_mode)
                            p_e = ds.p_hat(Me,q_e)
                            _ps_e[nboot].append(p_e)

                            ds.set_data(_Cvz)
                            ds.set_iC({key1:_Cv[k1][gp1][POL],key2:_Cv[k2][gp2][POL]})
                            q_v = ds.q_hat(key1,key2,use_cov=False);
                            _qs_v[nboot].append(q_v)
                            Fv = ds.get_F(key1,key2)
                            (Mv,Wv) = ds.get_MW(Fv,mode=normalize_mode)
                            p_v = ds.p_hat(Mv,q_v)
                            _ps_v[nboot].append(p_v)
        #tmp_c.append(np.log10(np.linalg.cond(ds.C(k))))
        #iC_r = ds.iC(k)
        #ds.set_data({k:r.T}); q_r = ds.q_hat(k,k); tmp_qs_r.append(q_r)
        #F = ds.get_F(k,k)
        #(M,W) = ds.get_MW(F,mode=normalize_mode)
        #p_r = ds.p_hat(M,q_r); tmp_ps_r.append(p_r)

        #ds.set_data({k1:e.T});
        #iC_e = ds.iC(k1);
        #I_e = np.identity(iC_e.shape[0]); 
        #ds.set_iC({k1:I_e})
        #q_e = ds.q_hat(k1,k1); tmp_qs_e.append(q_e)
        #F = ds.get_F(k1,k1)
        #(M,W) = ds.get_MW(F,mode=normalize_mode)
        #p_e = ds.p_hat(M,q_e); tmp_ps_e.append(p_e)

        #ds.set_data({k2:v.T});
        #iC_v = ds.iC(k2); 
        #I_v = np.identity(iC_v.shape[0]); 
        #ds.set_iC({k2:I_v})
        #q_v = ds.q_hat(k2,k2); tmp_qs_v.append(q_v)
        #F = ds.get_F(k2,k2)
        #(M,W) = ds.get_MW(F,mode=normalize_mode)
        #p_v = ds.p_hat(M,q_v); tmp_ps_v.append(p_v)

        ##ds.set_iC({k:iC_r})
        #ds.set_data({k1:e.T, k2:v.T})
        #ds.set_iC({k1:iC_e, k2:iC_v})
        #q_ev = ds.q_hat(k1,k2); tmp_qs_ev.append(q_ev)
        #F = ds.get_F(k1,k2)
        #(M,W) = ds.get_MW(F,mode=normalize_mode)
        #p_ev = ds.p_hat(M,q_ev); tmp_ps_ev.append(p_ev)
        bar.update(cnt*NRUN +run +1)
        tmp_c.append(np.mean(_tmp_c))
        tmp_qs_e.append(np.sum(_qs_e,axis=1))
        tmp_qs_v.append(np.sum(_qs_v,axis=1))
        tmp_qs_r.append(np.sum(_qs_r,axis=1))
        tmp_ps_e.append(np.sum(_ps_e,axis=1))
        tmp_ps_v.append(np.sum(_ps_v,axis=1))
        tmp_ps_r.append(np.sum(_ps_r,axis=1))
    c_nums.append(tmp_c)
    qs_e.append(tmp_qs_e)
    qs_v.append(tmp_qs_v)
    qs_r.append(tmp_qs_r)
    qs_ev.append(tmp_qs_ev)
    ps_e.append(tmp_ps_e)
    ps_v.append(tmp_ps_v)
    ps_r.append(tmp_ps_r)
    #ps_ev.append(tmp_ps_ev)

qs_e, qs_v, qs_r = np.array(qs_e), np.array(qs_v), np.array(qs_r)
#, qs_ev , np.array(qs_ev)
ps_e, ps_v, ps_r = np.array(ps_e), np.array(ps_v), np.array(ps_r)
np.savez('sigloss_calc_nsamp_{0}_ch_{1}_ngps_{2}_nboot_{3}.npz'.format(NSAMP,CH0,NGPS,NBOOT),p_in=ps_v,p_noise=ps_e,p_outs=ps_r)
#, ps_ev , np.array(ps_ev)
qs_e, qs_v, qs_r = np.array(qs_e).reshape( NSCALE, NRUN ,NSAMP*NCHAN*NBOOT), np.array(qs_v).reshape( NSCALE , NRUN , NSAMP*NCHAN*NBOOT ), np.array(qs_r).reshape( NSCALE, NRUN, NSAMP*NCHAN*NBOOT)
##,qs_ev , np.array(qs_ev).reshape( qs_ev.shape[0],qs_ev.shape[1],NSAMP*NCHAN*NBOOT)
ps_e, ps_v, ps_r = np.array(ps_e).reshape( NSCALE, NRUN,NSAMP*NCHAN*NBOOT), np.array(ps_v).reshape(NSCALE, NRUN,NSAMP*NCHAN*NBOOT), np.array(ps_r).reshape( NSCALE, NRUN, NSAMP*NCHAN*NBOOT)
#,ps_ev , np.array(ps_ev).reshape( ps_ev.shape[0],ps_ev.shape[1],NSAMP*NCHAN)

c_nums = np.array(c_nums)
tic = time()
print '\nBoth in Log 10 Units'
print 'EoR Scale', 'Condition_number'
for sc,num in zip(v_scale,c_nums.mean(axis=1)):
    print '{0:.2f}'.format(np.log10(sc)), '\t\t{0:.2f}'.format(num)
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

print "\nTime taken: {0:.3f}min".format((tic - toc)/60.)

p_ins = np.median(np.abs(ps_v),axis=-1)
p_noises = np.median(np.abs(ps_e), axis=-1)
p_outs= np.median(np.abs(ps_r),axis=-1)# - np.median(np.abs(ps_e),axis=-1) 

ks = np.arange(len(v_scale)) - int(len(v_scale)/2.)

p_in = np.median(p_ins,axis=1)
p_in_err = errorbars(p_ins,axis=1)
p_out = np.median(p_outs,axis=1)
#xis = axisp_out_err = p_outs.std(axis=1)
p_out_err = errorbars(abs(p_outs),axis=1)
p_noise  = np.median(p_noises,axis=1)
p_noise_err = errorbars(p_noises, axis=1)

plt.figure()
plt.errorbar(p_in,p_out,xerr=p_in_err,yerr=p_out_err,fmt='b.')
plt.plot(p_in,p_in,'k--')
plt.errorbar(p_in,p_noise,p_noise_err,fmt='g-.')
plt.yscale('log')
plt.xscale('log')
plt.xlabel('input eor signal')
plt.ylabel('output power spectrum')
plt.grid()
plt.savefig('sigloss_errorbar_nsamp_{0}_ch_{1}_ngps_{2}_nboot_{3}.png'.format(NSAMP,CH0,NGPS,NBOOT),format='png')

plt.figure()
#hist,xbins,ybins = np.histogram2d(np.log10(p_ins.flatten()),np.log10(p_outs.flatten()),bins=(30,100),range=[[-20,20],[-20,20]])
H,xbins,ybins = np.histogram2d(np.log10(np.abs(ps_v).flatten()),np.log10(np.abs(ps_r).flatten()),bins=(500,500),range=[[-20,20],[-20,20]])
xb,yb = (xbins[1:]+xbins[:-1])/2.,(ybins[1:]+ybins[:-1])/2.

H_noise,xbins_n,ybins_n = np.histogram2d(np.log10(np.abs(ps_v).flatten()),np.log10(np.abs(ps_e).flatten()),bins=(500,500),range=[[-20,20],[-20,20]])
xb_n,yb_n = (xbins_n[1:]+xbins_n[:-1])/2.,(ybins_n[1:]+ybins_n[:-1])/2.

#H[H==0]=np.nan
H = H.T

#H_noise[H_noise==0]=np.nan
H_noise= H_noise.T
#exten = [xbins[0],xbins[-1],ybins[0],ybins[-1]]
#plt.imshow(hist,extent=exten,aspect='auto',interpolation='none')
plt.pcolor(10**(xb), 10**(yb), H)
plt.contour(10**(xb_n),10**(yb_n),H_noise,alpha=.8 )
plt.plot(p_in,p_in,'k--')
plt.plot(p_in,p_noise,'g-.')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('input eor signal')
plt.ylabel('output power spectrum')
plt.grid()
plt.savefig('sigloss_distribution_nsamp_{0}_ch_{1}_ngps_{2}_nboot_{3}.png'.format(NSAMP,CH0,NGPS,NBOOT),format='png')

p_eor, lower, upper = sf.get_pspec_from_sigloss(H,10**(xb),10**(yb), 3e-6,2e-5)
gn = sf.gauss(10**(yb),3e-6,2e-5)
gn.shape+=(1,)
plt.pcolor(10**(xb), 10**(yb), H*gn)
plt.savefig('sigloss_distribution_marginal_nsamp_{0}_ch_{1}_ngps_{2}_nboot_{3}.png'.format(NSAMP,CH0,NGPS,NBOOT),format='png')

plt.figure()
plt.plot(10**(xb), (H*gn).sum(0)/(H*gn).sum(0).sum())
plt.savefig('Peor_marginal_distribution_nsamp_{0}_ch_{1}_ngps_{2}_nboot_{3}.png'.format(NSAMP,CH0,NGPS,NBOOT),format='png')
print 'P_eor, lower, upper:', p_eor, lower, upper
if PLOT: plt.show()
