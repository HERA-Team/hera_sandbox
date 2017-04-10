#! /usr/bin/env python
import aipy as a, numpy as np, capo as C, pylab as p
from joblib import Parallel, delayed
import sys
import multiprocessing
num_cores = multiprocessing.cpu_count()
#@p.ion()
#fqs = np.linspace(.1,.2,203)
fq = .15
bl1, bl2 = (0,103),(0,95)
N = 2400   #number of universes to average over

VIS = False
MAXCORR_EQUIV = 1228687.26108#1261990. #MAXCORR_EQUIV is the peak of baseline_ton

aa = a.cal.get_aa('psa6622_v001',np.array([fq])) #128
#aa = a.cal.get_aa('psa6240_v003', np.array([fq]))
#h = a.healpix.HealpixMap(nside=256)
h = a.healpix.HealpixMap(nside=64)
h.set_interpol(False)
ex,ey,ez = h.px2crd(np.arange(h.map.size), ncrd=3)

#IM = a.img.Img(size=200)
if VIS:
    p.ion()
    img = a.img.Img(200,.5)
    tx2,ty2,tz2 = img.get_top((200,200))
    SH = tx2.shape
    tx2,ty2,tz2 = tx2.flatten(),ty2.flatten(),tz2.flatten()
    top2 = np.array([tx2,ty2,tz2], dtype=tx2.dtype)
eq = np.array([ex,ey,ez], dtype=ex.dtype)

TT = np.arange(2455700.3,2455700.7,0.001)
NORM = float(TT.size)/1000.
def prepare(TT):
    aa.set_jultime(2455700.5)
    bl1x,bl1y,bl1z = aa.get_baseline(bl1[0],bl1[1],'z')
    bl2x,bl2y,bl2z = aa.get_baseline(bl2[0],bl2[1],'z')
    bm_fngs = None
    for j, jd in enumerate(TT):
        # convert tx,ty,tz to ex,ey,ez (time dependent)
        aa.set_jultime(jd)
        m = aa.eq2top_m
        tx,ty,tz = np.dot(m, eq)
        bl1_prj = tx*bl1x + ty*bl1y + tz*bl1z
        bl2_prj = tx*bl2x + ty*bl2y + tz*bl2z
        fng1 = np.exp(-2j*np.pi*bl1_prj*fq)
        fng2 = np.exp(-2j*np.pi*bl2_prj*fq)
        bm = aa[0].bm_response((tx,ty,tz),pol='I')[0]**2#/np.abs(tz)  #baseline beam is antenna beam squared.
        bm = np.where(tz > 0.001, bm, 0)   
        # bm /= bm.sum()   #avoiding the monopole
        #print bm.shape
        if bm_fngs is None:
            bm_fngs = np.zeros((2,TT.size,bm.shape[0]),dtype=np.dtype('c8'))
        bm_fng1 = bm * fng1
        bm_fng2 = bm * fng2
        bm_fngs[:,j,:] = [bm_fng1, bm_fng2]
    return bm_fngs

def find_corr(i, bm_fngs):
    print i
    sys.stdout.flush()
    sky = np.random.normal(size=h.map.size)
    h.map = sky # assume sky is in eq coord
    #import IPython; IPythonp.embed()
    sky_prj = h.map
    #sky_prj = np.ones_like(h.map)
    ###################
    if VIS:
        ex2,ey2,ez2 = np.dot(m, top2)
        sky_prj2 = h[ex2,ey2,ez2]
        sky_prj2.shape = SH
        if plt is None: plt = p.imshow(sky_prj2, vmax=2, vmin=-2,origin='lower')
        else:
           plt.set_data(sky_prj2)
           p.draw()
    ####################
    bm_fng1, bm_fng2 = bm_fngs[0], bm_fngs[1]
    vis1 = np.sum((sky_prj * bm_fng1),axis=1)
    vis2 = np.sum((sky_prj * bm_fng2),axis=1)

    _vis1,_vis2 = np.fft.fft(vis1), np.fft.fft(vis2)
    #print _vis1.shape
    #each universe has one corr
    temp = np.fft.ifftshift(np.fft.ifft(_vis2*np.conj(_vis1)))
    #p.plot(np.abs(corr[i]))
    #print temp.shape
    return temp

if __name__ == '__main__':
    print 'preparing bfs'
    bfs = prepare(TT)
    print 'done'

    corr = Parallel(n_jobs=8)(delayed(find_corr)(i, bfs) for i in xrange(N))
    corr = np.array(corr)
    print 'shape of corr:',corr.shape

    corr = corr/NORM
    #cuedict = {'26_23':0.125, '26_38': 0.042, '26_50': 0.083,'26_26':0., '50_57':0.122}
    cuedict = {'26_26':0.,'26_38': 0.032557,'26_50':0.073557,'13_32':0.030557,'13_14':0.066557,'50_59':0.071557}
    try: ver = cuedict[str(bl1[1])+'_'+str(bl2[1])]
    except(KeyError): ver = 0.
    meancorr = np.mean(corr,axis=0)
    meancorr = meancorr/1.e6        #MAXCORR_EQUIV is the peak of baseline_ton
    maxind = np.argmax(np.abs(meancorr))
    absmax = np.abs(meancorr[maxind])
    print '############## baseline_toff RESULT for', bl1, bl2, '#####################'
    print "Peak: ", meancorr[maxind], 'abs=',absmax, 'at dT = ', TT[np.argmax(np.abs(meancorr))]-2455700.5
    #meancorr = meancorr/absmax
    #import IPython; IPythonp.embed()
    blstr = str(bl1[0])+'_'+str(bl1[1])+'_'+str(bl2[0])+'_'+str(bl2[1])
    np.savez('blout_'+blstr, TT-2455700.5,meancorr)

    p.figure()
    p.subplot(211)
    p.plot(TT-2455700.5,np.real(meancorr))
    p.plot(TT-2455700.5,np.imag(meancorr))
    p.axvline(ver,color='k',alpha=0.5,linewidth=3)
    p.grid()
    p.subplot(212)
    p.plot(TT-2455700.5,np.abs(meancorr))
    p.axvline(ver,color='k',alpha=0.5,linewidth=3)
    p.grid()
    p.savefig('toff_24000')
#p.show()


