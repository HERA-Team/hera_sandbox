#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p
from joblib import Parallel, delayed
import sys
import multiprocessing
from baseline_toff_para import find_corr, prepare
num_cores = multiprocessing.cpu_count()

#@p.ion()
#fqs = n.linspace(.1,.2,203)
fq = .15
bl1, bl2 = (0,103),(0,103)
N = 240   #number of universes to average over

VIS = False
REDNORM = 1.
aa = a.cal.get_aa('psa6622_v001',n.array([fq])) #128
#aa = a.cal.get_aa('psa6240_v003', n.array([fq]))
#h = a.healpix.HealpixMap(nside=256)
h = a.healpix.HealpixMap(nside=64)
h.set_interpol(False)
ex,ey,ez = h.px2crd(n.arange(h.map.size), ncrd=3)

#IM = a.img.Img(size=200)
if VIS:
    p.ion()
    img = a.img.Img(200,.5)
    tx2,ty2,tz2 = img.get_top((200,200))
    SH = tx2.shape
    tx2,ty2,tz2 = tx2.flatten(),ty2.flatten(),tz2.flatten()
    top2 = n.array([tx2,ty2,tz2], dtype=tx2.dtype)
eq = n.array([ex,ey,ez], dtype=ex.dtype)

plt = None
#TT = n.arange(2455700,2455701,1/a.const.sidereal_day*42.9*0.5) #*5 for speed
TT = n.arange(2455700.3,2455700.7,0.001)
NORM = float(TT.size)/1000.
#for i in xrange(N):


print 'preparing bfs'
bfs = prepare(TT)
print 'done'

corr = Parallel(n_jobs=12)(delayed(find_corr)(i, bfs) for i in xrange(N))
corr = n.array(corr)
print 'shape of corr:',corr.shape
#import IPython; IPython.embed()
#p.show()

corr = corr/NORM
#cuedict = {'26_23':0.125, '26_38': 0.042, '26_50': 0.083,'26_26':0., '50_57':0.122}
cuedict = {'26_26':0.,'26_38': 0.032557,'26_50':0.073557,'13_32':0.030557,'13_14':0.066557,'50_59':0.071557}
try: ver = cuedict[str(bl1[1])+'_'+str(bl2[1])]
except(KeyError): ver = 0.
meancorr = n.mean(corr,axis=0)
#meancorr = meancorr/REDNORMm
maxind = n.argmax(n.abs(meancorr))
absmax = n.abs(meancorr[maxind])
print '############## baseline_toff RESULT for', bl1, bl2, '#####################'
print "Peak: ", meancorr[maxind], 'abs=',absmax/TT.size, 'at dT = ', TT[n.argmax(n.abs(meancorr))]-2455700.5
meancorr = meancorr/absmax
#import IPython; IPython.embed()
blstr = str(bl1[0])+'_'+str(bl1[1])+'_'+str(bl2[0])+'_'+str(bl2[1])
n.savez('blout_'+blstr, TT-2455700.5,meancorr)

p.figure()
p.subplot(211)
#p.plot(TT-2455700.5,n.abs(n.mean(corr,axis=0)))
p.plot(TT-2455700.5,n.real(meancorr))
p.plot(TT-2455700.5,n.imag(meancorr))
p.axvline(ver,color='k',alpha=0.5,linewidth=3)
p.grid()
p.subplot(212)
#p.plot(TT-2455700.5,n.abs(n.mean(corr,axis=0)))
p.plot(TT-2455700.5,n.abs(meancorr))
p.axvline(ver,color='k',alpha=0.5,linewidth=3)
p.grid()
p.savefig('bl_ton_10000')


