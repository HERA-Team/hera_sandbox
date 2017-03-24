#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p
from scipy import signal
import seaborn as sns
sns.set_context("paper")

#@p.ion()
#fqs = n.linspace(.1,.2,203)
fq = .15
a128_dict1 = {'10':(103,26),'11':(103,38),'12':(103,50),'13':(103,23),'14':(103,59),'1-1':(103,46), '-11':(103,12),'-12':(103,54)}
a128_dict2 = {'20':(0,26),'21':(0,38),'22':(0,50),'23':(0,23)}
a128_dict3 = {'30':(102,26),'31':(102,38),'32':(102,50)}
a128_dict4 = {'40':(44,26),'41':(44,38)}
bl1 = a128_dict1['11']
bl2 = a128_dict1['1-1']
#bl1, bl2 = (0,103),(0,95)
#cuedict = {'26_23':0.125, '26_38': 0.042, '26_50': 0.083,'26_26':0., '50_57':0.122}
cuedict = {'26_26':0.,'26_38': 0.032557,'26_46': -0.034, '26_50':0.073557,'13_32':0.030557,'13_14':0.066557,'50_59':0.071557}
#MAXRES_EQUIV = 0.000150070063408  #peak of bl1=bl2
MAXRES_EQUIV = 1260.7792508                #peak of equivalent Opp
#MAXRES_EQUIV = 1.  #to compute MAXRES_EQUIV
#MAXRES_EQUIV = 2.51664842232e-05 #nside=128
#BLUENORM=0.18755
COMPARE = False
try: ver = cuedict[str(bl1[1])+'_'+str(bl2[1])]
except(KeyError): ver = 0.
print 'DelT = ', ver
T0 = 2455700.5
T1 = n.arange(2455700.3,2455700.701,0.001)
#############################################################################
if COMPARE:
    fname = 'blout_'+str(bl1[0])+'_'+str(bl1[1])+'_'+str(bl2[0])+'_'+str(bl2[1])+'.npz'
    print 'Reading', fname
    file = n.load(fname)
    TT, meancorr = file['arr_0'],file['arr_1']
    #try: meancorr = meancorr/BLUENORM
    #except: pass 
#############################################################################
aa = a.cal.get_aa('psa6622_v003', n.array([fq]))
aa.set_jultime(T0)

h = a.healpix.HealpixMap(nside=64)
#h = a.healpix.HealpixMap(nside=64)
tx,ty,tz = h.px2crd(n.arange(h.map.size), ncrd=3)
bl1x, bl1y, bl1z = aa.get_baseline(bl1[0],bl1[1],'z')
bl2x, bl2y, bl2z = aa.get_baseline(bl2[0],bl2[1],'z')
import IPython; IPython.embed()
#fng1=1
bm = aa[0].bm_response((tx,ty,tz),pol='I')[0]**2#/n.abs(tz)   #tz is the Jacobian
#bm = n.ones_like(tx)
tzsave = tz
#bm = n.where(tz > 0, bm, 0)
#import IPython; IPython.embed()
bm = n.where(tz > 0.001, bm, 0)
# -1<tx<1
#import IPython; IPython.embed()
#print bm.sum()
#bm /= bm.sum()



#Create equatorial coordinates of the first frame T0
top = n.array([tx,ty,tz], dtype=tx.dtype)

aa.set_jultime(T0)
m = n.linalg.inv(aa.eq2top_m)
ex,ey,ez = n.dot(m, top)
eq = n.array([ex,ey,ez], dtype=ex.dtype)

V1,V2 = [],[]
print 'Computing Opp'
for t1 in T1:
    aa.set_jultime(t1)
    m = aa.eq2top_m
    tx,ty,tz = n.dot(m, eq)
    #import IPython; IPython.embed()
    #capo.arp.plothmap_ortho

    bl2_prj = tx*bl2x + ty*bl2y + tz*bl2z
    bl1_prj = tx*bl1x + ty*bl1y + tz*bl1z
    fng1 = n.exp(-2j*n.pi*bl1_prj*fq)
    fng2 = n.exp(-2j*n.pi*bl2_prj*fq)
    #fng2=1
    bm = aa[0].bm_response((tx,ty,tz),pol='I')[0]**2#/n.abs(tz)#*n.abs(tzsave)

    #bm = n.ones_like(tx)
    #bm = n.where(tz > 0, bm, 0)
    bm = n.where(tz > 0.0001, bm, 0)
    #import IPython; IPython.embed()
    #print bm.sum()
    #bm /= bm.sum()
    bm_fng2 = bm * fng2
    bm_fng1 = bm * fng1
    #import IPython; IPython.embed()
    #print bm_fng1.shape, bm_fng2.shape
    #print n.sum(bm_fng1.conj()*bm_fng2)
    #res.append(n.sum(bm_fng1.conj()*bm_fng2))
    V1.append(bm_fng1)
    V2.append(bm_fng2)
V1,V2 = n.array(V1), n.array(V2) 
_V1,_V2 = n.fft.fft(V1,axis=0),n.fft.fft(V2,axis=0)
#import IPython; IPython.embed()
res = n.fft.ifftshift(n.fft.ifft(n.sum(_V2*n.conj(_V1),axis=1)))
###################
res = res/MAXRES_EQUIV/T1.size
#res = res/maxres
###################
maxind = n.argmax(n.abs(res))
maxres = n.abs(res[maxind])

T1ac = -T0+T1[maxind]
print '############## OPP RESULT for', bl1, bl2, '#####################'
print 'max, max adjusted for NT, dT(max)'
print maxres,maxres, T1ac
T1 = T1-T0
p.figure()
p.plot(T1,res.real,'b',label='real')
p.plot(T1,res.imag,'g',label='imag')
p.plot(T1,n.abs(res),'r',alpha=0.5,linewidth=1,label='Theory(Opp)')
if COMPARE:
    meancorr = meancorr/1.257
    p.plot(TT,meancorr.real,'b--')
    p.plot(TT,meancorr.imag,'g--')
    p.plot(TT,n.abs(meancorr),'r--',alpha=0.5,linewidth=1,label='Simulation')
#p.axvline(T1ac,color='k',alpha=0.5,linewidth=2)
p.xlabel('dT (Julian Day)')
p.title('Correlation Normalized to Equivalent Baseline Peak')
p.legend()

p.show()
#import IPython; IPython.embed()