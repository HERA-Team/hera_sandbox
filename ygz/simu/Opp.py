#! /usr/bin/env python
import aipy as a, numpy as n, capo as C, pylab as p
from scipy import signal

#@p.ion()
#fqs = n.linspace(.1,.2,203)
fq = .15
bl1, bl2 = (0,26),(0,26)
#cuedict = {'26_23':0.125, '26_38': 0.042, '26_50': 0.083,'26_26':0., '50_57':0.122}
cuedict = {'26_26':0.,'26_38': 0.032557,'26_46': -0.034, '26_50':0.073557,'13_32':0.030557,'13_14':0.066557,'50_59':0.071557}
#REDNORM = 0.000150070063408  #peak of bl1=bl2
REDNORM = 9.E-05 #nside=64
#REDNORM = 1.  #to compute REDNORM
#REDNORM = 2.51664842232e-05 #nside=128
BLUENORM=0.18755
COMPARE = True
try: ver = cuedict[str(bl1[1])+'_'+str(bl2[1])]
except(KeyError): ver = 0.
print 'DelT = ', ver
T0 = 2455700.5
#dT = n.arange(-0.01,0.01,0.0001)
dT = n.arange(-0.03,0.03,1/a.const.sidereal_day*42.9*0.5)
#dT = n.arange(-0.1,0.1,0.001)
#dT = n.array([0])
T1 = T0+ver+dT
#############################################################################
if COMPARE:
    fname = 'blout_'+str(bl1[0])+'_'+str(bl1[1])+'_'+str(bl2[0])+'_'+str(bl2[1])+'.npz'
    print 'Reading', fname
    file = n.load(fname)
    TT, meancorr = file['arr_0'],file['arr_1']
    try: meancorr = meancorr/BLUENORM
    except: pass 
#############################################################################
aa = a.cal.get_aa('psa6240_v003', n.array([fq]))
aa.set_jultime(T0)

h = a.healpix.HealpixMap(nside=64)
#h = a.healpix.HealpixMap(nside=64)
tx,ty,tz = h.px2crd(n.arange(h.map.size), ncrd=3)
bl1x, bl1y, bl1z = aa.get_baseline(bl1[0],bl1[1],'z')
bl2x, bl2y, bl2z = aa.get_baseline(bl2[0],bl2[1],'z')

bl1_prj = tx*bl1x + ty*bl1y + tz*bl1z
fng1 = n.exp(-2j*n.pi*bl1_prj*fq)
#bm = aa[0].bm_response((tx,ty,tz),pol='x')[0]**2#/n.abs(tz)   #tz is the Jacobian
bm = n.ones_like(tx)
tzsave = tz
#bm = n.where(tz > 0, bm, 0)
bm = n.where(tz > 0.001, bm, 0)
# -1<tx<1
#import IPython; IPython.embed()
#print bm.sum()
bm /= bm.sum()
bm_fng1 = bm * fng1


top = n.array([tx,ty,tz], dtype=tx.dtype)

aa.set_jultime(T0)
m = n.linalg.inv(aa.eq2top_m)
ex,ey,ez = n.dot(m, top)
eq = n.array([ex,ey,ez], dtype=ex.dtype)

res = []
print 'Computing Opp'
for t1 in T1:
    aa.set_jultime(t1)
    m = aa.eq2top_m
    tx,ty,tz = n.dot(m, eq)
    #capo.arp.plothmap_ortho

    bl2_prj = tx*bl2x + ty*bl2y + tz*bl2z
    fng2 = n.exp(-2j*n.pi*bl2_prj*fq)
    #bm = aa[0].bm_response((tx,ty,tz),pol='x')[0]**2#/n.abs(tz)#*n.abs(tzsave)
    bm = n.ones_like(tx)
    #bm = n.where(tz > 0, bm, 0)
    bm = n.where(tz > 0.001, bm, 0)
    #print bm.sum()
    bm /= bm.sum()
    bm_fng2 = bm * fng2

    #import IPython; IPython.embed()
    #print bm_fng1.shape, bm_fng2.shape
    #print n.sum(bm_fng1.conj()*bm_fng2)
    res.append(n.sum(bm_fng1.conj()*bm_fng2))
res = n.array(res)
###################
res = res/REDNORM
###################
maxind = n.argmax(n.abs(res))
T1ac = -T0+T1[maxind]
print '############## OPP RESULT for', bl1, bl2, '#####################'
print 'max, abs(max), dT(max)'
print res[maxind],n.abs(res[maxind]), T1ac
T1 = T1-T0
p.figure()
p.plot(T1,res.real,'b',label='real')
p.plot(T1,res.imag,'g',label='imag')
p.plot(T1,n.abs(res),'r',alpha=0.5,linewidth=1,label='Theory(Opp)')
if COMPARE:
    p.plot(TT,meancorr.real,'b--')
    p.plot(TT,meancorr.imag,'g--')
    p.plot(TT,n.abs(meancorr),'r--',alpha=0.5,linewidth=1,label='Simulation')
p.axvline(T1ac,color='k',alpha=0.5,linewidth=2)
p.xlabel('dT (Julian Day)')
p.title('Correlation Normalized to Equivalent Baseline Peak')
p.legend()
p.grid()

#recovered, remainder = signal.deconvolve(n.abs(meancorr), n.abs(res))
#recovered = n.where(n.abs(recovered)<100, recovered, 0 )
#p.figure()
#p.plot(recovered)
#p.grid()

p.show()
#import IPython; IPython.embed()