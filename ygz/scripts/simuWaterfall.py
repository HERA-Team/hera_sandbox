#loads the data and match times
import numpy as n, aipy as a, capo, os
DIR1 = '/Users/yunfanzhang/local/simuDATA/64_Deltac/0_26/'
DIR2 = '/Users/yunfanzhang/local/simuDATA/64_Deltac/0_38/'
#DIR2 = DIR1
F1 = os.listdir(DIR1)
F2 = os.listdir(DIR2)
for i in range(len(F1)): F1[i] = DIR1+F1[i]
for i in range(len(F2)): F2[i] = DIR2+F2[i]
t026 = 2456249.2666900107
#t038 = 2456249.3086900176
t038 = 2456249.300
#t038 = t026
dt = t038-t026
df = 100./203
bl1, bl2 = (0,26),(0,38)
bl1c, bl2c = a.miriad.ij2bl(*bl1),a.miriad.ij2bl(*bl2)
#freql = n.arange(100,200,df)*1.E6
#phsfac = n.exp(-2*n.pi*3000*n.sin(dt*2*n.pi)*freql/a.phs.const.c*1.j)
T1, dat1, flg1 = capo.arp.get_dict_of_uv_data(F1,antstr='_'.join(map(str,bl1)),polstr='xx')
T2, dat2, flg2 = capo.arp.get_dict_of_uv_data(F2,antstr='_'.join(map(str,bl2)),polstr='xx')
print dat1[283]['xx'].shape

data1 = dat1[bl1c]['xx']
data2 = dat2[bl2c]['xx']  #check 295
print len(T1), len(T2)

while T1[0]+dt>T2[0]: 
    T2 = T2[1:]
    data2 = data2[1:]
while len(T1)>len(T2): 
    T1 = T1[:-1]
    data1 = data1[:-1]

while T1[0]+dt<T2[0]: 
    T1 = T1[1:]
    data1 = data1[1:]
while len(T2)>len(T1): 
    T2 = T2[:-1]
    data2 = data2[:-1]
#for i in range(data2.shape[0]):
#    data2[i] = n.multiply(data2[i],phsfac)
dataave = n.mean(data1*data2.conj(),axis=0)
print dataave.shape

#Phase data to original source
uv1 = a.miriad.UV(DIR1+'pspec_2456249.26465.uv/')
#uv2 = a.miriad.UV(DIR2+'pspec_2456249.30536.uv/')
uv2 = uv1
polstr = 'xx'
nchan = 203
schan = 0
#freqlist = n.array((uv1['sfreq'] + uv1['sdf']*schan + uv1['sdf']*n.arange(nchan)))  #GHz
freqlist = a.cal.get_freqs(uv1['sdf'],uv1['sfreq'],uv1['nchan'])  #GHz
aa = a.cal.get_aa('psa6240_v003',freqlist)
aa.set_active_pol(polstr)
aa.set_jultime(t026)      #must set to exact time determined with the same src and aa
src = a.phs.RadioFixedBody(0, aa.lat)
#phs1, phs2 = n.zeros(data1.shape,dtype='complex64'), n.zeros(data2.shape,dtype='complex64')
#ind = 0
src.compute(aa)
#print aa.gen_uvw(0,26,src)
#print aa.gen_uvw(0,26,'z')
#import IPython; IPython.embed()
phs1 = aa.gen_phs(src, *bl1); phs1.shape = (1,nchan)
aa.set_jultime(t038)
src.compute(aa)
phs2 = aa.gen_phs(src, *bl2); phs2.shape = (1,nchan)
print data2.shape
print phs2.shape
dapa1, dapa2 = data1*phs1,data2*phs2
#dapa1, dapa2 = n.multiply(data1,phs1.conj()),n.multiply(data2,phs2.conj())
#print data2

import matplotlib.pyplot as P

fig = P.figure()
ax = fig.add_subplot(331)
capo.arp.waterfall(data1,mode='real')
ax.set_title("data1")
#ax.set_xlabel("channel")
ax = fig.add_subplot(332)
ax.set_title("data2")
capo.arp.waterfall(data2,mode='real')
ax = fig.add_subplot(334)
capo.arp.waterfall(phs1,mode='phs')
ax.set_title("phs1")
#ax.set_xlabel("channel")
ax = fig.add_subplot(335)
ax.set_title("phs2")
capo.arp.waterfall(phs2,mode='phs')
ax = fig.add_subplot(336)
ax.set_title("phs1phs2")
capo.arp.waterfall(phs1*phs2.conj(),mode='phs')
#P.colorbar()
ax = fig.add_subplot(337)
capo.arp.waterfall(dapa1,mode='real')
ax.set_title("d1phs1")
#ax.set_xlabel("channel")
ax = fig.add_subplot(338)
ax.set_title("d2phs2")
capo.arp.waterfall(dapa2,mode='real')
ax = fig.add_subplot(333)
capo.arp.waterfall(data1*data2.conj(),mode='real')
P.colorbar()
#ax = fig.add_subplot(255)
#capo.arp.waterfall(data1*data2.conj())
ax = fig.add_subplot(339)
capo.arp.waterfall(dapa1*dapa2.conj(),mode='phs')
P.colorbar()
#ax = fig.add_subplot(2510)
#capo.arp.waterfall(dapa1*dapa2.conj())
P.show()
