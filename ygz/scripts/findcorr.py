#loads the data and match times
import numpy as n, aipy as a, capo, os


def get_corr(t1,t2,F1,F2,bl1,bl2,chan=None,phsd=False):
    #t038 = t026
    dt = t2-t1
    df = 100./203
    bl1c, bl2c = a.miriad.ij2bl(*bl1),a.miriad.ij2bl(*bl2)
    #freql = n.arange(100,200,df)*1.E6
    #phsfac = n.exp(-2*n.pi*3000*n.sin(dt*2*n.pi)*freql/a.phs.const.c*1.j)
    T1, dat1, flg1 = capo.arp.get_dict_of_uv_data(F1,antstr='_'.join(map(str,bl1)),polstr='xx')
    T2, dat2, flg2 = capo.arp.get_dict_of_uv_data(F2,antstr='_'.join(map(str,bl2)),polstr='xx')
    if T1[0]+dt>T2[-1] or T2[0]>T1[-1]+dt:
        print "input ignored, outside %f <dt< %f" % (T2[0]-T1[-1],T2[-1]-T1[0])
        return None
    data1 = dat1[bl1c]['xx']
    data2 = dat2[bl2c]['xx']  #check 295
    if chan != None: data1, data2 = data1[:,chan],data2[:,chan]
    #print data2.shape
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
    #print len(T1), len(T2)
    #print data2.shape
    if phsd:
        #Phase data to original source
        #uv1 = a.miriad.UV(DIR1+'pspec_2456249.26525.uv/')
        #uv2 = a.miriad.UV(DIR2+'pspec_2456249.30497.uv/')
        uv1 = a.miriad.UV(DIR1+'pspec_2456249.26465.uv/')
        #uv2 = a.miriad.UV(DIR2+'pspec_2456249.30536.uv/')
        #uv2 = uv1
        polstr = 'xx'
        nchan = 203
        schan = 0
        freqlist = a.cal.get_freqs(uv1['sdf'],uv1['sfreq'],uv1['nchan'])  #GHz
        if chan != None:
            freqlist = n.array([freqlist[chan]])
            nchan = 1
        aa = a.cal.get_aa('psa6240_v003',freqlist)
        aa.set_active_pol(polstr)
        aa.set_jultime(t026)      #must set to exact time determined with the same src and aa
        #src = a.phs.RadioFixedBody(aa.sidereal_time(), aa.lat, epoch=aa.epoch)
        src = a.phs.RadioFixedBody(0, aa.lat)
        src.compute(aa)
        phs1 = aa.gen_phs(src, *bl1); phs1.shape = (1,nchan)
        aa.set_jultime(t038)
        src.compute(aa)
        phs2 = aa.gen_phs(src, *bl2); phs2.shape = (1,nchan)
        data1, data2 = data1*phs1,data2*phs2
    #corr = n.mean(data1*data2.conj())
    #print data2.shape
    #print (data1*data2.conj()).shape
    corr = data1*data2.conj()
    #print corr.shape
    corr = n.mean(corr)
    return corr

import matplotlib.pyplot as P
#DIR1 = '/Users/yunfanzhang/local/simuDATA/Deltatest/0_26/'
#DIR2 = '/Users/yunfanzhang/local/simuDATA/Deltatest/0_38/'
DIR1 = '/Users/yunfanzhang/local/simuDATA/64_Deltac/0_26/'
DIR2 = '/Users/yunfanzhang/local/simuDATA/64_Deltac/0_38/'
#single bl test
#DIR2 = DIR1
bl1, bl2 = (0,26),(0,38)

F1 = os.listdir(DIR1)
F2 = os.listdir(DIR2)
for i in range(len(F1)): F1[i] = DIR1+F1[i]
for i in range(len(F2)): F2[i] = DIR2+F2[i]

t026 = 2456249.2666900107
corl,corp = [],[]
Trang = n.arange(2456249.15,2456249.35,0.01)
for t038 in Trang:
    co,cp = get_corr(t026,t038,F1,F2,bl1,bl2,chan=140),get_corr(t026,t038,F1,F2,bl1,bl2,chan=140,phsd=True)
    if co ==None or cp==None:
        n.delete(Trang,t038)
        break
    #co= get_corr(t026,t038,F1,F2,bl1,bl2)
    print t038, co,cp
    corl.append(co)
    corp.append(cp)
corl,corp = n.array(corl),n.array(corp)
#print Trang[n.argmax(corp.real)],n.max(corp.real)
P.figure()
P.subplot(121)
P.plot(Trang,n.abs(corl),label='raw')
P.plot(Trang,n.abs(corp),label='phsd')
P.legend()
P.subplot(122)
P.plot(Trang, n.angle(corl),label='raw')
P.plot(Trang, n.angle(corp),label='phsd')
P.legend()
P.show()



