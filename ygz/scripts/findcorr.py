#loads the data and match times
import numpy as n, aipy as a, capo, os

def fftshift2(arr):
    res = []
    for i in range(arr.shape[0]): res.append(n.fft.fftshift(arr[i]))
    return n.array(res)
def ifftshift2(arr):
    res = []
    for i in range(arr.shape[0]): res.append(n.fft.ifftshift(arr[i]))
    return n.array(res)

def wind(arr,wi='hamming'):
    w = a.dsp.gen_window(arr.shape[0], wi)
    arrw = []
    for i in range(arr.shape[1]): arrw.append(w*arr[:,i])
    return n.array(arrw).T

def comb(T,dat):
    tp = 0
    TP, data = [],[]
    for i in range(len(T)):
        if T[i] != tp:
            tp = T[i]
            TP.append(tp); data.append(dat[i])
    #print len(TP)
    return n.array(TP),n.array(data)

def get_corr_fourier(F1,F2,bl1,bl2,chan=None):
    bl1c, bl2c = a.miriad.ij2bl(*bl1),a.miriad.ij2bl(*bl2)
    T1, dat1, flg1 = capo.arp.get_dict_of_uv_data(F1,antstr='_'.join(map(str,bl1)),polstr='xx')
    T2, dat2, flg2 = capo.arp.get_dict_of_uv_data(F2,antstr='_'.join(map(str,bl2)),polstr='xx')
    #print T2.shape
    data1 = dat1[bl1c]['xx']
    data2 = dat2[bl2c]['xx']  #check 295
    #print data2.shape
    T2,data2 = comb(T2,data2)
    #print data2.shape
    if chan != None:
        try: data1, data2 = data1[:,chan],data2[:,chan]
        except(TypeError):
            try:
                shape0 = (data1.shape[0],len(chan))
                data1, data2 = n.array([data1[:,c] for c in chan]),n.array([data1[:,c] for c in chan])
                data1.shape,data2.shape = shape0, shape0
            except(TypeError):
                print "chan must be integeter or a list"
                return
    data1,data2 = wind(data1), wind(data2)
    d1f, d2f = n.fft.fft(data1,axis=0),n.fft.fft(data2,axis=0)
    #d1f, d2f = n.fft.fftshift(d1f),n.fft.fftshift(d2f)
    #print d1f.shape
    data = n.fft.ifft(d2f*d1f.conj(),axis=0)
    data = ifftshift2(data.T).T
    #print data.shape
    #import IPython; IPython.embed()
    print data.shape
    return data




import matplotlib.pyplot as P

bl1, bl2 = (0,26),(0,26)

DIR1 = '/Users/yunfanzhang/local/simuDATA/64_Deltac/0_'+str(bl1[1])+'/'
DIR2 = '/Users/yunfanzhang/local/simuDATA/64_Deltac/0_'+str(bl2[1])+'/'
#clist = [50,80,100,120,150]
clist = [50,80,100,120,150]
F1 = os.listdir(DIR1)
F2 = os.listdir(DIR2)
for i in range(len(F1)): F1[i] = DIR1+F1[i]
for i in range(len(F2)): F2[i] = DIR2+F2[i]

CORL= get_corr_fourier(F1,F2,bl1,bl2,chan=clist)
Trang = n.arange(-0.5,0.5,1./2014)
drang = range(0,CORL.shape[0])
#drang = range(800,1200)
#print Trang.size
P.figure()
P.subplot(111)

cuedict = {'26_23':0.125, '26_38': 0.042, '26_50': 0.083,'26_26':0.}
ver = cuedict[str(bl1[1])+'_'+str(bl2[1])]
P.axvline(ver,color='k',alpha=0.5,linewidth=5)    #for 26 23
for i in range(len(clist)):
    corl = CORL[drang,i]/n.max(n.abs(CORL[:,i]))
    print clist[i], n.max(n.abs(CORL[:,i]))
    P.plot(Trang[drang],n.abs(corl),label=str(clist[i]))
P.legend()
#import IPython; IPython.embed()
P.show()



