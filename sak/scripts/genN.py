from matplotlib import pyplot as plt
import numpy as np, sys, os, aipy, capo as C, glob
import subprocess

def get_varspec(d1,d2,f1,f2,phase=True,mweight=False,single=False,plot=False):
    """
    Given visibilities d1 and d2 [ndarrays complex128] for two redundant
    baselines and their associated flags f1 and f2 [ndarrays boolean] 
    calculate the variance of their phase product over time.
    
    Input: Four ndarrays [time x freq]
    Output: ndarray [freq]
    """
    d12 = d2*np.conj(d1)
    w12 = f1*f2 #still not comfortable with multiplying bools
    if not phase: d12_var = np.var(d12,axis=0)
    else: d12_var = np.var(np.angle(d12),axis=0)
    d12_wgt = np.sum(np.logical_not(w12),axis=0)
    d12_var *= d12_wgt/d12_wgt.max() #Apply weighting
    d12_mean = np.mean(np.abs(d12),axis=0)
    d12_mean *= d12_wgt/d12_wgt.max() #Apply weighting
    if plot:
        f,axarr=plt.subplots(1,3)
        axarr[0].imshow(np.angle(d12),aspect='auto',interpolation=None)
        axarr[1].plot(range(203),d12_wgt)
        axarr[2].plot(range(203),d12_var)
        plt.show()
        plt.close()
    if not single: 
        if not mweight: return d12_var
        else: return d12_var/d12_mean
    else:
        if not mweight: return np.nanmean(d12_var)
        else: return np.nanmean(d12_var)/np.nanmean(d12_mean)

proc = subprocess.Popen('grid2ant.py -C psa6622_v003 --seps="0,1;0,2;0,3;0,4;0,5;0,6;0,7;0,8;0,9;0,10;0,11;0,12;0,13;0,14;0,15"', shell=True, stdout=subprocess.PIPE)
bls = proc.communicate()[0].rstrip()
sbls = bls.split(',')

file='/home/saulkohn/PAPER_TEST_DATA/S1E1/zen.2456632.50090.xx.uvcRREA'
pol = 'xx'
#bas = [33,38,42,43,47,89,90,107,15,2,31,58,6,88,91,10,105,22,49,64,72,73,97]
datadict = {}

#get data
print 'Getting %i visibilities'%len(sbls)
_,d,f = C.arp.get_dict_of_uv_data([file],antstr=bls,polstr=pol)
print '...Done'

for bl1 in sbls:
    for bl2 in sbls:
        if bl1==bl2: continue
        #every code wants it's own sorta key...
        tbl1=tuple(map(int,bl1.split('_')))
        tbl2=tuple(map(int,bl2.split('_')))
        #if a1 in TESTbas or a2 in TESTbas or b1 in TESTbas or b2 in TESTbas: continue 
        datadict[(bl1,bl2)] = get_varspec(d[tbl1][pol],d[tbl2][pol],\
                                          f[tbl1][pol],f[tbl2][pol],\
                                          mweight=True,single=True)
import IPython;IPython.embed()
