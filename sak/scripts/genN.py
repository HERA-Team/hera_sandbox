from matplotlib import pyplot as plt
import numpy as np, sys, os, aipy, capo as C, glob
import subprocess

def get_varspec(d1,d2,f1,f2,phase=True,single=False,plot=False):
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
    d12_var *= (d12_wgt/d12_wgt.max()) #Apply weighting
    if plot:
        f,axarr=plt.subplots(1,3)
        axarr[0].imshow(np.angle(d12),aspect='auto',interpolation=None)
        axarr[1].plot(range(203),d12_wgt)
        axarr[2].plot(range(203),d12_var)
        plt.show()
        plt.close()
    if not single: return d12_var
    else: return np.nanmean(d12_var)

proc = subprocess.Popen('grid2ant.py -C psa6622_v003 --seps="0,1;0,2;0,3;0,4;0,5;0,6;0,7;0,8;0,9;0,10;0,11;0,12;0,13;0,14;0,15"', shell=True, stdout=subprocess.PIPE)
bls = proc.communicate()[0]
bls = bls.split(',')

TESTfile='/home/saulkohn/PAPER_TEST_DATA/S1E1/zen.2456632.50090.xx.uvcRREA'
TESTpol = 'xx'
TESTbas = [33,38,42,43,47,89,90,107,15,2,31,58,6,88,91,10,105,22,49,64,72,73,97]
TESTbls = bls[:int(sys.argv[1])]
datadict = {}
for bl1 in TESTbls:
    for bl2 in TESTbls:
        if bl1==bl2: continue
        #every code wants it's own sorta key...
        a1,b1=tbl1=tuple(map(int,bl1.split('_')))
        a2,b2=tbl2=tuple(map(int,bl2.split('_')))
        #if a1 in TESTbas or a2 in TESTbas or b1 in TESTbas or b2 in TESTbas: continue
        
        antstr = ','.join([bl1,bl2])
        _,d,f = C.arp.get_dict_of_uv_data([TESTfile],antstr=antstr,polstr=TESTpol)
        datadict[(bl1,bl2)] = get_varspec(d[tbl1][TESTpol],d[tbl2][TESTpol],\
                                          f[tbl1][TESTpol],f[tbl2][TESTpol],single=True)

