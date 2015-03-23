#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import glob, sys, os

CHAN = 100
#filelist = glob.glob('*xR'); filelist.sort()
filelist = sys.argv[1:]
bls = '41_49,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_47,3_25,1_58,35_61,42_63,2_33,15_21,8_45,19_47,25_48,1_4,18_35,37_42,6_33,15_16,8_11,19_29,24_48,4_17,5_18,37_40,6_52,16_62,11_36,28_29,24_55,13_17,5_32,14_40,7_52,44_62,36_60,28_34,27_55,13_56,30_32,14_54,7_12,0_44,39_60,34_51,27_57,56_59,23_30,50_54,12_38,0_26,39_46'

for f in filelist:
    npz = {}
    print 'Reading', f
    outfile = os.path.basename(f) + '.npz'
    uv = a.miriad.UV(f)
    a.scripting.uv_selector(uv, ants=bls, pol_str='I')
    for (crd,t,(i,j)),d,_f in uv.all(raw=True):
        bl = a.miriad.ij2bl(i,j)
        sbl = str(bl)
        npz[sbl] = npz.get(sbl,[]) + [d[CHAN]]
        npz['t'+sbl] = npz.get('t'+sbl,[]) + [uv['lst']]
    print 'Writing', outfile
    n.savez(outfile, **npz)
        

#dsort = {}
#for k in data:
#    print k, len(data[k])
#    i = n.argsort(time[k])
#    dsort[str(k)] = n.array(data[k])[i]
#    dsort['t'+str(k)] = n.array(time[k])[i]
    


