#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import glob, sys, os,optparse

o = optparse.OptionParser()
o.add_option('-c', "--chan", help="Channel to pull [default 100]",type=int)
o.add_option('--bls',help='Baselines, comma delim eg "41_64,65_66" default is all bls of type 0N,2E in the psa128 array')
o.add_option('-p','--pol',help='polarization. For now, you must pick one of: xx,yy,xy or yx')
opts, args = o.parse_args()
CHAN = opts.chan
#filelist = glob.glob('*xR'); filelist.sort()
filelist = sys.argv[1:]
#bls = '41_49,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_47,3_25,1_58,35_61,42_63,2_33,15_21,8_45,19_47,25_48,1_4,18_35,37_42,6_33,15_16,8_11,19_29,24_48,4_17,5_18,37_40,6_52,16_62,11_36,28_29,24_55,13_17,5_32,14_40,7_52,44_62,36_60,28_34,27_55,13_56,30_32,14_54,7_12,0_44,39_60,34_51,27_57,56_59,23_30,50_54,12_38,0_26,39_46'
#psa128 0,2 spacings
if not opts.bls is None:
    bls=opts.bls
else:
    bls='49_64,65_66,72_73,80_81,88_89,96_97,104_105,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_49,66_67,73_74,81_82,89_90,97_98,105_106,3_25,1_58,35_61,42_63,2_33,15_21,8_45,19_41,47_67,74_75,82_83,90_91,98_99,106_107,25_48,1_4,18_35,37_42,6_33,15_16,8_11,19_29,47_68,75_76,83_84,91_92,99_100,107_108,24_48,4_17,5_18,37_40,6_52,16_62,11_36,28_29,68_69,76_77,84_85,92_93,100_101,108_109,24_55,13_17,5_32,14_40,7_52,44_62,36_60,28_34,69_70,77_78,85_86,93_94,101_102,109_110,27_55,13_56,30_32,14_54,7_12,0_44,39_60,34_51,70_71,78_79,86_87,94_95,102_103,110_111,27_57,56_59,23_30,50_54,12_38,0_26,39_46'

print "using bls:",bls
for f in args:
    npz = {}
    print 'Reading', f
    outfile = os.path.basename(f) + '.npz'
    if os.path.exists(outfile):
        print outfile, "exists, skipping"
    uv = a.miriad.UV(f)
    a.scripting.uv_selector(uv, ants=bls, pol_str=opts.pol)
    for (crd,t,(i,j)),d,_f in uv.all(raw=True):
        bl = a.miriad.ij2bl(i,j)
        sbl = str(bl)
        if len(d)<CHAN: print "ERROR Channel not found in data range: 0",len(d)
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
    


