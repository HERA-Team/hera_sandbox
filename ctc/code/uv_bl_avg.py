#! /usr/bin/env python

import aipy
import numpy
import capo
import matplotlib.pyplot as plt
import os,sys
import optparse

#####
# Averages similar types of baselines together in a UV file
#####


### Options ###
o = optparse.OptionParser()
o.set_usage('uv_bl_avg.py *uvcRREO')
o.set_description(__doc__)
o.add_option('--psa64',dest='psa64',default=False,action="store_true",
            help="Application is for PSA64 data.")
o.add_option('--psa128',dest='psa128',default=False,action="store_true",
            help="Application is for PSA128 data.")
opts,args = o.parse_args(sys.argv[1:])

### seps 0,2 ### 
if opts.psa128:
    ANTS='49_64,65_66,72_73,80_81,88_89,96_97,104_105,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_49,66_67,73_74,81_82,89_90,97_98,105_106,3_25,1_58,35_61,42_63,2_33,15_21,8_45,19_41,47_67,74_75,82_83,90_91,98_99,106_107,25_48,1_4,18_35,37_42,6_33,15_16,8_11,19_29,47_68,75_76,83_84,91_92,99_100,107_108,24_48,4_17,5_18,37_40,6_52,16_62,11_36,28_29,68_69,76_77,84_85,92_93,100_101,108_109,24_55,13_17,5_32,14_40,7_52,44_62,36_60,28_34,69_70,77_78,85_86,93_94,101_102,109_110,27_55,13_56,30_32,14_54,7_12,0_44,39_60,34_51,70_71,78_79,86_87,94_95,102_103,110_111,27_57,56_59,23_30,50_54,12_38,0_26,39_46' #PSA128 30-m baselines
    pol='xx'
    aa = aipy.cal.get_aa('psa6622_v003',0.001,0.1,203)
if opts.psa64:
    ANTS='41_49,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_47,3_25,1_58,35_61,42_63,2_33,15_21,8_45,19_47,25_48,1_4,18_35,37_42,6_33,15_16,8_11,19_29,24_48,4_17,5_18,37_40,6_52,16_62,11_36,28_29,24_55,13_17,5_32,14_40,7_52,44_62,36_60,28_34,27_55,13_56,30_32,14_54,7_12,0_44,39_60,34_51,27_57,56_59,23_30,50_54,12_38,0_26,39_46' #PSA64 30-m baselines
    pol = 'I'
    aa = aipy.cal.get_aa('psa6240_v003',0.001,0.1,203)
if opts.psa128==False and opts.psa64==False:
    print "Use option --psa64 or --psa128"
    sys.exit()

def mfunc(uv,p,d):
    mask = d.mask
    d = numpy.ma.masked_array(avg_d[p[1]],mask) #d needs to be a masked array
    return p,d

bl_str,bl_conj,bl2sep_str = capo.zsa.grid2ij(aa.ant_layout)

for file in args:
    print file, '->', file+'V'
    if os.path.exists(file+'V'):
        print '    File exists... skipping.'
        continue
    uvi = aipy.miriad.UV(file)
    t,d,f = capo.arp.get_dict_of_uv_data([file],antstr=ANTS,polstr=pol)
    avg_d = {}
    for tt,time in enumerate(t):
        d_bls = [] #data for the baselines
        for bl in d.keys():
            if bl_conj[bl]: #conjugate data if needed
                d_bl = d[bl][pol][tt].conj()
            else:
                d_bl = d[bl][pol][tt]
            d_bls.append(d_bl)
        avg_d[time] = numpy.mean(numpy.array(d_bls),axis=0) #avg vis over baselines
    aipy.scripting.uv_selector(uvi, ants='2_33') #only one baseline
    uvo = aipy.miriad.UV(file+'V', status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc,append2hist='UV_BL_AVG:'+' '.join(sys.argv)+'\n')

