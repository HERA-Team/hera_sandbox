#! /usr/bin/env python

import aipy
import numpy
import capo
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pylab
import optparse
import glob 

### Options ###
o = optparse.OptionParser()
o.set_usage('compare_uvs.py [options]')
o.set_description(__doc__)
o.add_option('--lst',dest='lst',default='1-4',type='string',
            help='LST range to extract from data. Default is 1-4.')
o.add_option('--plot',dest='plot',default=False,action="store_true")
o.add_option('--abscal',dest='abscal',default=False,action="store_true")
o.add_option('--mode',dest='mode',default='log',type='string',
            help='Mode to plot data. Default is "log".')
opts,args = o.parse_args(sys.argv[1:])


### Save Options ###
lstmin,lstmax = opts.lst.split('-')
lstmin=numpy.float(lstmin)*numpy.pi/12 #hours to rad
lstmax=numpy.float(lstmax)*numpy.pi/12
print 'lstmin =',lstmin,'rad'
print 'lstmax=',lstmax,'rad'


### Read Files ###
aa = aipy.cal.get_aa('psa898_v003',0.001,0.1,203) #parameters don't matter... only used to find LSTs
data128 = numpy.sort(glob.glob('/data4/paper/2014EoR/pot3/data2/2456928/*xx.uvcRRE')) #compressed 128
#data128 = numpy.sort(glob.glob('/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/good/zen.2456962.*xx.uvcRREOO')) #omnicaled 128
data64 = glob.glob('/data4/paper/2012EoR/psa_live/forlstbinning_omnical_2/zen.2456247.*.uvcRREcACO')
simdata = glob.glob('/home/cacheng/capo/ctc/simpawz_sims/08_21_2015/gsm_Jy.uv')


### Open 128 Data ###
good_files = []
lst_files = []
for file in data128:
    jd = file.split('.')[-4]+'.'+file.split('.')[-3]
    aa.set_jultime(float(jd))
    lst = aa.sidereal_time() #[rad] when used but [hrs] when printed
    if lst > lstmin and lst < lstmax:
        good_files.append(file)
        lst_files.append(lst)
lst_files = numpy.array(lst_files)
print 'Reading 128 Data in LST Range:',len(numpy.array(good_files)),'files'
t1,d1,f1 = capo.arp.get_dict_of_uv_data(good_files,antstr='64_49',polstr='xx',return_lsts=True)
d1 = d1[aipy.miriad.ij2bl(64,49)]['xx']
plt.subplot(2,2,1)
capo.arp.waterfall(d1,drng=3,mode=opts.mode,extent=(0,202,lst_files.max()*12/numpy.pi,lst_files.min()*12/numpy.pi))
plt.colorbar()
plt.ylabel('LST Hours')
plt.xlabel('Freq Chan')
plt.title('128 Data')

### Open 64 Data ###
good_files = []
lst_files = []
for file in data64:
    longjd = file.split('/')[-1]
    jd = longjd.split('.')[-3]+'.'+longjd.split('.')[-2]
    aa.set_jultime(float(jd))
    lst = aa.sidereal_time() #[rad] when used but [hrs] when printed
    if lst > lstmin and lst < lstmax:
        good_files.append(file)
        lst_files.append(lst)
lst_files = numpy.array(lst_files)
print 'Reading 64 Data in LST Range:',len(numpy.array(good_files)),'files'
bl = '24_48'
t2,d2,f2 = capo.arp.get_dict_of_uv_data(good_files,antstr=bl,polstr='xx',return_lsts=True)
d2 = d2[aipy.miriad.ij2bl(int(bl.split('_')[0]),int(bl.split('_')[1]))]['xx'] #24_48
#d2 = numpy.nan_to_num(d2)
plt.subplot(2,2,2)
capo.arp.waterfall(d2,drng=3,mx=3.8,mode=opts.mode,extent=(0,202,lst_files.max()*12/numpy.pi,lst_files.min()*12/numpy.pi))
plt.colorbar()
plt.ylabel('LST Hours')
plt.xlabel('Freq Chan')
plt.title('64 Data')

### Open GSM Simulation ###
print 'Reading Sim Data'
t3,d3,f3 = capo.arp.get_dict_of_uv_data(simdata,antstr='cross',polstr='xx',return_lsts=True)
d3 = d3[aipy.miriad.ij2bl(64,49)]['xx']
d3 = d3[numpy.where(numpy.logical_and(t3<lstmax,t3>lstmin))] #select LST range
lsts = t3[numpy.where(numpy.logical_and(t3<lstmax,t3>lstmin))] #LSTs in range
plt.subplot(2,2,3)
capo.arp.waterfall(d3,drng=3,mx=3.8,mode=opts.mode,extent=(0,202,lsts.max()*12/numpy.pi,lsts.min()*12/numpy.pi))
plt.colorbar()
plt.ylabel('LST Hours')
plt.xlabel('Freq Chan')
plt.title('GSM Sim')


if opts.plot == True and opts.abscal == False:
    plt.tight_layout()
    plt.show()

"""
### Average 64 Data for Good Calibration Model ###
print 'Reading Calfile'
aa = aipy.cal.get_aa('psa6240_v003',numpy.array([.15]))
sep, conj = capo.red.group_redundant_bls(aa.ant_layout)
bl_str = sep['0,1'] #change in row, change in col (ex: '0,1' is all the E/W 30 m bls)
bl_str = ['%d_%d' % aipy.miriad.bl2ij(bl) for bl in bl_str] #list comprehension
print bl_str
bls = ""
for i in bl_str:
    bls += i+','
bls = bls[:-1]
t2,d2,f2 = capo.arp.get_dict_of_uv_data(good_files,bls,polstr='xx',return_lsts=True)
data_allbls = []
for b in range(len(d2.keys())):
    bl = d2.keys()[b]
    if conj[bl] == False:
        dat = numpy.conj(d2[bl]['xx'])
    else:
        dat = d2[bl]['xx']
    data_allbls.append(dat)
data_allbls = numpy.array(data_allbls)
d2 = numpy.mean(data_allbls,axis=0)
"""

### Calibration Against 64-Data ###
if opts.abscal == True:
    print 'Absolute Calibrating...'
    ref_t1_index = len(t1)/2
    ref_t1 = t1[ref_t1_index] 
    print '...calibrating to LST =',ref_t1*12/numpy.pi,'hours =',ref_t1,'rad'
    diff = numpy.abs(t2-ref_t1)
    close_lst1_index = numpy.argmin(diff)
    close_lst1 = t2[close_lst1_index]
    test2_1 = diff[close_lst1_index-1]
    test2_2 = diff[close_lst1_index+1]
    if test2_1 < test2_2:
        close_lst2_index = close_lst1_index-1
        close_lst2 = t2[close_lst1_index-1]
    elif test2_1 > test2_2:
        close_lst2_index = close_lst1_index+1
        close_lst2 = t2[close_lst1_index+1]
    if close_lst1 < close_lst2: #interpolation between LSTs
        dinterp = interp1d([close_lst1,close_lst2],[d2[close_lst1_index],d2[close_lst2_index]],kind='linear',axis=0)
    else:
        dinterp = interp1d([close_lst2,close_lst1],[d2[close_lst2_index],d2[close_lst1_index]],kind='linear',axis=0)
    d2_interp = dinterp(ref_t1)
    d_128 = numpy.ma.masked_where(d1[ref_t1_index]==0,d1[ref_t1_index])
    d_64 = numpy.ma.masked_where(d1[ref_t1_index]==0,d2_interp)
    factors = d_64/d_128
    print 'Applying Sample Calibration...'
    d1_cal = numpy.zeros_like(d1)
    for t in range(len(d1)):
        d1_cal[t] = d1[t]*factors
    plt.subplot(2,2,4)
    capo.arp.waterfall(d1_cal,drng=3,mode=opts.mode,extent=(0,202,lstmax*12/numpy.pi,lstmin*12/numpy.pi),mx=3.8)
    plt.colorbar()
    plt.ylabel('LST Hours')
    plt.xlabel('Freq Chan')
    plt.title('128 Data Calibrated to 64 Data')

if opts.plot == True and opts.abscal == True:
    plt.tight_layout()
    plt.show()

"""
### Absolute Cal on 60-m Baseline ###
print args
t,d,f = capo.arp.get_dict_of_uv_data(args,antstr='64_41',polstr='xx',return_lsts=True)
print t
d = d[aipy.miriad.ij2bl(64,41)]['xx']
d_cal = numpy.zeros_like(d)
for t in range(len(d)):
    d_cal[t] = d[t]*factors
capo.arp.waterfall(d_cal,drng=3,mode=opts.mode,mx=3.8)
plt.colorbar()
plt.show()
"""
