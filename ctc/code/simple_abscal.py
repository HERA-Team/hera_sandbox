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
o.set_usage('simple_abscal.py *uvcRREO')
o.set_description(__doc__)
o.add_option('--lst',dest='lst',default='1-4',type='string',
            help='LST range to extract from data. Default is 1-4.')
o.add_option('--plot',dest='plot',default=False,action="store_true")
o.add_option('--factor',dest='factor',default=None,type='float',
            help='Factor multiplied by data for absolute calibration.')
o.add_option('--abscal',dest='abscal',default=False,action="store_true")
opts,args = o.parse_args(sys.argv[1:])


### Read in data to find matching LST-ranges and factor ###

if opts.factor == None:
    lstmin,lstmax = opts.lst.split('-')
    lstmin=numpy.float(lstmin)*numpy.pi/12 #hours to rad
    lstmax=numpy.float(lstmax)*numpy.pi/12
    print 'lstmin =',lstmin,'rad'
    print 'lstmax=',lstmax,'rad'

    aa = aipy.cal.get_aa('psa898_v003',0.001,0.1,203) #parameters don't matter... only used to find LSTs
    data128 = numpy.sort(glob.glob('/data4/paper/2014EoR/Analysis/ProcessedData/epoch3/omni_v2/lstbin_noxtalk/even/*uvL')) #LST-binned, FRF data
    data64 = numpy.sort(glob.glob('/data4/paper/2012EoR/psa_live/forlstbinning_omnical_2/lstbin_even_noxtalk/sep0,1/*uvGL')) #LST-binned, FRF, abscal Data

    good_files_128 = [] #read in 128 data
    lst_files_128 = []
    for file in data128:
        jd = file.split('.')[1]+'.'+file.split('.')[2]
        aa.set_jultime(float(jd))
        lst = aa.sidereal_time() #[rad] when used but [hrs] when printed
        if lst > lstmin and lst < lstmax:
            good_files_128.append(file)
            lst_files_128.append(lst)
    lst_files_128 = numpy.array(lst_files_128)
    print 'Reading 128 Data in LST Range:',len(numpy.array(good_files_128)),'files'
    t1,d1,f1 = capo.arp.get_dict_of_uv_data(good_files_128,antstr='64_49',polstr='xx',return_lsts=True)
    d1 = d1[aipy.miriad.ij2bl(64,49)]['xx']
    plt.subplot(2,2,1)
    capo.arp.waterfall(d1,mode='log',extent=(0,202,lst_files_128.max()*12/numpy.pi,lst_files_128.min()*12/numpy.pi))
    plt.colorbar()
    plt.ylabel('LST Hours')
    plt.xlabel('Freq Chan')
    plt.title('128 Data (Not absolute calibrated)')

    good_files_64 = [] #read in 64 data 
    lst_files_64 = []
    for file in data64:
        longjd = file.split('/')[-1]
        jd = longjd.split('.')[-3]+'.'+longjd.split('.')[-2]
        aa.set_jultime(float(jd))
        lst = aa.sidereal_time() #[rad] when used but [hrs] when printed
        if lst > lstmin and lst < lstmax:
            good_files_64.append(file)
            lst_files_64.append(lst)
    lst_files_64 = numpy.array(lst_files_64)
    print 'Reading 64 Data in LST Range:',len(numpy.array(good_files_64)),'files'
    bl = '24_48'
    t2,d2,f2 = capo.arp.get_dict_of_uv_data(good_files_64,antstr=bl,polstr='I',return_lsts=True)
    d2 = d2[aipy.miriad.ij2bl(int(bl.split('_')[0]),int(bl.split('_')[1]))]['I'] #24_48
    plt.subplot(2,2,2)
    capo.arp.waterfall(d2,mx=3,drng=4,mode='log',extent=(0,202,lst_files_64.max()*12/numpy.pi,lst_files_64.min()*12/numpy.pi))
    plt.colorbar()
    plt.ylabel('LST Hours')
    plt.xlabel('Freq Chan')
    plt.title('64 Data (Absolute calibrated)')

    freq_chan = 125 #hard-coded frequency channel to get factor from
    d2_f = d2[:,freq_chan]
    d1_f = d1[:,freq_chan]
    factors = d2_f/d1_f
    factor = numpy.mean(factors)
    factor = numpy.real(factor)
    print 'YOUR FACTOR IS:',factor

    print 'Applying Sample Calibration...'
    d1_cal = numpy.zeros_like(d1)
    d1_cal = d1*factor
    plt.subplot(2,2,3)
    capo.arp.waterfall(d1_cal,mx=3,drng=4,mode='log',extent=(0,202,lstmax*12/numpy.pi,lstmin*12/numpy.pi))
    plt.colorbar()
    plt.ylabel('LST Hours')
    plt.xlabel('Freq Chan')
    plt.title('128 Data Calibrated to 64 Data')

if opts.plot == True:
    plt.tight_layout()
    plt.show()

if opts.factor != None: #if factor is given in command-line
    factor = opts.factor


# Absolute calibrate

if opts.abscal == True:

    def mfunc(uv,p,d):
        d *= factor
        return p,d

    for file in args:
        uvi = aipy.miriad.UV(file)
        uvo = aipy.miriad.UV(file+'G', status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi,mfunc=mfunc)
        print file, '->', file+'G'
    





