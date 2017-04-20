#! /usr/bin/env python

import os
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
o.add_option('--plot',dest='plot',default=False,action="store_true")
o.add_option('--factor',dest='factor',default=None,type='string',
            help='Name of npz file containing bandpass.')
o.add_option('--abscal',dest='abscal',default=False,action="store_true")
o.add_option('--poly',dest='poly',default=False,action="store_true",
            help="Fit polynomial vs. frequency instead of single number.")
o.add_option('--pols',dest='pols',default='xx,xx',
            help='Pol of 64 data, pol of 128 data.')
o.add_option('--bls',dest='bls',default='1_4,1_4',
            help='Baseline of 64 data, baseline of 128 data.')
opts,args = o.parse_args(sys.argv[1:])


if opts.factor == None:

    ### FILES ###
    aa = aipy.cal.get_aa('psa898_v003',0.001,0.1,203) #parameters don't matter... only used to find LSTs
    data128 = numpy.sort(glob.glob('/data4/paper/2013EoR/Analysis/ProcessedData/epoch1/omni_v4_xtalk/lstbin_balanced_fg/even/lst*I.uv')) #S1 PSA-128, FG-containing
    #data128 = numpy.sort(glob.glob('/data4/paper/2013EoR/Analysis/ProcessedData/epoch2/omni_v2_xtalk/lstbin_fg/even/lst*I.uv')) #LST-binned, FG-containing
    #data64 = numpy.sort(glob.glob('/home/jacobsda/storage/psa128/2014_epoch3/v5_xtalksub_omni/lstbin_June2_v1/even/sep0,2/*uvAS')) #LST-binned, abscal 128 !! 
    #data64 = numpy.sort(glob.glob('/data4/paper/2012EoR/psa_live/forlstbinning_omnical_2/lstbin_even_noxtalk/*uvG')) #LST-binned
    data64 = numpy.sort(glob.glob('/data4/paper/2012EoR/psa_live/forlstbinning_omnical_2/lstbin_fg_even/lst.*.uvA')) #LST-binned, FG-containing    
    
    ### 64-DATA ###    
    print 'Reading 64 Data:',len(numpy.array(data64)),'files'
    t1,d1,f1 = capo.arp.get_dict_of_uv_data(data64,antstr=opts.bls.split(',')[0],polstr=opts.pols.split(',')[0]) 
    i,j = int(opts.bls.split(',')[0].split('_')[0]),int(opts.bls.split(',')[0].split('_')[1])
    d1 = d1[(i, j)][opts.pols.split(',')[0]]
    t1 = t1['lsts']
    order = numpy.argsort(t1)
    t1 = t1[order]
    d1 = d1[order]

    ### 128-DATA ###
    print 'Reading 128 Data:',len(numpy.array(data128)),'files'
    t2,d2,f2 = capo.arp.get_dict_of_uv_data(data128,antstr=opts.bls.split(',')[1],polstr=opts.pols.split(',')[1]) 
    i,j = int(opts.bls.split(',')[1].split('_')[0]),int(opts.bls.split(',')[1].split('_')[1])
    d2 = d2[(i, j)][opts.pols.split(',')[1]] 
    t2 = t2['lsts']
    order = numpy.argsort(t2) # sort if LSTs are mixed up
    t2 = t2[order]
    d2 = d2[order]
    
    ### Make sure LST-range for 64-data is larger than 128-data ###
    if numpy.min(t2) < numpy.min(t1):
        print 'Data range for 128 data too large: clipping LSTs below',numpy.min(t1)
        t2_clip = [t for t in t2 if t > numpy.min(t1)]
        d2_clip = [d for i,d in enumerate(d2) if t2[i] > numpy.min(t1)]
        t2 = t2_clip
        d2 = d2_clip
    if numpy.max(t2) > numpy.max(t1):
        print 'Data range for 128 data too large: clipping LSTs above',numpy.max(t1)
        t2_clip = [t for t in t2 if t < numpy.max(t1)]
        d2_clip = [d for i,d in enumerate(d2) if t2[i] < numpy.max(t1)]
        t2 = t2_clip
        d2 = d2_clip
    ### FIND MATCHING LSTs and FACTORS ###
    factors = []
    factors_complex = []
    d64_all = []
    d128_all = []
    for i2,t in enumerate(t2):
        i1 = numpy.argmin(numpy.abs(t1-t)) #XXX could interpolate here
        d64_complex = d1[i1]
        d64_complex[numpy.abs(d64_complex) == 0.] = numpy.nan
        d64 = numpy.abs(d1[i1])
        d64[d64 == 0.] = numpy.nan #make 0 values in PSA64 nan's
        #d64[:30] = numpy.nan #chop off low band edge in polyfit
        #d64[-25:] = numpy.nan #chop off high band edge in polyfit
        d128_complex = d2[i2]
        d128_complex[numpy.abs(d128_complex) == 0.] = numpy.nan
        d128 = numpy.abs(d2[i2])
        d64_all.append(d64_complex)
        d128_all.append(d128_complex)
        factors.append(d64/d128)
        factors_complex.append(d64_complex/d128_complex)
    factors = numpy.ma.masked_invalid(factors) #mask inf and nan
    factors_complex = numpy.ma.masked_invalid(factors_complex)
    d64_all = numpy.ma.masked_invalid(d64_all)
    d128_all = numpy.ma.masked_invalid(d128_all)
    if opts.poly == True:
        print "Fitting polynomial..."
        freqs = numpy.linspace(0,202,203)
        capo.arp.waterfall(factors,drng=4)
        plt.title('PSA64/PSA128 (Abs) Pre-Calibration')
        plt.xlabel('freq')
        plt.ylabel('time')
        plt.colorbar()
        if opts.plot == True:
            plt.show()
        plt.imshow(numpy.angle(factors_complex),aspect='auto')
        plt.title('PSA64/PSA128 (Phase) Pre-Calibration')
        plt.xlabel('freq')
        plt.ylabel('time')
        plt.colorbar()
        if opts.plot == True:
            plt.show()
        factors = numpy.ma.median(factors,axis=0)
        #factors[numpy.where(numpy.isfinite(factors.data) == False)] = 0.0
        factors_complex = numpy.ma.median(factors_complex,axis=0)
        #factors_complex[numpy.where(numpy.isfinite(factors_complex.data) == False)] = 0.0
        freqs = numpy.ma.masked_where(factors.mask == True,freqs)
        factor = numpy.polyval(numpy.ma.polyfit(freqs,factors,8),freqs)
        plt.plot(factors,'r-',label='PSA64/PSA128 (Abs)')
        plt.plot(factor,'k-',label='fit')
        plt.title('polynomial fit (abs)')
        plt.legend()
        if opts.plot == True:
            plt.show()
        starti = 20
        endi = 170 #to fit middle section of band only
        phase_factor = numpy.polyval(numpy.ma.polyfit(freqs[starti:endi],numpy.angle(factors_complex)[starti:endi],1),freqs) #linear fit for phase term
        plt.plot(numpy.angle(factors_complex),'r-',label='PSA64/PSA128 (Phase)')
        plt.plot(phase_factor,'k-',label='fit')
        plt.title('polynomial fit (phase)')
        plt.legend()
        if opts.plot == True:
            plt.show()
        plt.imshow(numpy.angle(d64_all/(d128_all*numpy.exp(1j*phase_factor))),aspect='auto')
        plt.title('PSA64/PSA128 (Phase) After-Calibration')
        plt.xlabel('freq')
        plt.ylabel('time')
        plt.colorbar()
        if opts.plot == True:
            plt.show()
        factor = factor*numpy.exp(1j*phase_factor) #final correction factor
    else:
        factors = numpy.mean(factors,axis=1) #average over freq
        factor = numpy.median(factors) #median over time
        print "Factor =",factor

### PLOT ###
if opts.plot == True:
    plt.subplot(2,3,1)
    capo.arp.waterfall(d2,mx=0,drng=4,mode='log')
    plt.colorbar()
    plt.title('128 Data')
    plt.subplot(2,3,2)
    t1_clip = [t for t in t1 if t > numpy.min(t2) and t < numpy.max(t2)]
    d1 = [d for i,d in enumerate(d1) if t1[i] > numpy.min(t2) and t1[i] < numpy.max(t2)]
    capo.arp.waterfall(d1,mx=4,drng=4,mode='log')
    plt.colorbar()
    plt.title('64 Data')
    plt.subplot(2,3,3)
    capo.arp.waterfall(numpy.array(d2)*factor,mx=4,drng=4,mode='log')
    plt.colorbar()
    plt.title('Calibrated 128 Data')
    plt.subplot(2,3,4)
    plt.imshow(numpy.angle(d2),aspect='auto')
    plt.colorbar()
    plt.subplot(2,3,5)
    plt.imshow(numpy.angle(d1),aspect='auto')
    plt.colorbar()
    plt.subplot(2,3,6)
    plt.imshow(numpy.angle(numpy.array(d2)*factor),aspect='auto')
    plt.colorbar()
    plt.tight_layout()
    plt.show()


"""
    good_files_128 = [] #read in 128 data
    lst_files_128 = []
    for file in data128:
        jd = file.split('.')[1]+'.'+file.split('.')[2]
        aa.set_jultime(float(jd))
        lst = aa.sidereal_time() #[rad] when used but [hrs] when printed
        if lst > lstmin*numpy.pi/12 and lst < lstmax*numpy.pi/12:
            good_files_128.append(file)
            lst_files_128.append(lst)
    lst_files_128 = numpy.array(lst_files_128)
    print 'Reading 128 Data in LST Range:',len(numpy.array(good_files_128)),'files'
    t1,d1,f1 = capo.arp.get_dict_of_uv_data(good_files_128,antstr='1_4',polstr='xx')
    d1 = d1[(1, 4)]['xx']
    plt.subplot(2,2,1)
    capo.arp.waterfall(d1,mx=0,drng=4,mode='log',extent=(0,202,lst_files_128.max()*12/numpy.pi,lst_files_128.min()*12/numpy.pi))
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
        if lst > lstmin*numpy.pi/12 and lst < lstmax*numpy.pi/12:
            good_files_64.append(file)
            lst_files_64.append(lst)
    lst_files_64 = numpy.array(lst_files_64)
    print 'Reading 64 Data in LST Range:',len(numpy.array(good_files_64)),'files'
    t2,d2,f2 = capo.arp.get_dict_of_uv_data(good_files_64,antstr='1_4',polstr='I')
    d2 = d2[(1, 4)]['I'] #24_48
    plt.subplot(2,2,2)
    capo.arp.waterfall(d2,mx=2,drng=4,mode='log',extent=(0,202,lst_files_64.max()*12/numpy.pi,lst_files_64.min()*12/numpy.pi))
    plt.colorbar()
    plt.ylabel('LST Hours')
    plt.xlabel('Freq Chan')
    plt.title('64 Data (Absolute calibrated)')

    freq_chans = numpy.arange(115,135,5) #hard-coded frequency channel to get factor from
    #freq_chans = [125]
    factor_all = []
    for f,freq_chan in enumerate(freq_chans):
        d2_f = d2[:,freq_chan]
        d1_f = d1[:,freq_chan]
        factors = d2_f/d1_f
        factor = numpy.real(numpy.mean(factors))
        factor_all.append(factor)
    print 'YOUR FACTOR IS:',numpy.mean(factor_all)
    factor = 1000
    print 'Applying Sample Calibration...'
    d1_cal = numpy.zeros_like(d1)
    d1_cal = d1*factor
    plt.subplot(2,2,3)
    capo.arp.waterfall(d1_cal,mx=2,drng=4,mode='log',extent=(0,202,lstmax*12/numpy.pi,lstmin*12/numpy.pi))
    plt.colorbar()
    plt.ylabel('LST Hours')
    plt.xlabel('Freq Chan')
    plt.title('128 Data Calibrated to 64 Data')

if opts.plot == True:
    plt.tight_layout()
    plt.show()


"""

if opts.factor != None and opts.abscal == True and opts.poly == False: #if factor is given in command-line
    #factor = opts.factor
    factor = numpy.load(opts.factor)['bandpass']

# Absolute calibrate
if opts.factor == None:
    print 'Saving bandpass.npz'
    numpy.savez('bandpass.npz',bandpass=factor)

if opts.abscal == True:

    def mfunc(uv,p,d):
        d *= factor
        try:
            uv['var'] = uv['var']* factor**2
        except: pass
        return p,d

    for file in args:
        uvi = aipy.miriad.UV(file)
        newfile = file+'G'
        if os.path.exists(newfile): 
            print '   %s exists. Skipping...' % newfile
            continue
        uvo = aipy.miriad.UV(newfile, status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc=mfunc, append2hist='SIMPLE_ABSCAL:' + ' '.join(sys.argv) + '\n')
        print file, '->', newfile
    





