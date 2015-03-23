#! /usr/bin/env python
import numpy as n
import aipy as a
import pylab as p
import optparse,sys,capo,glob

POL = 'xx'
aa32 = a.cal.get_aa('psa898_v003', n.linspace(.1,.2,203))
aa64 = a.cal.get_aa('psa6240_v003', n.linspace(.1,.2,203))

ant32 = '0_16,1_17,2_18,3_19,8_16,9_17,10_18,11_19,8_24,9_25,10_26,11_27,4_24,5_25,6_26,7_27,4_20,5_21,6_22,7_23,12_20,13_21,14_22,15_23,12_28,13_29,14_30,15_31'
ant64 = '41_49,3_10,9_58,22_61,41_47,3_25,1_58,35_61,19_47,25_48,1_4,18_35,19_29,24_48,4_17,5_18,28_29,24_55,13_17,5_32,28_34,27_55,13_56,30_32,34_51,27_57,56_59,23_30,20_63,2_43'
#ant64 = '20_63,2_43,21_53,31_45,42_63,2_33,15_21,8_45,37_42,6_33,15_16,8_11,37_40,6_52,16_62,11_36,14_40,7_52,44_62,36_60,14_54,7_12,0_44,39_60,50_54,12_38,0_26,39_46'
#ant64 = '41_49,3_10,9_58,22_61,41_47,3_25,1_58,35_61,19_47,25_48,1_4,18_35,19_29,24_48,4_17,5_18,28_29,24_55,13_17,5_32,28_34,27_55,13_56,30_32,34_51,27_57,56_59,23_30,20_63,2_43,21_53,31_45,42_63,2_33,15_21,8_45,37_42,6_33,15_16,8_11,37_40,6_52,16_62,11_36,14_40,7_52,44_62,36_60,14_54,7_12,0_44,39_60,50_54,12_38,0_26,39_46'

time32, d32, f32 = capo.arp.get_dict_of_uv_data(['/data4/paper/arp/psa903/zen.2455928.38051.uvcRREzC'], ant32, POL, verbose=True)
time64, d64, f64 = capo.arp.get_dict_of_uv_data(['/data2/home/zakiali/psa_live/psa6247/zen.2456247.52186.uvcRREcAC'], ant64, POL, verbose=True)
time64_O, d64O, f64O = capo.arp.get_dict_of_uv_data(['/data2/home/zakiali/psa_live/forlstbinning_omnical_2/zen.2456247.52186.uvcRREcACO'], ant64, POL, verbose=True)


#get lsts of the two datasets.
lsts32 = []
lsts64 = []
lsts64_O = []
for jd in time32:
    aa32.set_jultime(jd)
    lsts32.append(aa32.sidereal_time())
for jd in time64:
    aa64.set_jultime(jd)
    lsts64.append(aa64.sidereal_time())
for jd in time64_O:
    aa64.set_jultime(jd)
    lsts64_O.append(aa64.sidereal_time())

lsts64_O = n.array(lsts64_O)
lsts64 = n.array(lsts64)
lsts32 = n.array(lsts32)
lst32mask = n.logical_and(lsts32>=lsts64[0], lsts32<=lsts64[-1]) 
lst64mask = lsts64>0; lst64mask[-1] = False #has extra lst at end
lst64_Omask = lsts64_O>0; lst64_Omask[-1] = False #has extra lst at end
freqs = n.linspace(.1,.2,203)
print lsts64[lst64mask]
print lsts64_O[lst64_Omask]
print lsts32[lst32mask]

blkey32 = d32.keys()[0]
total32 = n.zeros_like(d32[blkey32][POL][lst32mask,:])
blkey64 = d64.keys()[0]
total64 = n.zeros_like(d64[blkey64][POL][lst64mask,:])
total64_O = n.zeros_like(d64[blkey64][POL][lst64_Omask,:])
cnt32 = 0
cnt64 = 0
print total32.shape
print total64.shape
print total64_O.shape

for bl in d32.keys():
    #32 data
#    total32 += d32[bl][POL][lst32mask,:] * (364./424) #to correct for flux scale of 32 data. Correction is exact for 160MHz.
    total32 += d32[bl][POL][lst32mask,:] * (363.6/424) * (freqs/.160)**(-.76+.95) #to correct for flux scale of 32 data. Correction is exact for 160MHz.
    cnt32 += 1
for bl in d64.keys():
    #64 data
    total64 += d64[bl][POL][lst64mask,:]
    total64_O += d64O[bl][POL][lst64_Omask,:]
    cnt64 += 1

total32 /= cnt32
total64 /= cnt64
total64_O /= cnt64

#rms over time.
total32_spec = n.sqrt(n.mean(n.abs(total32)**2, axis=0))
total64_spec = n.sqrt(n.mean(n.abs(total64)**2, axis=0))
total64_Ospec = n.sqrt(n.mean(n.abs(total64_O)**2, axis=0))

ratio32_64 = total32_spec/total64_spec
ratio64_64O = total64_spec/total64_Ospec

chanrange = range(110,150)
print '32 / 64 data:        '
print '\tFrequency range for stats: ', freqs[chanrange][0], freqs[chanrange][-1]
print '\tMean: ',n.mean(ratio32_64[chanrange])
print '\tMedian: ',n.median(ratio32_64[chanrange])
print '\tVariance: ',n.std(ratio32_64[chanrange])**2
print '64 / 64_omnical data:        '
print '\tFrequency range for stats: ', freqs[chanrange][0], freqs[chanrange][-1]
print '\tMean: ',n.mean(ratio64_64O[chanrange])
print '\tMedian: ',n.median(ratio64_64O[chanrange])
print '\tVariance: ',n.std(ratio64_64O[chanrange])**2

p.semilogy(freqs, total32_spec, label='32')
p.semilogy(freqs, total64_spec, label='64')
p.semilogy(freqs, total64_Ospec, label='64O')

p.semilogy(freqs, ratio32_64, label='32/64')
p.semilogy(freqs, ratio64_64O, label='64/64_O')
p.semilogy(freqs, 382*(freqs/.15)**-.76, label='Pic Spec')
p.ylabel('Jy')
p.xlabel('GHz')
p.ylim(.3,5000)
p.legend(ncol=6, fontsize=10)

#p.plot(n.abs(dtotal[:,channel]), color='b')
#p.plot(n.abs(ddtotal[:,channel]), color='r')
#p.figure(2)
#p.plot(n.abs(dtotal[:,channel]/ddtotal[:,channel]))
p.show()
