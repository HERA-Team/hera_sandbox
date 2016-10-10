import aipy as a
import numpy as n
from pylab import *
import optparse, sys, os
from capo.arp import get_dict_of_uv_data
from scipy.interpolate import interp2d
from astropy.time import Time
o=optparse.OptionParser()
o.set_usage("ratio_two_files.py [options]")
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True,cal=True)
o.add_option('-v',action='store_true',help='turn on more verbs')
o.add_option('--pols',default='xx,xx',type=str,help='pols to compare, ex xx,I compares xx in file 1 to I in file 2 [default xx,xx]')
#o.add_option('--plot',dest='plot',default=False, action='store_true',\
#    help='Outputs plot to X-Window and saves plot')
opts,args=o.parse_args(sys.argv[1:])
#start by assuming things are lst aligned (they aren't)
filename = args[0]
pol_A = opts.pols.split(',')[0]
print "reading",filename
t_A,dat_A,flg_A = get_dict_of_uv_data([filename],opts.ant,pol_A,verbose=opts.v)
bls_A = dat_A.keys()

uv = a.miriad.UV(filename)
aa = a.cal.get_aa(opts.cal,uv['sdf'], uv['sfreq'], uv['nchan'])
freqs_A = aa.get_afreqs()
t_A = Time(t_A,scale='utc',format='jd',location=(aa.lat,aa.long))
lst_A = t_A.sidereal_time('apparent')
if len(bls_A)==0: 
    print "no data found in file ",filename[0]
    sys.exit()
D_A = n.ma.masked_where([flg_A[bl][pol_A] for bl in bls_A],[dat_A[bl][pol_A] for bl in bls_A])

filename = args[1]
print "reading",filename
pol_B = opts.pols.split(',')[1]
t_B,dat_B,flg_B = get_dict_of_uv_data([filename],opts.ant,pol_B,verbose=opts.v)
bls_B = dat_B.keys()
uv = a.miriad.UV(filename)
aa = a.cal.get_aa(opts.cal,uv['sdf'], uv['sfreq'], uv['nchan'])
freqs_B = aa.get_afreqs()
t_B = Time(t_B,scale='utc',format='jd',location=(aa.lat,aa.long))
lst_B = t_B.sidereal_time('apparent')
if len(bls_B)==0: 
    print "no data found in file ",filename[0]
    sys.exit()

D_B = n.ma.masked_where([flg_B[bl][pol_B] for bl in bls_B],[dat_B[bl][pol_B] for bl in bls_B])

for i,bl in enumerate(bls_A):
    try:
        subplot(131)
        A_amp = n.ma.masked_invalid(n.log10(n.abs(dat_A[bl][pol_A])))
        print A_amp.shape
        imshow(A_amp,aspect='auto')
        subplot(132)
        B_amp = n.ma.masked_invalid(n.log10(n.abs(dat_B[bl][pol_B])))
        print B_amp.shape
        imshow(B_amp,aspect='auto')
        subplot(133)
        #A_amp_model = interp2d(freqs_A,lst_A,A_amp)
        #A_amp_Bmatch = A_amp_model(freqs_B,lst_B)
        #A_amp_Bmatch.shape = (len(lst_B),len(freqs_B))
        imshow(A_amp-B_amp,aspect='auto') 
        colorbar()
        show()
    except(KeyError):continue
    
