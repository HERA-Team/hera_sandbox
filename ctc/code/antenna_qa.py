#!/usr/bin/env python

import aipy, capo as C, optparse, numpy as np, matplotlib.pyplot as plt, sys
from collections import Counter

o = optparse.OptionParser()
o.set_description(__doc__)
o.set_usage('antenna_qa.py *xx.uvcRRE (only use XX pol, script will find the others).')
aipy.scripting.add_standard_options(o, cal=True)
o.add_option('-A','--startant',dest='startant',default='0',help='antenna to start relative to (default is all antennae, i.e. start with zero)')
o.add_option('-b','--badants',dest='ba',default=None,help='bad antennae to remove, separated by commas. e.g. "-b 1,2,3"')
o.add_option('-v','--verbose',dest='verb',action='store_true',help='Toggle verbosity')
opts, args = o.parse_args(sys.argv[1:])

#parse files: assumes all pols are in the same directory
if 'xx' in args[0]:
    xxfiles = args
    xyfiles = [w.replace('xx','xy') for w in args]
    yxfiles = [w.replace('xx','yx') for w in args]
    yyfiles = [w.replace('xx','yy') for w in args]
else:
    print 'Hand this script the xx files. Assumes all pols are in the same directory.'
    raise Exception

#parse array data
if not opts.ba is None: badants = map(int,opts.ba.split(','))
else: badants = []
print 'Reading %s'%opts.cal
exec("import {calfile} as cal".format(calfile=opts.cal))
antpos = cal.prms['antpos']
nants = len(antpos.keys())
aa = aipy.cal.get_aa(opts.cal, np.array([0.15]))
sep2ij, blconj, bl2sep = C.zsa.grid2ij(aa.ant_layout)

#read data
tpxx,dpxx,fpxx = C.arp.get_dict_of_uv_data(xxfiles,antstr='cross',polstr='xx',verbose=True)
tpxy,dpxy,fpxy = C.arp.get_dict_of_uv_data(xyfiles,antstr='cross',polstr='xy',verbose=True)
tpyx,dpyx,fpyx = C.arp.get_dict_of_uv_data(yxfiles,antstr='cross',polstr='yx',verbose=True)
tpyy,dpyy,fpyy = C.arp.get_dict_of_uv_data(yyfiles,antstr='cross',polstr='yy',verbose=True)

nants = 112 #XXX

### Bad Antennas ###
antpowers = []
print 'Looping over reference antennas...'
for anchor_ant in range(int(opts.startant),nants):
    if anchor_ant in badants:
        avg_ratios.append(np.nan)
        continue
    #data analysis
    avg_ratios = []
    for ant in range(nants):
        if anchor_ant == ant: #neglect auto
            avg_ratios.append(np.nan)
            continue
        elif anchor_ant < ant: tup = (anchor_ant,ant)
        else: tup = (ant,anchor_ant)
        #if str(tup[0])+'_'+str(tup[1]) not in sep2ij['0,2'].split(','): #only use sep 0,2 to find bad ants
        #    continue
        hor_index1 = 86
        hor_index2 = 115
        avg_xx = np.nanmean(np.absolute(np.fft.fftshift(C.arp.clean_transform(dpxx[tup]['xx'],fpxx[tup]['xx']),axes=1))[:,hor_index1:hor_index2]) / np.nanmean(np.delete(np.absolute(np.fft.fftshift(C.arp.clean_transform(dpxx[tup]['xx'],fpxx[tup]['xx']),axes=1)),np.s_[hor_index1:hor_index2],axis=1)) #dynamic range
        avg_xy = np.nanmean(np.absolute(np.fft.fftshift(C.arp.clean_transform(dpxy[tup]['xy'],fpxy[tup]['xy']),axes=1))[:,hor_index1:hor_index2]) / np.nanmean(np.delete(np.absolute(np.fft.fftshift(C.arp.clean_transform(dpxy[tup]['xy'],fpxy[tup]['xy']),axes=1)),np.s_[hor_index1:hor_index2],axis=1)) #dynamic range
        avg_yx = np.nanmean(np.absolute(np.fft.fftshift(C.arp.clean_transform(dpyx[tup]['yx'],fpyx[tup]['yx']),axes=1))[:,hor_index1:hor_index2]) / np.nanmean(np.delete(np.absolute(np.fft.fftshift(C.arp.clean_transform(dpyx[tup]['yx'],fpyx[tup]['yx']),axes=1)),np.s_[hor_index1:hor_index2],axis=1)) #dynamic range
        avg_yy = np.nanmean(np.absolute(np.fft.fftshift(C.arp.clean_transform(dpyy[tup]['yy'],fpyy[tup]['yy']),axes=1))[:,hor_index1:hor_index2]) / np.nanmean(np.delete(np.absolute(np.fft.fftshift(C.arp.clean_transform(dpyy[tup]['yy'],fpyy[tup]['yy']),axes=1)),np.s_[hor_index1:hor_index2],axis=1)) #dynamic range
        ratio = (avg_xy + avg_yx)/(avg_xx + avg_yy) #identify bad antennas
        avg_ratios.append(ratio)    
    baddies = np.where(avg_ratios > np.nanmean(avg_ratios)+2*np.nanstd(avg_ratios))[0]
    if opts.verb: print anchor_ant,':',baddies
    antpowers.append(np.nanmean(avg_ratios))

print "   FINAL BAD ANTENNAS:" 
print np.where(antpowers > np.mean(antpowers)+2*np.std(antpowers))[0]
import IPython;IPython.embed()
