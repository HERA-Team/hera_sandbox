#! /usr/bin/env python
"""
Calibrate one file to another.
usage:

"""
import omnical
import aipy
import pylab
import numpy
import capo
import pickle
import optparse
import os, sys

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_match.py <options> <reference_data> <data_to_calibrate>"')
o.set_description(__doc__)
o.add_option('--calpar',dest='calpar',type='string',
            help='Path and name of .p calpar file.')
o.add_option('--redinfo',dest='redinfo',type='string',
            help='Path and name of .bin redundant info file.')
o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
            help='Path and name of calfile.')
#o.add_option('--xtalk',dest='xtalk',default=False,action="store_true",
#            help='Option to use xtalk command when performing lincal. Default is False.')
o.add_option('--omniruntag',dest='omniruntag',default='',type='string',
            help='Tag for omni run, if wanted. Default is empty.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .omni_output npz files. Include final / in path.')
#o.add_option('--ubls',dest='ubls',default=None,
#            help='Unique baselines to include. Ex: [(64,49),(64,10)]')
#o.add_option('--ex_ubls',dest='ex_ubls',default=None,
#            help='Unique baselines to exclude. Ex: [(64,49),(64,10)]')
#o.add_option('--ants',dest='ants',default=None,
#            help='Antennas to include. Ex: [64,49,10]')
#o.add_option('--ex_ants',dest='ex_ants',default=None,
#            help='Antennas to exclude. Ex: [64,49,10]')
opts,args = o.parse_args(sys.argv[1:])


### Save Options ###
CALFILE = opts.cal
CALPAR = opts.calpar
RED_INFO = opts.redinfo


### Read Firstcal Info ###

print 'Reading calpar and redinfo files...'

info = omnical.info.RedundantInfoLegacy() #reading old firstcal files
info.fromfile(RED_INFO)
reds = info.get_reds()
print reds
sys.exit()
d_calpar = pickle.load(open(CALPAR,'rb'))
gains = {} #dictionary indexed by pol
for k in d_calpar.keys():
    k_1 = k[0] #XXX just taking first letter (works on 'xx' or 'yy' only)
    gains[k_1] = {} 
    for i in xrange(d_calpar[k].shape[1]):
        gains[k_1][i] = d_calpar[k][:,i]

### calibrate file number 2 against file number 1
reference_file = args[0] #the file to calibrate against
subject_file = args[1]   #the file which we are going to calibrate

rt,rd,rf = capo.arp.get_dict_of_uv_data([reference_file],antstr='cross')
st,sd,sf = capo.arp.get_dict_of_uv_data([subject_file],antstr='cross')

tag = opts.omnipath+subject_file+'c'+reference_file

print "calibrating ",subject_file,"against ",reference_file

m2,g2,v2 = omnical.calib.redcal(data,info,xtalk=xtalk,gains=g,vis=v,uselogcal=False,removedegen=True)





        
