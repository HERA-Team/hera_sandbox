#! /usr/bin/env python

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
o.set_usage('omini_apply.py [options] *uvcRRE')
o.set_description(__doc__)
#o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
#            help='Path and name of calfile.')
opts,args = o.parse_args(sys.argv[1:])

### Save Options ###
pol = 'xx' #XXX shouldn't be hard-coded in the future

### mfunc ###
def mfunc(uv,p,d): #loops over time and baseline
    a1,a2 = p[2]
    key1 = pol+',g_lin,'+str(a1)
    key2 = pol+',g_lin,'+str(a2)
    time_index = numpy.where(t_file == p[1])[0][0]
    try: g_ij = npz[key1][time_index]*npz[key2][time_index].conj()
    except: g_ij = numpy.zeros_like(d) #if no antenna solution, set gain to zero
    g_ij = g_ij.conj() #Omnical conjugation is backwards
    d /= g_ij #calibration happens here
    return p,d

### Read Data and Solutions ###
for f in range(len(args)): #loop over files
    file = args[f]
    print str(f+1)+'/'+str(len(args))+': '+'Reading '+str(file)
    tag = '.'.join(file.split('.')[:-1])
    npz = numpy.load(tag+'.omni_output.npz')
    print '   calibrating'
    t_file,d_file,f_file = capo.arp.get_dict_of_uv_data([file],antstr='cross',polstr=pol)
    uvi = aipy.miriad.UV(file)
    uvo = aipy.miriad.UV(file+'O',status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc)
