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
o.add_option('--omniruntag',dest='omniruntag',default='',type='string',
            help='Tag for omni run on npz file. Default is empty.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to .omni_output npz files. Include final / in path.')
opts,args = o.parse_args(sys.argv[1:])

### Save Options ###
pol = 'xx' #XXX shouldn't be hard-coded in the future
pol1 = pol[0]
pol2 = pol[1]

### mfunc ###
def mfunc(uv,p,d): #loops over time and baseline
    a1,a2 = p[2]
    #XXX To do: get xtalk for (a1,a2) and subtract it from d
        #if d is (a2,a1), then conjugate xtalk
        #if no xtalk, keep going without subtracting
    key1 = pol1+',gains,'+str(a1)
    key2 = pol1+',gains,'+str(a2)
    time_index = numpy.where(t_file == p[1])[0][0]
    try: g_ij = npz[key1][time_index]*npz[key2][time_index].conj()
    #XXX To do: split up gi and gj and if there's no solution, set gain to 1
        #dictionary.get(key,1) #returns 1 if no key
    except: g_ij = numpy.ones_like(d) #if no antenna solution, set gain to one
    g_ij = g_ij.conj() #Omnical conjugation is backwards
    d /= g_ij #calibration happens here
    return p,d

### Read Data and Solutions ###
for f in range(len(args)): #loop over files
    file = args[f]
    print str(f+1)+'/'+str(len(args))+': '+'Reading '+str(file)
    #tag = '.'.join(file.split('.')[:-1])
    tag = 'zen.'+'.'.join(file.split('.')[:-1][-3:])
    npz = numpy.load(opts.omnipath+tag+'.omni_output'+opts.omniruntag+'.npz')
    print '   calibrating'
    t_file,d_file,f_file = capo.arp.get_dict_of_uv_data([file],antstr='cross',polstr=pol)
    uvi = aipy.miriad.UV(file)
    newfile = tag+'.'+file.split('.')[-1]+'O'+opts.omniruntag
    uvo = aipy.miriad.UV(newfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc)
    print '   saving',newfile
