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
o.set_usage('omni_apply.py [options] *uvcRRE')
o.set_description(__doc__)
#o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
#            help='Path and name of calfile.')
o.add_option('--omniruntag',dest='omniruntag',default='',type='string',
            help='Tag for omni run on npz file. Default is empty.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to .omni_output npz files. Include final / in path.')
o.add_option('--xtalk',dest='xtalk',default=False,action='store_true',
            help='Apply xtalk solutions to data.')
opts,args = o.parse_args(sys.argv[1:])

### Save Options ###
pol = 'xx' #XXX shouldn't be hard-coded in the future
pol1 = pol[0]
pol2 = pol[1]

### mfunc ###
def mfunc(uv,p,d): #loops over time and baseline
    a1,a2 = p[2] #data always has bl of (lower num, higher num)
    if opts.xtalk:
        try: #xtalk isn't always like that
            xtalk = npz[pol+',xtalk,('+str(a1)+', '+str(a2)+')'] 
        except:
            try: 
                xtalk = npz[pol+',xtalk,('+str(a2)+', '+str(a1)+')']
                xtalk = xtalk.conj() #conjugate if bl is backwards
            except:
                xtalk = numpy.zeros_like(d) #some bls don't have xtalk soln
    else:
        xtalk = numpy.zeros_like(d)
    d = d-xtalk #subtract xtalk
    key1 = pol1+',gains,'+str(a1)
    key2 = pol1+',gains,'+str(a2)
    time_index = numpy.where(t_file == p[1])[0][0]
    try: g_i = npz[key1][time_index]
    except: g_i = numpy.ones_like(d)
    try: g_j = npz[key2][time_index]
    except: g_j = numpy.ones_like(d)
    g_ij = g_i*g_j.conj()

    #try: g_ij = npz[key1][time_index]*npz[key2][time_index].conj()
    #except: g_ij = numpy.ones_like(d) #if no antenna solution, set gain to one
    
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
    if opts.xtalk:
        print '   and subtracting xtalk'
    t_file,d_file,f_file = capo.arp.get_dict_of_uv_data([file],antstr='cross',polstr=pol)
    uvi = aipy.miriad.UV(file)
    newfile = tag+'.'+file.split('.')[-1]+'O'+opts.omniruntag
    uvo = aipy.miriad.UV(newfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc)
    print '   saving',newfile
