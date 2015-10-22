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

times = []

### mfunc ###
def mfunc(uv,p,d,f): #loops over time and baseline
    global times
    _,t,(a1,a2) = p
    p1,p2 = pol = aipy.miriad.pol2str[uv['pol']]
    if len(times) == 0 or times[-1] != t: times.append(t)
    if opts.xtalk:
        try: d -= xtalk[pol][(a1,a2)]
        except(KeyError):
            #print (a2,a1) in xtalk[pol] # Why should this happen XXX
            try: d -= xtalk[pol][(a2,a1)].conj()
            except(KeyError): pass
    time_index = len(times) - 1
    try: d /= gains[p1][a1][time_index].conj() # omnical conj is backwards
    except(KeyError): pass
    try: d /= gains[p2][a2][time_index] # omnical conj is backwards
    except(KeyError): pass
    return p,d,f

### Read Data and Solutions ###
for f,filename in enumerate(args):
    print 'Reading '+str(filename)
    #tag = '.'.join(filename.split('.')[:-1])
    tag = 'zen.'+'.'.join(filename.split('.')[:-1][-3:])
    newfile = tag+'.'+filename.split('.')[-1]+'O'+opts.omniruntag
    if os.path.exists(newfile):
        print newfile, 'exists.  Skipping...'
        continue
    npz = numpy.load(opts.omnipath+tag+'.omni_output'+opts.omniruntag+'.npz')
    gain_keys = [f for f in npz.files if f.find('gain') != -1]
    xtalk_keys = [f for f in npz.files if f.find('xtalk') != -1]
    gains, xtalk = {}, {}
    for k in gain_keys:
        pol,_,i = k.split(',')
        if not gains.has_key(pol): gains[pol] = {}
        gains[pol][int(i)] = npz[k]
    for k in xtalk_keys:
        pol,_,i,j = k.split(',')
        bl = (int(i[1:]),int(j[:-1]))
        if not xtalk.has_key(pol): xtalk[pol] = {}
        xtalk[pol][bl] = npz[k]
    print '   calibrating'
    if opts.xtalk:
        print '   and subtracting xtalk'
    #t_file,d_file,f_file = capo.arp.get_dict_of_uv_data([file],antstr='cross',polstr=pol)
    uvi = aipy.miriad.UV(filename)
    uvo = aipy.miriad.UV(newfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True)
    print 'Done.'
