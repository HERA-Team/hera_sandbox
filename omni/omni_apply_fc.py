#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_apply.py [options] *uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,pol=True)
o.add_option('--xtalk',dest='xtalk',default=False,action='store_true',
            help='Toggle: apply xtalk solutions to data. Default=False')
o.add_option('--omnipath',dest='omnipath',default='%s.npz',type='string',
            help='Format string (e.g. "path/%s.npz", where you actually type the "%s") which converts the input file name to the omnical npz path/file.')
opts,args = o.parse_args(sys.argv[1:])


#File Dictionary
pols = opts.pol.split(',')
files = {}
for filename in args:
    files[filename] = {}
    for p in pols:
        fn = filename.split('.')
        fn[3] = p
        files[filename][p] = '.'.join(fn)

firstfile = numpy.load('fcgains.xx.npz')
gains = {'x':{}}
for ant in firstfile.keys():
    if ant.isdigit():
        gains['x'][int(ant)] = firstfile[ant]

### Read Data and Solutions ###
for f,filename in enumerate(args):
        newfile = files[filename][p].split('/')[-1]+'F' #saves in cwd
#        omnifile = opts.omnipath % '.'.join(filename.split('/')[-1].split('.')[0:3])
        if os.path.exists(newfile):
            print '    %s exists.  Skipping...' % newfile
            continue
        times = []
    
        def mfunc(uv,p,d,f): #loops over time and baseline
            global times #global list
            _,t,(a1,a2) = p
            p1,p2 = pol = aipy.miriad.pol2str[uv['pol']]
            if len(times) == 0 or times[-1] != t: times.append(t) #fill times list
            if opts.xtalk: #subtract xtalk
                try: d -= xtalk[pol][(a1,a2)]
                except(KeyError):
                    try: d -= xtalk[pol][(a2,a1)].conj()
                    except(KeyError): pass
            ti = len(times) - 1 #time index
            try: d*=gains[p1][a1]
            except(KeyError): pass
            try: d*=gains[p1][a2].conj()
            except(KeyError): pass
            #try: d /= gains[p1][a1][ti] #apply gains
            #except(KeyError): pass
            #try: d /= gains[p2][a2][ti].conj() 
            #except(KeyError): pass
            return p, numpy.where(f,0,d), f
    
        if opts.xtalk: print '    Calibrating and subtracting xtalk'
        else: print '    Calibrating'
        uvi = aipy.miriad.UV(files[filename][p])
        uvo = aipy.miriad.UV(newfile,status='new')
        uvo.init_from_uv(uvi)
        print '    Saving', newfile
        uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='OMNICAL: ' + ' '.join(sys.argv) + '\n')
        
