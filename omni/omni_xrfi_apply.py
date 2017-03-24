#! /usr/bin/env python
import aipy as a, capo, sys, optparse, numpy as np
import os

o = optparse.OptionParser()
o.add_option('--omnipath',dest='omnipath',default='%s.npz',type='string',
            help='Format string (e.g. "path/%s.npz", where you actually type the "%s") which converts the input file name to the omnical npz path/file.')
opts,args = o.parse_args(sys.argv[1:])

#File Dictionary 
flags = {}
for f,filename in enumerate(args):
    newfile = filename+'R'
    if os.path.exists(newfile): 
        print newfile, 'exists, skipping...'
        continue
    omnifile = opts.omnipath % '.'.join(filename.split('/')[-1].split('.')[0:4])
    print '    Omnical npz:', omnifile
    flags = np.load(omnifile) 
    times = []

    def mfunc(uv,p,d,f):
        global times
        _,t,(a1,a2) = p
        p1,p2 = pol = a.miriad.pol2str[uv['pol']]
        if len(times) == 0 or times[-1] != t: times.append(t)
        return p, np.where(np.logical_or(f,flags[str(t)]),0,d), np.logical_or(f,flags[str(t)])
    
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(newfile, status='new')
    uvo.init_from_uv(uvi)
    print '    Saving', newfile
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='OMNI_XRFI: ' + ' '.join(sys.argv) + '\n')
    
    
            
    
    
        
    
     
