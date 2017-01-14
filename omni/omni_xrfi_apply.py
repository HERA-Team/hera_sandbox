#! /usr/bin/env python
import aipy as a, capo, sys, optparse, numpy as np
import os

o = optparse.OptionParser()
o.add_option('--to_npz', action='store_true',
            help='Write out flags to npz file.')
o.add_option('--omnipath',dest='omnipath',default='%s.npz',type='string',
            help='Format string (e.g. "path/%s.npz", where you actually type the "%s") which converts the input file name to the omnical npz path/file.')
o.add_option('--per_ant',action='store_true',
             help='Calculate and apply flags for each antenna seperately.')
o.add_option('--boxside', type='int', default=8, 
              help='Size of the side of a box in flagger.')
o.add_option('--sig', type='float', default=6.,
             help='Number of sigma to flag above.')
o.add_option('--sigl', type='float', default=2.,
             help='Number of sigma for secondary cut in flagger.')
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
    m,_,_,_ = capo.omni.from_npz(omnifile)
#    if opts.per_ant:
#        for chi in m.keys() if chi.endswith(opts.pol[0]):
    f=capo.xrfi.omni_chisq_to_flags(m['chisq'], K=opts.boxside, sigma=opts.sig, sigl=opts.sigl)
    if opts.to_npz:
        flags = dict(zip(map(str,m['jds']),f))
        np.savez(filename+'R', **flags)
        continue

    times = []
    flags = dict(zip(m['jds'],f))
    def mfunc(uv,p,d,f):
        global times
        _,t,(a1,a2) = p
        p1,p2 = pol = a.miriad.pol2str[uv['pol']]
        if len(times) == 0 or times[-1] != t: times.append(t)
        return p, np.where(np.logical_or(f,flags[t]),0,d), np.logical_or(f,flags[t])
    
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(newfile, status='new')
    uvo.init_from_uv(uvi)
    print '    Saving', newfile
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist='OMNI_XRFI: ' + ' '.join(sys.argv) + '\n')
    
    
            
    
    
        
    
     
