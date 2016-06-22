#! /usr/bin/env python
import numpy as n
import aipy as a
import sys, os, optparse 

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('--scale', type=float, 
            help='multiply the data by this number')
opts, args = o.parse_args(sys.argv[1:])
if opts.scale is None:
   print "ERROR: no scale supplied."
   sys.exit()


def mfunc(uv, p, d, f):
    uvw,t,(ij) = p
    d = d*opts.scale
    try:
        uv['var'] = uv['var']* opts.scale**2
    except(KeyError):
        pass
    
    #for stat in stats:
    #    if stat in uv.vartable.keys():uv[stat] = (uv[stat]* opts.scale).astype(n.double)
    return p, d, f
     

for infile in args:
    outfile = infile+'S'
    print infile, '-->', outfile
    if os.path.exists(outfile):
        print 'File exists, skipping...'
        continue
    
    uvi = a.miriad.UV(infile)
    uvo = a.miriad.UV(outfile, status='new')
    uvo.init_from_uv(uvi)
    stats = ['median','max','min']
    for stat in stats:
        uvo.add_var(stat,'d')
    for p,d,f in uvi.all(raw=True):
        d *= opts.scale
        uvo['var'] =uvi['var']*opts.scale**2
        uvo['cnt'] = uvi['cnt']
    #uvo.pipe(uvi, mfunc, raw=True, append2hist='SCALE_DATA: ' + ' '.join(sys.argv) + '\n')
        for stat in stats:
            if stat in uvi.vartable.keys(): uvo[stat] =uvi[stat]*opts.scale
        uvo.write(p,d,f)
    uvo['history'] += 'SCALE_DATA: ' + ' '.join(sys.argv) + '\n'
    del(uvo)
