#! /usr/bin/env python

import aipy as a, sys, optparse, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
o.add_option('-t', '--time', dest='time', default=1000000, type='int',
    help='Number of integrations to crop file too.  Default is a huge number.') 
opts,args = o.parse_args(sys.argv[1:])

#cntr = 0	
global cntr

def mfunc(uv, p, d,f):
    global cntr
    uvw,t,(i,j) = p
    cntr = cntr + 1
    #print cntr
    if cntr > (opts.time + 1): return None, None, None
    else: return p, d, f
 
for filename in args:
    print filename, '->', filename+'T'
    uvi = a.miriad.UV(filename)
    if os.path.exists(filename+'T'):
        print '    File exists... skipping.'
        continue
    uvo = a.miriad.UV(filename+'T',status='new')
    uvo.init_from_uv(uvi)
    cntr = 1
    uvo.pipe(uvi,mfunc=mfunc,raw=True,append2hist='Reduced number of integrations to %d \n' % opts.time)
