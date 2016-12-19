#! /usr/bin/env python

import aipy as a, sys, optparse, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, chan=True)
opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    outfile = os.path.basename(filename)+'A'
    print filename, '->', outfile
    uvi = a.miriad.UV(filename)
    if os.path.exists(outfile):
        print '    File exists... skipping.'
        continue
    chans = a.scripting.parse_chans(opts.chan,uvi['nchan'])
    print "extracting {nchan} channels".format(nchan=len(chans))
    def mfunc(uv,preamble,data):
        return preamble,data[chans]
        
        
    uvo = a.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi,mfunc=mfunc,append2hist='PULL CHANS:'+' '.join(sys.argv)+'\n')
