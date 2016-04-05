#! /usr/bin/env python

"""

NAME:
    combine_times.py.
PURPOSE:
    Combines a bunch of UV files of continuous time chunks into one UV file
EXAMPLE CALL:
    ./combine_times.py --uvnew <path>
AUTHOR:
    Carina Cheng

"""

import aipy
import numpy
import optparse
import os, sys

o = optparse.OptionParser()
o.set_usage('combine_times.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--uvnew', default='/Users/carinacheng/capo/ctc/tables/test.uv',
             help='Path and name of outputted UV file.')
opts,args = o.parse_args(sys.argv[1:])

uvo = aipy.miriad.UV(opts.uvnew, status='new')
#uvo.init_from_uv(aipy.miriad.UV(args[0]))

for uvfile in args:

    uvi = aipy.miriad.UV(uvfile)
    print uvfile,'->',opts.uvnew

    uvo.init_from_uv(uvi)
    uvo.pipe(uvi)

    #for p,d,f in uvi.all(raw=True):   #XXX: this way doesn't update UV variables each time unless explicitly told
        #print uvi['lst']
        #uvo.copyvr(uvi) 
        #uvo['lst']=uvi['lst']
        #uvo.write(p,d,f)
        #print uvo['lst']
        
    
del(uvo)
