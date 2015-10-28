#! /usr/bin/env python
import aipy as a
import numpy as n
from pylab import *
import optparse, sys, os
from capo.arp import get_dict_of_uv_data

o=optparse.OptionParser()
o.set_usage("print_blflagging.py [options]")
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True)
o.add_option('-v',action='store_true',help='turn on more verbs')
o.add_option('--detail',action='store_true',help='print stats for each bl')
opts,args=o.parse_args(sys.argv[1:])


t,dat,flg = get_dict_of_uv_data(args,opts.ant,opts.pol,verbose=opts.v)

bls = dat.keys()
nantot = 0
flagtot = 0
ntot = 0
if opts.detail:
    print 'bl pol nancount flagfrac'
for bl in bls:
    for pol in dat[bl].keys():
        if opts.detail:
            print a.miriad.bl2ij(bl),pol,n.sum(n.isnan(dat[bl][pol])),
            print n.sum(flg[bl][pol])/float(flg[bl][pol].size)
        nantot += n.sum(n.isnan(dat[bl][pol]))
        flagtot += n.sum(flg[bl][pol])
        ntot += n.float(flg[bl][pol].size)
print "summary"
print "nans = ",nantot
print "flagtot = ",flagtot,n.round(flagtot/ntot*100),'%'
print "total data = ",ntot


