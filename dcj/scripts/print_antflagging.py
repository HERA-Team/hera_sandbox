#! /usr/bin/env python
import aipy as a
import numpy as n
from pylab import *
import optparse, sys, os
from capo.arp import get_dict_of_uv_data
from capo.dcj import file2jd

o=optparse.OptionParser()
o.set_usage("print_blflagging.py [options]")
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True)
o.add_option('-v',action='store_true',help='turn on more verbs')
o.add_option('--detail',action='store_true',help='print stats for each bl')
o.add_option('--awesome',action='store_true',help='try to seperate the weights using linear algebra instead of just \
averaging over baselines')
opts,args=o.parse_args(sys.argv[1:])
uv = a.miriad.UV(args[0])
nant=uv['nants']

jds = []
for filename in args:
    F = n.zeros((nant,nant))
    C = n.zeros_like(F)
    jds.append(file2jd(filename))
    t,dat,flg = get_dict_of_uv_data([filename],opts.ant,opts.pol,verbose=opts.v)
    bls = dat.keys()
    for bl in bls:
        (i,j) = a.miriad.bl2ij(bl)
        pols = dat[bl].keys()
        for pol in pols:
            F[i,j] += n.sum(flg[bl][pol])
            F[j,i] += n.sum(flg[bl][pol])
            F[i,i] += n.sum(flg[bl][pol])
            F[j,j] += n.sum(flg[bl][pol])
            Size = flg[bl][pol].size
            C[i,j] +=1
            C[j,i] +=1 
            C[i,i] +=1
            C[j,j] +=1
    #imshow(n.log(F))
    #show()
    F[C>0] /= C[C>0]
    U,S,V = n.linalg.svd(F)
    per_ant_flagging = n.sqrt(n.dot(F,V[0])**2/S[0])
    print jds[-1],
    for i,antflag in enumerate(per_ant_flagging):
        if opts.awesome:
            print n.round(antflag**2/float(Size),2),
        else:
            print n.round(F[i,i]/float(Size),2),
    print
sys.exit()
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
