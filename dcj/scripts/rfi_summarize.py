#!/usr/bin/env python
#
#  rfi_summarize.py
#  
#
#  Created by Danny Jacobs on 10/4/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,ant=True,pol=True,chan=True,dec=True)
o.add_option('-r',default='R',
   help='Appendage for rfi. eg <filename>+R')
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
opts, args = o.parse_args(sys.argv[1:])

#difference two files to get flagged values A-B

D = []
for i,A in enumerate(args):
    uvA = a.miriad.UV(A)
    uvB = a.miriad.UV(A+opts.r)
    aa = a.cal.get_aa(opts.cal,uvA['sdf'],uvA['sfreq'],uvA['nchan'])
    print uvA['sfreq'],uvA['sdf'],uvA['nchan']
    chans = a.scripting.parse_chans(opts.chan, uvA['nchan'])
    aa.select_chans(chans)
    if i==0:
        dense = n.zeros_like(aa.get_afreqs())
        freqs = aa.get_afreqs()
    a.scripting.uv_selector(uvA, opts.ant, opts.pol)
    uvA.select('decimate', opts.decimate, opts.decphs)
    a.scripting.uv_selector(uvB, opts.ant, opts.pol)
    uvB.select('decimate', opts.decimate, opts.decphs)
#    afreqs = aa[0].beam.afreqs
    for (uvw,t,(i,j)),dA in uvA.all():
        pB,dB = uvB.read()
        d = (dB - dA)
#        d = d.take(chans)/(aa[i].amp*aa[j].amp)
        d = d.take(chans)/aa.passband(i,j)
        dense += n.logical_not(d.mask)
#        print n.abs(dB-dA),n.sum(d)
#        if (i,j)==(3,27) and n.sum(d)>0:
#            print i,j,t,len(n.argwhere(d==0)),n.max(n.abs(d[n.argwhere(d!=0)])),n.median(aa.passband(i,j))
#            p.plot(aa.passband(i,j))
#            p.show()
#            sys.exit()
#        d = n.abs(d.take(chans))
#        if n.ma.sum(d)==0:continue
        for c in n.argwhere(d>0):
            D.append(n.abs(d[c]))
#        D.append(d[n.argwhere(d>0)].squeeze())
#print D
#D = n.concatenate(D)
n.save(args[0]+'_flag_spec',dense)
print n.max(D),n.min(D)
dn,bins = n.histogram(n.log10(D))
print n.sum(dn),len(D)
s = bins[1:] - (bins[1]-bins[0])/2
#print dn.shape,bins.shape,s.shape
p.loglog(10**s,dn/((10**s)*n.diff(bins)),'.')
p.xlabel('flagged RFI flux [Jys]')
p.ylabel('$dn/ds [Jys^{-1}]$')


p.show()
        
        
        
    

#zip up the masked flux, time,frequency
