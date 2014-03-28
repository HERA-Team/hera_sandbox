#! /usr/bin/env python
"""
Takes derivatives of input uv files and spits out a spectrum (npz file).
Optionally plots.

usage: diff_uv.py -f -t -p *uv
"""

import aipy as a, numpy as n, sys, optparse, os 
from capo.pspec import jy2T
from pylab import *
o = optparse.OptionParser()
o.set_usage('plot_uv.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True,cal=True,pol=True)
o.add_option('-f', action='store_true',
    help='Do freq derivative.')
o.add_option('-t',action='store_true',
    help='Do time derivative.')
o.add_option('--plot',action='store_true',
    help='Plot.')
o.add_option('-i',action='store_true',
    help='clear saved npz and rerun on data')
o.add_option('-v',action='store_true',
    help='Verbose. print more debug stuff')
o.add_option('--bl',action='store_true',
    help='Difference across baselines!.')
opts, args = o.parse_args(sys.argv[1:])
if opts.cal != None:
    uv = a.miriad.UV(args[0])
    aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
    del(uv)
else: aa = None
uv = a.miriad.UV(args[0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan']).squeeze()
if opts.f:
    freqs = freqs[:-1]
if opts.v: print "print temp at f = ",freqs[len(freqs)/2],"GHz"
del(uv)
if not os.path.exists('diff_uv.npz') or opts.i:
    for uvfile in args:
        print 'Reading', uvfile
        uv = a.miriad.UV(uvfile)
        if not opts.cal is None:
            aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
        a.scripting.uv_selector(uv, opts.ant, opts.pol)
        curtime = 0
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            bl = '%d,%d,%d' % (i,j,uv['pol'])
            if opts.v: print t
            if opts.f:
                d = n.ma.diff(d)
            # if this a new time step, stop and add in the old data
            if curtime!=t and curtime!=0:
                if opts.bl:
                    D = n.sum([D[i]*(-1)**i for i in range(len(D))],axis=0)
                try:
                    dif += n.ma.array(D).squeeze() * S
                    w += n.logical_not(f).astype(n.float).squeeze()
                    if opts.v:print S
                except(NameError):
                    dif = n.ma.array(D).squeeze()
                    w = n.logical_not(f).astype(n.float).squeeze()
                    S = 1
                    if opts.v: print "RESET"
                if opts.t: S *= -1 #flop the sign every time
                if opts.v:
                    print "Trms = ",n.ma.mean(n.sqrt(dif*n.conj(dif))*jy2T(freqs)/w)
            #accumulate differen bls and pols
            if t != curtime:
                D = [d]
                curtime = t
            else:
                D.append(d) 
    dif = dif.squeeze()
    Trms = n.sqrt(n.real(dif*n.conj(dif))).squeeze()/w.squeeze()*jy2T(freqs)
    n.savez('diff_uv.npz',freqs=freqs,dif=dif.filled(0),weights=w,Trms=Trms.filled(0),mask=Trms.mask)
if opts.plot:
    #theoretical noise level
    F = n.load('diff_uv.npz')
    freqs = F['freqs']
    Trms = n.ma.masked_where(F['mask'],F['Trms'])
    print Trms.shape,freqs.shape
    plot(freqs*1e3,Trms,'k')
    savefig('diff_uv.png')
    grid()
    ylabel('mK')
    xlabel('MHz')
#    show()
