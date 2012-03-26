#!/usr/bin/env python
#
#  stack_fits.py
#  
#
#  Created by Danny Jacobs on 11/14/09.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse, pyfits,ephem
"""
Stack input fits files into a cube. 
stack_fits.py *.fits  -o stack.fits

NB: All unitary dimensions are left out of the frames.
"""
o = optparse.OptionParser()
a.scripting.add_standard_options(o)
o.add_option('-o', dest='outfile', default='newstack.fits',
    help='Output file.')
o.add_option('--name',dest='axes3',default='auto',
    help='New axis label.  Examples: auto (default), freq, channel')
o.add_option('--origin',dest='freq',default=1,
    help='Value of axis at origin of new axis. default=1')
o.add_option('--inc',dest='d_freq',default=1,
    help='Axis increment. eg channel width, or sample rate. Default=1')
o.add_option('-v',dest='verb',action='store_true',default=False,
    help='Print more')
opts, args = o.parse_args(sys.argv[1:])

cubedata = []
freqs = []
if opts.verb: print ' '.join(sys.argv)
curpos = []
for slice in args:
    print "loading: "+slice,
    data, kwds = a.img.from_fits(slice)
    print ephem.hours(kwds['ra'] * a.img.deg2rad),
    print ephem.degrees(kwds['dec'] * a.img.deg2rad)
    if curpos ==[]: curpos = [kwds['ra'],kwds['dec']]
    elif curpos[0]!=kwds['ra'] or curpos[1]!=kwds['dec']:
        print "rejecting channel with wrong pointing!" 
        continue
    cubedata.append(data)
    try: freqs.append(kwds['freq'])
    except(KeyError): opts.axes3 = 'channel'
cubedata = [cubedata[i] for i in n.argsort(freqs)]

#Add blank images for any missing channels.
freqs = n.sort(n.array(freqs))
chans = n.diff(n.sort(n.array(freqs)))
blank = n.ones(data.shape[:2])*n.nan
blank.shape = blank.shape + 2 * (1,)
if opts.verb: print "found frequencies",freqs/1e6
if opts.verb: print "channel widths",chans/1e6#n.diff(n.sort(n.array(freqs))/1e6)
if opts.verb: print "missing channels",chans/n.min(chans)-1
cubedata = n.dstack(cubedata)

if opts.verb: print "Adding blank channels..."
if opts.verb: print cubedata.shape, freqs.shape
for df,i in reversed(zip(chans,range(len(freqs)))):
    n_missing_chans = df/n.min(chans)-1
    if n_missing_chans>0:
        for j in range(n_missing_chans):
            if opts.verb: print df,i,j,n_missing_chans
            if opts.verb: print "inserting a blank channel at %3.1f MHz" %(freqs[i+j]/1e6)            
            cubedata = n.dstack([cubedata[:,:,:i+j+1,:],blank,cubedata[:,:,i+j+1:]])
            freqs = n.insert(freqs,i+j+1,freqs[i+j]+n.min(chans)*(j+1))
if opts.verb:
    for i,(f,v) in enumerate(zip(freqs/1e6,cubedata[0,0,:,0])):            
        if i==0: print "%3.1f,%2e"%(f,v)
        else: print "%3.1f,%2e,%2.1f"%(f,v,freqs[i]-freqs[i-1])
if opts.verb: print cubedata.shape, freqs.shape            
if opts.verb: print "... done adding blank channels."

if opts.axes3=='auto':
    print "finding frequency axis automagically"
    opts.freq =  n.min(freqs)
    opts.d_freq = n.abs(n.diff(freqs)[0])
    opts.axes3 = 'FREQ'
    print "origin = %s MHZ, d_freq = %e MHz"%(opts.freq/10**6, opts.d_freq/10**6)

#a.img.to_fits(opts.outfile,cubedata,axes=kwds['axes'],ra=kwds['ra'],
#    dec=kwds['dec'],d_ra=kwds['d_ra'],d_dec=kwds['d_dec'],freq=0,d_freq=1)
files = [l.split('/')[-1] for l in args]
history = "stack_fits.py:"+ '\n  '+'\n  '.join(files)
if opts.verb: print history
kwds['CTYPE3']=opts.axes3
kwds['CRVAL3']=opts.freq
kwds['CDELT3']=opts.d_freq
kwds['BLANK']=-1
a.img.from_fits_to_fits(args[0],opts.outfile,cubedata,kwds,
    history=history)
#file = pyfits.open(opts.outfile,mode='update')
print 'created: '+opts.outfile
if opts.verb: print "end stack_fits.py"
