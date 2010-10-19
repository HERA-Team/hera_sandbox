#!/usr/bin/env python
#
#  plot_chan.py
#  
#
#  Created by Danny Jacobs on 6/18/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#
"""
Plot single channels on the same axis (not appended). Choose a single channel, 
or a range of channels to be summed.
"""

import aipy as a, numpy as n, pylab as p, math, sys, optparse

o = optparse.OptionParser()
o.set_usage('plot_chan.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, pol=True,ant=True,cal=True,dec=True)
o.add_option('-f','--freq',dest='freq',
             help="""Edge of a desired channel [GHz]""")
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
o.add_option('-t','--time_units',dest='ttype',default='jd',
    help="""Select time axis units: jd (default), lst or index. """)
o.add_option('-o', '--out_file', dest='out_file', default='',
    help='If provided, will save the figure to the specified file instead of popping up a window.')

opts, args = o.parse_args(sys.argv[1:])

#convert from frequency to channel number
#opts.bw = float(opts.bw)/10**6 #convert to GHz
uv = a.miriad.UV(args[0])
#if opts.bw<=uv['sdf']:opts.bw=uv['sdf']
chan = str(round((float(opts.freq)-
            float(uv['sfreq']))/float(uv['sdf']))).split('.')[0]
if int(chan)<0:chan=str(0)
elif int(chan)>uv['nchan']:chan=str(uv['nchan'])
#chan = str(str(round((float(opts.freq)-
#            float(uv['sfreq']))/float(uv['sdf']))).split('.')[0]+','+
#         str(round((float(opts.freq)+float(opts.bw)-
#                    float(uv['sfreq']))/float(uv['sdf']))).split('.')[0])
print "plotting channel # %s at %s GHz (bw=%e)" %(chan,opts.freq,uv['sdf'])
del(uv)
plot_t = {}#{'jd':[], 'lst':[], 'cnt':[]}
plot_x = {}
for file in args:
    print "loading file :",file
    sys.stdout.flush()
    uv = a.miriad.UV(file)
    aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
    chans = a.scripting.parse_chans(chan,uv['nchan']) 
    a.scripting.uv_selector(uv,opts.ant, opts.pol)
    uv.select('decimate',opts.decimate,opts.decphs)
#    if not plot_t.has_key(file): 
#        plot_t[file] = {'jd':[], 'lst':[], 'cnt':[]}
#        plot_x[file] = {}
    times = []
    for (uvw,t,(i,j)),d in uv.all():
        bl = '%d,%d' % (i,j)        
        if not plot_t.has_key(file):
            plot_t[file] = {}#{'jd':[], 'lst':[], 'cnt':[]}
            plot_x[file] = {}
        times.append(t)
        d = n.ma.compressed(d.take(chans))
        if opts.mode.startswith('phs'): d = n.angle(d)
        elif opts.mode.startswith('lin'): d = n.ma.absolute(d)
        elif opts.mode.startswith('real'): d = d.real
        elif opts.mode.startswith('imag'): d = d.imag
        elif opts.mode.startswith('log'):
            d = n.ma.absolute(d)
            d = n.ma.masked_less_equal(d, 0)
            d = n.ma.log10(d)
        else: raise ValueError('Unrecognized plot mode.')
        if d!=n.NaN:
            aa.date = a.phs.juldate2ephem(t)
            if plot_t[file].has_key(bl):
                plot_t[file][bl]['lst'].append(uv['lst'])
                plot_t[file][bl]['jd'].append(t)
                plot_t[file][bl]['cnt'].append(len(times)-1)
            else: plot_t[file][bl] = {'jd':[t],
                  'lst':[uv['lst']],'cnt':[len(times)-1]}
            if plot_x[file].has_key(bl): plot_x[file][bl].append(d)
            else: plot_x[file][bl] = [d]
files = plot_t.keys()
if not (len(files)>0): print "No data found, please check your inputs"; sys.exit()
else: bls = plot_x[files[0]].keys()
if len(bls)<1: print "No data found, this channel might be flagged."; sys.exit()
m2 = int(math.sqrt(len(bls)))
m1 = int(math.ceil(float(len(bls)) / m2))

p.figure()
print bls
for i,bl in enumerate(bls):
    for file in args:
        p.subplot(m1,m2,i+1)
        p.title(bl)
        p.plot(plot_t[file][bl][opts.ttype],plot_x[file][bl],label=file)
        p.legend()
print "finished plotting"
# Save to a file or pop up a window
if opts.out_file != '': p.savefig(opts.out_file)
else: p.show()

    
