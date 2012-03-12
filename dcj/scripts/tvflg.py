#!/usr/bin/env python
#
#  tvflg.py
#  
#
#  Created by Danny Jacobs on 12/7/08.
#  Creative Commons (cc 2008) Some rights reserved
#
import aipy as a, numpy as n, math, sys, pylab as p,optparse


o = optparse.OptionParser()
o.set_usage('tvflg.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True, cal=True,chan=True)
o.add_option('--docalib',dest='docalib',action='store_true',default=False,
        help='Use all available calibration information in the location file')

opts, args = o.parse_args(sys.argv[1:])


def convert_arg_range(arg):
    """Split apart command-line lists/ranges into a list of numbers."""
    arg = arg.split(',')
    return [map(float, option.split('_')) for option in arg]
def gen_chans(chanopt, uv, coords, is_delay):
    """Return an array of active channels and whether or not a range of
    channels is selected (as opposed to one or more individual channels)
    based on command-line arguments."""
    is_chan_range = True
    if chanopt == 'all': chans = n.arange(uv['nchan'])
    else:
        chanopt = convert_arg_range(chanopt)
        if coords != 'index':
            if is_delay:
                def conv(c):
                    return int(n.round(c * uv['sdf'] * uv['nchan'])) \
                        + uv['nchan']/2
            else:
                def conv(c): return int(n.round((c - uv['sfreq']) / uv['sdf']))
        else:
            if is_delay:
                def conv(c): return int(c) + uv['nchan']/2
            else:
                def conv(c): return c
        chanopt = [map(conv, c) for c in chanopt]
        if len(chanopt[0]) != 1:
            chanopt = [n.arange(x,y, dtype=n.int) for x,y in chanopt]
        else: is_chan_range = False
        chans = n.concatenate(chanopt)
    return chans.astype(n.int), is_chan_range


def length3(v):
    return n.sqrt(v[0]**2+v[1]**2+v[2]**2)

m = int(math.ceil(len(args)))


for uvcnt, uvfile in enumerate(args):
    uv = a.miriad.UV(uvfile)
    chans, is_chan_range = gen_chans(opts.chan, uv, 'index', False)
    freqs = n.arange(uv['sfreq'], uv['sfreq']+uv['nchan']*uv['sdf'], uv['sdf'])
    freqs = freqs.take(chans)
    if opts.docalib: aa = a.loc.get_aa(opts.loc, uv['sdf'], uv['sfreq'], uv['nchan'])
    print "loading file: "+uvfile
    #Choose the desired antennae
    a.scripting.uv_selector(uv, 'cross', 'yy')
    #Read the data from the UV file
    tv = {}
    times = []
    for (uvw,t,(i,j)),d in uv.all():
        bl = '%d,%d' % (i,j)
        #if opts.docalib: 1
        ind = length3(uvw)
        times.append(t)
#        if not tv.has_key(ind): tv[ind] =[]
        if not tv.has_key(bl): tv[bl] =[]
        if opts.docalib: d /= aa.passband(i,j)
        if not opts.chan is None: d = d.take(chans)
        else:d = d.sum(axis=0)
#        tv[ind].append(d)
        tv[bl].append(d)
    del(uv)
bls = tv.keys()
def sort_func(a, b):
     ai,aj = map(int, a.split(','))
     bi,bj = map(int, b.split(','))
     if bi > ai or (bi == ai and bj > aj): return -1
     return 1
bls.sort(cmp=sort_func)
#bls.sort()
if len(bls) == 0:
    print 'No data to plot.'
    sys.exit(0)
del(d)
d = n.vstack((tv[bls[0]],tv[bls[1]]))
for cnt,bl in enumerate(bls[2:-1]):
    print cnt,bl
    d = n.vstack((d,tv[bl]))
t1,t2 = times[0],times[-1]
bl1,bl2 = bls[0],bls[-1]
print 't1 =',t1,', t2 = ',t2, ', bl1 = ',bl1,', bl2 = ',bl2
#p.subplot(1,m,uvcnt+1)
#p.xticks(n.linspace(0,cnt,cnt),(bls))
if opts.docalib: p.title(uvfile+'+bandpass cal',horizontalalignment='left')
else: p.title(uvfile,horizontalalignment='left')
#p.imshow(n.angle(d.transpose()),extent=(bl1,bl2,t1,t2),aspect='auto',interpolation='nearest')
p.imshow(n.angle(d.transpose()),extent=(0,len(tv),t1,t2),aspect='auto',interpolation='nearest')
p.xticks(arange(0,len(tv),len(tv)/16.0),arange(0,16))
p.colorbar()
p.show()
