#! /usr/bin/env python
'''Grab a set of visibilities at the same LST and then drop into an ipython
prompt, where all the data will be in a dictionary "data" with keys of LST
and values being a list of (baseline, pol, juldate, spectrum) tuples.'''
import aipy as a, numpy as n
import sys, optparse, ephem

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
o.add_option('-b', '--bin', dest='bin', type='float', default=0.001,
    help='Bin size in LST.  Default 0.001')
o.add_option('-l', '--lstrng', dest='lstrng', default='0_6.2832',
    help='Range of LSTs to bin.')
o.add_option('--tfile', dest='tfile', type='float', default=3600,
    help='Length of time spanned by each input file.  Helps in filtering out files that are not needed for the lst range being processed.')
opts, args = o.parse_args(sys.argv[1:])

lsts = n.arange(0, 2*n.pi, opts.bin)
# Select only some of these LST bins to analyze
lstrng = map(float, opts.lstrng.split('_'))
if lstrng[0] < lstrng[1]:
    lsts = lsts.compress(n.logical_and(lsts >= lstrng[0], lsts < lstrng[1]))
else:
    lsts = lsts.compress(n.logical_or(lsts >= lstrng[0], lsts < lstrng[1]))

tfile = 2 * n.pi * opts.tfile / (24. * 3600)
data = {}
for lst in lsts:
    data[lst] = {}
    for pol in opts.pol.split(','):
        data[lst][pol] = {}

print 'Filtering input files for LSTs of interest'
nargs = []
for f in args:
    uv = a.miriad.UV(f)
    start_t = uv['lst']
    end_t = (start_t + tfile) % (2*n.pi)
    if start_t < end_t:
        if lstrng[0] < lstrng[1]:
            if end_t < lstrng[0] or start_t > lstrng[1]: continue
        else: 
            if end_t < lstrng[0] and start_t > lstrng[1]: continue
    else:
        if lstrng[0] < lstrng[1]:
            if start_t > lstrng[1] and end_t < lstrng[0]: continue
        # Never bail if both wrap...
    nargs.append(f)

for f in nargs:
    uv = a.miriad.UV(f)
    print 'Reading', f
    sys.stdout.flush()
    if lstrng[0] < lstrng[1]:
        uv.select('lst', lstrng[0], lstrng[1])
    else:
        uv.select('lst', lstrng[1], lstrng[0], include=False)

    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    # Gather data from file
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        bl = a.miriad.ij2bl(i,j)
        pol = a.miriad.pol2str[uv['pol']]
        lst = n.round(uv['lst'] / opts.bin) * opts.bin
        # Only take this LST if we have a bin for it already allocated
        if not data.has_key(lst): continue
        try: data[lst][pol][bl].append((t,d,f))
        except(KeyError): data[lst][pol][bl] = [(t,d,f)]

cnt = n.array([sum([n.logical_not(x[2]).astype(n.int) for x in data[lst][pol][bl]]) \
    for lst in data \
    for pol in data[lst] \
    for bl in data[lst][pol]])
add = n.array([sum([x[1] for x in data[lst][pol][bl]]) \
    for lst in data \
    for pol in data[lst] \
    for bl in data[lst][pol]])
d0 = add/n.where(cnt==0, 1, cnt)

d = n.array([n.where(d0[i] == 0, 0, x[1] / d0[i]) \
    for i, lst in enumerate(data) \
    for pol in data[lst] \
    for bl in data[lst][pol] \
    for x in data[lst][pol][bl]])

d = n.ma.array(d, mask=n.where(d==0,1,0))
g = n.ma.median(d, axis=1)
g = n.ma.masked_where(n.abs(g.real - 1) > .5, g)
#g = g.real.clip(.5,1.5) + 1j * g.imag
g.shape = (g.size,1)
d /= g

dabs = n.abs(d)
dang = n.angle(d)

import pylab as p
p.subplot(231)
p.imshow(dabs, vmax=1.1, vmin=.9, aspect='auto')
p.title(r'$\left|V_{ij}/\langle V_{ij}\rangle\right| 1\pm 0.1$')
#p.colorbar(shrink=.5)
p.subplot(232)
p.imshow(dang, vmax=.2, vmin=-.2, aspect='auto')
p.title('Phase $\pm 0.2$')
#p.colorbar(shrink=.5)
p.subplot(233)
p.plot(g.real, g.imag, '.')
p.title('Gain')
p.xlabel('Real(Gain)')
p.ylabel('Imag(Gain)')

p.subplot(234)
p.plot(n.ma.average(dabs, axis=0), 'k.', label='avg')
p.plot(n.ma.std(dabs, axis=0), 'r.', label='std')
p.ylim(0,1.25)

p.subplot(235)
p.plot(n.average(dang, axis=0), 'k.', label='avg')
p.plot(n.std(dang, axis=0), 'r.', label='std')
p.legend()
p.ylim(-.25,.75)

#p.imshow(n.log10(n.abs(d)), vmax=0, vmin=-1, aspect='auto')
p.show()

#from IPython.Shell import IPShellEmbed; IPShellEmbed('')()
