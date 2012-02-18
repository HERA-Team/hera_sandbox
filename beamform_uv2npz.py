#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys, re

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True, dec=True, 
    chan=True)
o.add_option('--thresh',default=0)
o.add_option('--cat', dest='cat', default='helm,misc',
    help='A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aipy._src.  Some catalogs may require a separate data file to be downloaded and installed.')
o.add_option('--altmin', dest='altmin', type='float', default=0,
    help="Minimum allowed altitude for pointing, in degrees.  When phase center is lower than this altitude, data is omitted.  Default is 0.")
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
if opts.chan is None: opts.chan = 'all'
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa.select_chans(chans)
inputs = re.search('(?<=BEAMFORM: ).*',open(args[0]+'/history').read()).group(0).split()

for i in inputs:
    i = i.split('=')
    if i[0]=='src': srcname=i[1]

print 'Looking at source ',srcname
srclist,cutoff,catalogs, = a.scripting.parse_srcs(srcname, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)
src = cat[srcname]
src.compute(aa)
print src

del(uv)

# Need to stack: spec, times, mask, x, y, z, path name
# File name is srcname__*

spec = []
times = []
mask = []
x = []
y = []
z = []

#spec, swgt = 0, 0
for filename in args:

    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    uv.select('decimate', opts.decimate, opts.decphs)
    curtime = None
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if curtime != t:
            curtime = t
            aa.set_jultime(t)
            src.compute(aa)
            xi,yi,zi = src.get_crds('top')
        if opts.altmin*n.pi/180 > src.alt: continue
        x.append(xi)
        y.append(yi)
        z.append(zi)
        times.append(t)
        if n.mean(f)>opts.thresh:
            mask.append(0)
        else: mask.append(1)
        spec.append(d)

x,y,z = n.array(x), n.array(y), n.array(z)
spec = n.array(spec)
mask = n.array(mask)
times = n.array(times)

#        d,f = d.take(chans), f.take(chans)
#        wgt = aa.bm_response(i,j,pol=opts.pol).squeeze()
#        d = n.where(f, 0, d); wgt = n.where(f, 0, wgt)
#        # Optimal SNR: down-weight beam-attenuated data 
#        # by another factor of the beam response.
#        d *= wgt; wgt *= wgt
#        spec += d; swgt += wgt

#spec = spec.real / swgt
#spec = n.abs(spec) / swgt
#valid = n.logical_not(n.isnan(spec))
#spec = spec.compress(valid)

npzfile = '/data3/paper/jaguirre/beamdata/'+srcname+'__'+'srctrack_%.0f-%.0f.npz'%(n.min(times),n.max(times))

afreqs = aa.get_afreqs()#.compress(valid)
#npzfile = src.src_name + '_spec.npz'
print 'Writing spectrum to', npzfile
n.savez(npzfile, spec=spec, freq=afreqs, times=times, x=x, y=y, z=z, mask=mask)
