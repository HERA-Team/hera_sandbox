#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, optparse, sys
#import _beamuv as beamuv
import beamuv

o = optparse.OptionParser()
o.set_usage('srcpredict.py -C [calfile] -s [src] [npzfile]')
o.add_option('-b','--beam',dest='beam',default=None,
    help='The beam npz file to use.')
o.add_option('-n','--npz',dest='npz',default=None,
    help='The source data npz file to use.')
o.add_option('-f','--flux',dest='flux',type='float',default=None,
    help='Manually set the flux of beam track to this value.')
a.scripting.add_standard_options(o, cal=True,src=True)
             
opts,args = o.parse_args(sys.argv[1:])

_coeffs = beamuv.coeffs_from_file(opts.beam)

print _coeffs.shape

k = opts.src
beam = beamuv.BeamUV(_coeffs,.150,size=1000.,pol='y')
aa = a.cal.get_aa(opts.cal, n.array([.150]))
srclist, cutoff, catalogs = a.scripting.parse_srcs(k, opts.cat)
cat = a.cal.get_catalog(opts.cal,srclist)
cat.compute(aa)

for file in args:
    npz = n.load(file)
    alt,az,times = [],[],[]
    for t in npz['times']:
        aa.set_jultime(t)
        times.append(aa.sidereal_time())
        cat[k].compute(aa)
        alt.append(cat[k].alt)
        az.append(cat[k].az)
    alt = n.array(alt)
    az = n.array(az)
    x,y,z = a.coord.azalt2top((az, alt))
    spec = n.mean(npz['spec'].real,axis=1)
    p.semilogy(n.array(times),spec,',')
    peak = n.max(spec)

track = []
for ix,iy,iz in zip(x,y,z):
    track.append(beam.response(n.array([ix]),n.array([iy]),n.array([iz]))**2)

if opts.flux != None:
    track = n.array(track)*opts.flux
    flux = opts.flux
elif opts.npz == None:
    track = n.array(track)*peak/n.max(track)
    flux = peak/n.max(track)
else:
    srcnames = n.load(opts.npz)['srcnames']
    srcfluxes = 10**n.load(opts.npz)['srcfluxes']
    flux = srcfluxes[n.where(srcnames == k)]
    track = n.array(track)*flux
print flux

p.semilogy(times,track,'.',label='y')
p.legend()
p.show()
