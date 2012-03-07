#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, optparse, sys
import beamuv

o = optparse.OptionParser()
o.set_usage('srcpredict.py -C [calfile] -s [src] [npzfile]')
o.add_option('--b1',dest='beam1',default=None,
    help='The first beam npz file to use.')
o.add_option('--b2',dest='beam2',default=None,
    help='The second beam npz file to use.')
o.add_option('--n1',dest='npz1',default=None,
    help='The first source data npz file to use.')
o.add_option('--n2',dest='npz2',default=None,
    help='The second source data npz file to use.')
a.scripting.add_standard_options(o, cal=True,src=True)
             
opts,args = o.parse_args(sys.argv[1:])

_coeffs1 = beamuv.coeffs_from_file(opts.beam1)
_coeffs2 = beamuv.coeffs_from_file(opts.beam2)

k = opts.src
beam1 = beamuv.BeamUV(_coeffs1,.150,size=500,pol='y')
beam2 = beamuv.BeamUV(_coeffs2,.150,size=500,pol='y')

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

track1,track2 = [],[]
for ix,iy,iz in zip(x,y,z):
    track1.append(beam1.response(n.array([ix]),n.array([iy]),n.array([iz]))**2)
    track2.append(beam2.response(n.array([ix]),n.array([iy]),n.array([iz]))**2)

srcnames1 = n.load(opts.npz1)['srcnames']
srcfluxes1 = 10**n.load(opts.npz1)['srcfluxes']
flux1 = srcfluxes1[n.where(srcnames1 == k)][0]
track1 = n.array(track1)*flux1

srcnames2 = n.load(opts.npz2)['srcnames']
srcfluxes2 = 10**n.load(opts.npz2)['srcfluxes']
flux2 = srcfluxes2[n.where(srcnames2 == k)][0]
track2 = n.array(track2)*flux2

p.semilogy(times,track1,'.',label='1')
p.semilogy(times,track2,'.',label='2')
p.legend()
p.show()
