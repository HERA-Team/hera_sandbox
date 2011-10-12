#! /usr/bin/env python
import aipy as a, numpy as n, optparse, sys
import capo as C

o = optparse.OptionParser()
o.add_option('--nside', dest='nside', type='int', default=32,
    help='NSIDE parameter for HEALPix map of beam.')
o.add_option('-o', '--outfile', dest='outfile',
    help='The basename of the output files to create.  If not supplied, just plots source tracks.')
o.add_option('-b','--beam',dest='beam',default=None,
    help='The beam file to use (can be npz or fits).')
             
opts,args = o.parse_args(sys.argv[1:])

afreqs = n.load(args[0])['afreqs']
srctimes,srcfluxes,x,y,z = C.jcp.read_srcnpz(args, verbose=True)
for src in srcfluxes: srcfluxes[src] = n.mean(srcfluxes[src], axis=1)

srcs = srctimes.keys()

if opts.beam.endswith('npz'):
    import beamuv
    _coeffs = beamuv.coeffs_from_file(opts.beam)
    beam = beamuv.BeamUV(_coeffs,.150,size=1000)
    response = lambda x,y,z: beam.response(x,y,z)**2
else:
    h = a.map.Map(fromfits=opts.beam)
    h.set_interpol(True)
    #h.set_interpol(False)
    response = lambda x,y,z: h[x,y,z]

newflux = {}

for k in srcs:
    beamtrack = response(x[k], y[k], z[k])
    newflux[k] = n.abs(n.sum(beamtrack*srcfluxes[k])/n.sum(beamtrack**2))
    if not opts.outfile:
        import pylab as p
        p.semilogy(srctimes[k], beamtrack * newflux[k])
        p.semilogy(srctimes[k], srcfluxes[k])

fluxes = n.array([newflux[k] for k in srcs])

for src,flx in zip(srcs,fluxes):
    print 'flux', src, flx

if not opts.outfile:
    p.show()
    sys.exit(0)

print 'Making beams...'
rebeam = a.map.Map(nside=opts.nside,interp=True)

for j,k in enumerate(srcs):
    srcgains = srcfluxes[k]/fluxes[j]
    for i, meas in enumerate(srcgains):
        pcrd,pwgt = rebeam.crd2px(n.array([x[k][i]]),n.array([y[k][i]]),n.array([z[k][i]]),interpolate=True)
        ncrd,nwgt = rebeam.crd2px(n.array([-x[k][i]]),n.array([-y[k][i]]),n.array([z[k][i]]),interpolate=True)
        pcrd.shape,ncrd.shape = (4,),(4,)
        pwgt.shape,nwgt.shape = (4,),(4,)
        rebeam.add(pcrd,pwgt*(fluxes[j])**2,meas)
        rebeam.add(ncrd,nwgt*(fluxes[j])**2,meas)

print 'Saving source-tracks beam to', opts.outfile+'.fits'
rebeam.to_fits(opts.outfile+'.fits',clobber=True)

