#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, optparse, sys
import beamuv

o = optparse.OptionParser()
o.set_usage('srcpredict.py -C [calfile] -s [src] [npzfile]')
o.add_option('-b','--beam',dest='beam',default=None,
    help='The beam npz file to use.')
o.add_option('-n','--npz',dest='npz',default=None,
    help='The source data npz file to use.')
o.add_option('--nomask', dest='nomask', action='store_true',
    help='Do not use masked npz files.')
o.add_option('--fluxcal',dest='fluxcal',default='cyg',
    help='The source to use as a flux calibrator.')
o.add_option('--vis',dest='vis',action='store_true',
    help='Use this flag if you are looking at visibility simulations.')
o.add_option('--ver',dest='ver',default=None,
    help='The version of iterations we are on.')
a.scripting.add_standard_options(o, cal=True,src=True)
             
opts,args = o.parse_args(sys.argv[1:])

if opts.vis:
    srcs = [f.split('.s__')[0] for f in args]
else:
    srcs = [f.split('__')[0] for f in args]
srcstring = ''
for src in srcs:
    srcstring+=src+','
srcstring=srcstring.rstrip(',')

afreqs = n.load(args[0])['afreqs']
_coeffs = beamuv.coeffs_from_file(opts.beam)
fluxcal = opts.fluxcal

beam = beamuv.BeamUV(_coeffs,.150,size=500)
aa = a.cal.get_aa(opts.cal, afreqs)
srclist, cutoff, catalogs = a.scripting.parse_srcs(srcstring, opts.cat)
cat = a.cal.get_catalog(opts.cal,srclist)
cat.compute(aa)

#srcnames = n.load(opts.npz)['srcnames']
#srcfluxes = 10**n.load(opts.npz)['srcfluxes']
srctimes,srctrack = {},{}
beamtrack,newflux = {},{}
x,y,z, = {},{},{}

for k, npz in zip(srcs,args):
    try: f = n.load(npz)
    except:
        print 'Load file failed.'
        continue
    if opts.nomask:
        mask = n.zeros_like(f['times'])
    else: mask = f['mask']
    if not srctimes.has_key(k):
        ftimes = n.ma.array(f['times'],mask=mask)
        srctimes[k] = ftimes.compressed()
    if not srctrack.has_key(k):
        ftrack = f['spec']
        ftrack = n.ma.array(n.mean(ftrack.real,axis=1),mask=mask)
        srctrack[k] = ftrack.compressed()
        if k == fluxcal:
            calflux = n.reshape(cat[fluxcal].jys,(1,afreqs.size))
            calflux = n.mean(calflux.real,axis=1)

            
    alt, az = [],[]
    for t in srctimes[k]:
        aa.set_jultime(t)
        cat[k].compute(aa)
        alt.append(cat[k].alt)
        az.append(cat[k].az)
    alt = n.array(alt)
    az = n.array(az)
    x[k],y[k],z[k] = a.coord.azalt2top((az, alt))
    if not beamtrack.has_key(k):
        beamtrack[k] = []
        for ix,iy,iz in zip(x[k],y[k],z[k]):
            beamtrack[k].append(beam.response(n.array([ix])[0],n.array([iy])[0],n.array([iz])[0])**2)
        flux = srcfluxes[n.where(srcnames == k)]
        beamtrack[k] = n.array(beamtrack[k])
        newflux[k] = n.sum(beamtrack[k]*srctrack[k])/n.sum(beamtrack[k]**2)
        #p.plot(srctrack[k]*beamtrack[k]-(beamtrack[k]**2*newflux[k]))
        #p.plot(srctrack[k])
        #p.plot(beamtrack[k] * newflux[k])
        #p.title(k)
        #p.show()

fluxes = []
for k in srcs:
    fluxes.append(newflux[k])

fluxes = n.array(fluxes)
npzfluxes = n.log10(fluxes)


outnpz = 'sbeam.vis.r3.fits.npz'
print 'saving to ', outnpz
n.savez(outnpz,srcnames=srcs,srcfluxes=npzfluxes)
rebeam = a.map.Map(nside=64,interp=True)

for j,k in enumerate(srcs):
    print k, fluxes[j], fluxes[j]/srcfluxes[j]
    srcgains = srctrack[k]/fluxes[j]
    for i, meas in enumerate(srcgains):
        pcrd,pwgt = rebeam.crd2px(n.array([x[k][i]]),n.array([y[k][i]]),n.array([z[k][i]]),interpolate=True)
        ncrd,nwgt = rebeam.crd2px(n.array([-x[k][i]]),n.array([-y[k][i]]),n.array([z[k][i]]),interpolate=True)
        pcrd.shape,ncrd.shape = (4,),(4,)
        pwgt.shape,nwgt.shape = (4,),(4,)
        rebeam.add(pcrd,pwgt*(fluxes[j])**2,meas)
        rebeam.add(ncrd,nwgt*(fluxes[j])**2,meas)
#outbeam = 're'+opts.beam.replace('.smoothe.npz','.fits')
outbeam = 'sbeam.vis.v3.fits'
print 'saving to ', outbeam
rebeam.to_fits(outbeam,clobber=True)

