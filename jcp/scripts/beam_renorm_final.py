#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, optparse, sys
import beamuv

o = optparse.OptionParser()
o.set_usage('srcpredict.py -C [calfile] -s [src] [npzfile]')
o.add_option('-n','--npz',dest='npz',default=None,
    help='The source data npz file to use.')
o.add_option('--nomask', dest='nomask', action='store_true',
    help='Do not use masked npz files.')
o.add_option('--fluxcal',dest='fluxcal',default='cyg',
    help='The source to use as a flux calibrator.')
o.add_option('-b','--beam',dest='beam',default='resbeam.vis.smoothe.npz',
    help='The beam to use.')
a.scripting.add_standard_options(o, cal=True,src=True)
             
opts,args = o.parse_args(sys.argv[1:])

srcs = [f.split('.s__')[0] for f in args]
srcstring = ''
for src in srcs:
    srcstring+=src+','
srcstring=srcstring.rstrip(',')

afreqs = n.load(args[0])['afreqs']
_coeffs = beamuv.coeffs_from_file(opts.beam)
fluxcal = opts.fluxcal

beam = beamuv.BeamUV(_coeffs,.150,size=500,pol='y')
aa = a.cal.get_aa(opts.cal, afreqs)
srclist, cutoff, catalogs = a.scripting.parse_srcs(srcstring, opts.cat)
cat = a.cal.get_catalog(opts.cal,srclist)
cat.compute(aa)

srctimes,srctrack = {},{}
beamtrack,newflux = {},{}
x,y,z, = {},{},{}
scale = 1.

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
        beamtrack[k] = n.array(beamtrack[k])
        beamtrack[k] = n.where(beamtrack[k] < beamtrack[k].max()/10,0,beamtrack[k])
        newflux[k] = n.sum(beamtrack[k]*srctrack[k])/n.sum(beamtrack[k]**2)
        #newflux = n.sum(srctrack[k])/n.sum(beamtrack[k])
        #p.semilogy(srctimes[k],beamtrack[k]*flux,'.')
        #p.semilogy(srctimes[k],srctrack[k],'.')
        p.show()
        if k == fluxcal:
            scale = calflux/newflux[fluxcal]
            print scale

fluxes = []
for k in srcs:
    fluxes.append(newflux[k])
    
n.savez('finalfluxes_sim.npz',srcnames=srcs,srcfluxes=fluxes)

