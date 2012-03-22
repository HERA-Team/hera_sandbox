#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse, random

o = optparse.OptionParser()
o.set_usage('beamcal.py [options]')
a.scripting.add_standard_options(o, cal=True, pol=True, cmap=True, src=True)
o.add_option('-f', '--freq', dest='freq', type='float', default=.150,
    help='Frequency to plot beam for.  Default .150 GHz')
o.add_option('--nside', dest='nside', type='int', default=32,
    help='NSIDE parameter for HEALPix map of beam.')
o.add_option('-o', '--outfile', dest='outfile',
    help='The name of the fits file to create.')
             
opts,args = o.parse_args(sys.argv[1:])

cmap = p.get_cmap(opts.cmap)
print 'Modelling beam at %f GHz' % (opts.freq)

afreqs = n.load(args[0])['afreqs']
srcs = [f.split('__')[0] for f in args]
srcstring = ''
for src in srcs:
    srcstring+=src+','
srcstring=srcstring.rstrip(',')

aa = a.cal.get_aa(opts.cal, afreqs)
srclist, cutoff, catalogs = a.scripting.parse_srcs(srcstring, opts.cat)
cat = a.cal.get_catalog(opts.cal,srclist)
cat.compute(aa)

beam = a.map.Map(opts.nside,interp=True)

fluxcal = 'cyg'
srctimes = {}
srcfluxes = {}
srcgains = {}
for src, npz in zip(srcs,args):
    print 'Reading:', npz
    try: f = n.load(npz)
    except:
        print 'Load file failed.'
        continue
    if not srctimes.has_key(src):
        srctimes[src] = f['times']
    if not srcfluxes.has_key(src):
        srcfluxes[src] = f['spec']
        srcfluxes[src] = n.sum(srcfluxes[src].real,axis=1)
        if src == fluxcal:
            srcgains[src] = f['spec']/n.reshape(cat[src].jys,(1,afreqs.size))
            srcgains[src] = n.sum(srcgains[src].real,axis=1)
            offset = max(srcgains[src])
            #print offset
        else:
            srcgains[src] = f['spec']
            #srcgains[src] = f['spec']/n.reshape(cat[src].jys,(1,afreqs.size))
            srcgains[src] = n.sum(srcgains[src].real,axis=1)

alt = {}
az = {}
ha = {}
dec = {}
x,y,z = {},{},{}
#if fluxcal in srclist:
#    srclist.remove(fluxcal)
#else:
#    print 'You do not have a flux calibrator!'
#    sys.exit()

srclist0 = [s for s in ['cyg','crab','vir','cas'] if s in srcs]
for src in srclist0: srcs.remove(src)
random.shuffle(srcs)
srclist0 += srcs
    
for k in srclist0:
    fluxes = srcfluxes[k]/offset
    gains = srcgains[k]/offset
    if not alt.has_key(k):
        alt[k], az[k] = [],[]
        x[k],y[k],z[k] = [],[],[]
    #if not dec.has_key(k): dec[k] = []
    #if not ha.has_key(k): ha[k] = []
    for i,t in enumerate(srctimes[k]):
        aa.set_jultime(t)
        cat[k].compute(aa)
        #ha[k].append(aa.sidereal_time() - cat[k].ra)
        #dec[k].append(cat[k].dec)
        alt[k].append(cat[k].alt)
        az[k].append(cat[k].az)
        #print fluxes[i],aa.sidereal_time(),cat[k].alt,cat[k].az
    alt[k] = n.array(alt[k])
    az[k] = n.array(az[k])
    #ha[k] = n.array(ha[k])
    #dec[k] = n.array(dec[k])
    #x,y,z = a.coord.radec2eq((ha[k],dec[k]))
    x[k],y[k],z[k] = a.coord.azalt2top((az[k], alt[k]))
    resp = beam.map[x[k],y[k],z[k]]
    weight = beam.wgt[x[k],y[k],z[k]]
    crossing = n.where(weight != 0)[0]
    if n.any(crossing):
        #print crossing
        beamgain = n.sum(beam.map[x[k][crossing],y[k][crossing],z[k][crossing]])/n.sum(beam.wgt[x[k][crossing],y[k][crossing],z[k][crossing]])
        #beamgain = n.average(srcgains['cyg'][crossing])/offset
        print k, beamgain
        f = n.sum((gains[crossing]*beam.wgt[x[k][crossing],y[k][crossing],z[k][crossing]])/((beam.map[x[k][crossing],y[k][crossing],z[k][crossing]])/beam.wgt[x[k][crossing],y[k][crossing],z[k][crossing]]))/n.sum(beam.wgt[x[k][crossing],y[k][crossing],z[k][crossing]])
        print f, n.max(gains)
        gains /= f
    #wgts = n.where(fluxes > 10,n.log10(fluxes),0)
    #wgts = n.where(fluxes > 10, fluxes, 0)
    #wgts = n.ones_like(fluxes)
    wgts = fluxes**2
    #fluxes = n.log10(fluxes)
    beam.add((x[k],y[k],z[k]), wgts, gains)
    if False: #k == 'cyg':
        presp = beam.map[x[k],y[k],z[k]]
        pweight = beam.wgt[x[k],y[k],z[k]]
        p.subplot(211)
        p.plot(presp)
        p.subplot(212)
        p.plot(pweight)
        p.show()
    if True:
        xx,yy,zz = -x[k],-y[k],z[k]
        beam.add((xx,yy,zz),wgts,gains)

beam.to_fits(opts.outfile, clobber=True)
