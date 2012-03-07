#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse

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
        srcfluxes[src] = n.mean(srcfluxes[src].real,axis=1)
        if src == fluxcal:
            calflux = n.reshape(cat[fluxcal].jys,(1,afreqs.size))
            calflux = n.mean(calflux.real,axis=1)

print calflux

alt, az = {},{}
x,y,z = {},{},{}

print 'Calculating source tracks...'
track = {}
tracks = n.array([]) 
for k in srcs:
    if not alt.has_key(k):
        alt[k], az[k] = [],[]
        x[k],y[k],z[k] = [],[],[]
    for i,t in enumerate(srctimes[k]):
        aa.set_jultime(t)
        cat[k].compute(aa)
        alt[k].append(cat[k].alt)
        az[k].append(cat[k].az)
    alt[k] = n.array(alt[k])
    az[k] = n.array(az[k])
    x[k],y[k],z[k] = a.coord.azalt2top((az[k], alt[k]))
    track[k] = n.append(n.unique(beam.crd2px(x[k],y[k],z[k])), n.unique(beam.crd2px(-x[k],-y[k],z[k])))
    #track[k] = n.unique(n.append(beam.crd2px(x[k],y[k],z[k]), beam.crd2px(-x[k],-y[k],z[k])))
    tracks = n.append(tracks,track[k])

tracks = tracks.astype(n.long)

print 'Determining crossing points...'
cnt = {}
for i in xrange(max(tracks)+1):
    cnt[i] = n.where(tracks == i)[0].shape[0]
    crossing_pixels = n.where(n.array(cnt.values()) > 1)[0]

C,S,M,wM = [],[],[],[]
for k in srcs:
    for i,meas in enumerate(srcfluxes[k]):
        pcrd = beam.crd2px(n.array([x[k][i]]),n.array([y[k][i]]),n.array([z[k][i]]))
        ncrd = beam.crd2px(n.array([-x[k][i]]),n.array([-y[k][i]]),n.array([z[k][i]]))
        if pcrd in crossing_pixels:
            S.append(k)
            wM.append(n.log10(meas)*(meas**2))
            M.append(meas)
            C.append(pcrd[0])
        if ncrd in crossing_pixels:
            S.append(k)
            wM.append(n.log10(meas)*(meas**2))
            M.append(meas)
            C.append(ncrd[0])

print 'Constructing matrices...'
dC,dS = {},{}
for i,c in enumerate(crossing_pixels):
    dC[c] = i
for i,k in enumerate(srcs):
    if not dS.has_key(k): dS[k] = i

nmeas = len(M)
npix = len(crossing_pixels)
nsrcs = len(srcs)
A = n.zeros((len(M)+1,npix+nsrcs))


for k in srcs:
    for ind,meas in enumerate(M):
        A[ind,dC[C[ind]]] = 1.*meas**2
        A[ind,npix+dS[S[ind]]] = 1.*meas**2
    if k == fluxcal:
        A[-1,npix+dS[fluxcal]] = 1e16

wM.append(1e16*n.log10(calflux))
wM = n.array(wM)

#p.imshow(A,aspect='auto',interpolation='nearest',vmax=1)
#p.show()

print 'Solving equation...'
#B = n.linalg.inv(A) * M
B = n.linalg.lstsq(A,wM)
n.savez(opts.outfile+'.npz',B[0])
#print B
print 10**B[0]

print 'Making beam...'
bm = 10**(B[0][:-nsrcs])
#bm = 10**(B[0][:-nsrcs]-n.max(B[0][:-nsrcs]))
beam.add(crossing_pixels,n.ones_like(bm),bm)
beam.to_fits(opts.outfile, clobber=True)

