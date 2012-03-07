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

beama = a.map.Map(opts.nside,interp=True)
beamb = a.map.Map(opts.nside,interp=True)
beamc = a.map.Map(opts.nside,interp=True)

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
track,wgt = {},{}
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
    track[k] = n.append(n.unique(beama.crd2px(x[k],y[k],z[k],interpolate=True)[0]), n.unique(beama.crd2px(-x[k],-y[k],z[k],interpolate=True)[0]))
    #track[k] = n.unique(n.append(beama.crd2px(x[k],y[k],z[k]), beam.crd2px(-x[k],-y[k],z[k])))
    tracks = n.append(tracks,track[k])

tracks = tracks.astype(n.long)

print 'Determining crossing points...'
cnt = {}
for i in xrange(max(tracks)+1):
    cnt[i] = n.where(tracks == i)[0].shape[0]
    crossing_pixels = n.where(n.array(cnt.values()) > 1)[0]

C,S,M,w,wM = [],[],[],[],[]
for k in srcs:
    for i,meas in enumerate(srcfluxes[k]):
        pcrd,pwgt = beama.crd2px(n.array([x[k][i]]),n.array([y[k][i]]),n.array([z[k][i]]),interpolate=True)
        ncrd,nwgt = beama.crd2px(n.array([-x[k][i]]),n.array([-y[k][i]]),n.array([z[k][i]]),interpolate=True)
        for ind,crd in enumerate(pcrd[0]):
            if crd in crossing_pixels:
                S.append(k)
                wM.append(n.log10(meas)*(meas**2)*(pwgt[0][ind]**2))
                M.append(meas)
                w.append(pwgt[0][ind])
                C.append(crd)
        for ind,crd in enumerate(ncrd[0]):
            if crd in crossing_pixels:
                S.append(k)
                wM.append(n.log10(meas)*(meas**2)*(nwgt[0][ind]**2))
                w.append(nwgt[0][ind])
                M.append(meas)
                C.append(crd)

print 'Constructing matrices...'
dC,dS = {},{}
for i,c in enumerate(crossing_pixels):
    dC[c] = i
for i,k in enumerate(srcs):
    if not dS.has_key(k): dS[k] = i

nmeas = len(M)
npix = len(crossing_pixels)
nsrcs = len(srcs)


print nmeas,npix,nsrcs

A = n.zeros((len(M)+1,npix+nsrcs),dtype=n.float32)


for k in srcs:
    for ind,meas in enumerate(M):
        A[ind,dC[C[ind]]] = (1.*(meas**2)*(w[ind]**2))
        A[ind,npix+dS[S[ind]]] = 1.*(meas**2)*(w[ind]**2)
        #print type(A[ind,npix+dS[S[ind]]])
    if k == fluxcal:
        A[-1,npix+dS[fluxcal]] = 1e16

wM.append(1e16*n.log10(calflux))
wM = n.array(wM,dtype=n.float32)

#p.imshow(A,aspect='auto',interpolation='nearest')
#p.colorbar()
#p.show()

print 'Solving equation...'
#B = n.linalg.inv(A) * M
B = n.linalg.lstsq(A,wM)
n.savez(opts.outfile+'.npz',B[0])
#print B
print 10**B[0]

print 'Making beams...'
bm = 10**(B[0][:-nsrcs])
#bm = 10**(B[0][:-nsrcs]-n.max(B[0][:-nsrcs]))
beama.add(crossing_pixels,n.ones_like(bm),bm)
if beama[0,0,1] != 0:
    print 'Normalizing...'
    for pix,gain in enumerate(beama):
        beama.put(n.array([pix]),beama.wgt[pix],gain/beama[0,0,1])
basename = opts.outfile.split('.')
outnamea = basename[0]+'a.'+basename[1]
beama.to_fits(outnamea, clobber=True)

if True:
    outnameb = basename[0]+'b.'+basename[1]
    outnamec = basename[0]+'c.'+basename[1]
    fluxes = 10**(B[0][npix:])
    for j,k in enumerate(srcs):
        srcgains = srcfluxes[k]/fluxes[j]
        beamb.add((x[k],y[k],z[k]),n.ones_like(srcgains),srcgains)
        beamb.add((-x[k],-y[k],z[k]),n.ones_like(srcgains),srcgains)
        beamc.add((x[k],y[k],z[k]),srcfluxes[k]**2,srcgains)
        beamc.add((-x[k],-y[k],z[k]),srcfluxes[k]**2,srcgains)
    beamb.to_fits(outnameb,clobber=True)
    beamc.to_fits(outnamec,clobber=True)
