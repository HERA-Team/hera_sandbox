#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse, pickle
import capo as C

o = optparse.OptionParser()
o.set_usage('beamcal.py [options]')
a.scripting.add_standard_options(o, cal=True)
o.add_option('--cat', dest='cat', default='helm,misc',
    help='A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aip y._src.  Some catalogs may require a separate data file to be downloaded and installed.')
o.add_option('--nside', dest='nside', type='int', default=32,
    help='NSIDE parameter for HEALPix map of beam.')
o.add_option('-o', '--outfile', dest='outfile',
    help='The basename of the output files to create.')
o.add_option('--nomask', dest='nomask', action='store_true',
    help='Do not use masked npz files.')
o.add_option('--pickle',dest='pickle',action='store_true',
    help='Store measurement dictionaries in a pickle file.')
o.add_option('--fluxcal',dest='fluxcal',default='cyg',
    help='The source to use as a flux calibrator.')
             
opts,args = o.parse_args(sys.argv[1:])

afreqs = n.load(args[0])['afreqs']
aa = a.cal.get_aa(opts.cal, afreqs)
srctimes,srcfluxes = C.jcp.read_srcnpz(args, verbose=True)
for src in srcfluxes: srcfluxes[src] = n.mean(srcfluxes[src], axis=1)

srcs = srctimes.keys()
srclist, cutoff, catalogs = a.scripting.parse_srcs(','.join(srcs), opts.cat)
cat = a.cal.get_catalog(opts.cal,srclist)
cat.compute(aa)

beama = a.map.Map(opts.nside,interp=True)
beamb = a.map.Map(opts.nside,interp=True)

calflux = n.reshape(cat[opts.fluxcal].jys,(1,afreqs.size))
calflux = n.mean(calflux.real,axis=1)

srcgains = {}

x,y,z = {},{},{}

print 'Calculating source tracks...'
track = {}
tracks = n.array([]) 
for k in srcs:
    for t in srctimes[k]:
        aa.set_jultime(t)
        cat[k].compute(aa)
        xi,yi,zi = cat[k].get_crds('top')
        x[k] = x.get(k,[]) + [xi]
        y[k] = y.get(k,[]) + [yi]
        z[k] = z.get(k,[]) + [zi]
    x[k],y[k],z[k] = n.array(x[k]), n.array(y[k]), n.array(z[k])
    #valid = n.where(z[k] > 0, 1, 0)
    #x[k], y[k], z[k] = n.compress(valid, x[k]), n.compress(valid,y[k]), n.compress(valid, z[k])
    track[k] = n.append(n.unique(beama.crd2px(x[k],y[k],z[k],interpolate=True)[0]), n.unique(beama.crd2px(-x[k],-y[k],z[k],interpolate=True)[0]))
    #track[k] = n.unique(n.append(beama.crd2px(x[k],y[k],z[k]), beam.crd2px(-x[k],-y[k],z[k])))
    tracks = n.append(tracks,track[k])
tracks = tracks.astype(n.long)

print 'Determining crossing points...'
cnt = {}
for i in xrange(max(tracks)+1):
    cnt[i] = n.where(tracks == i)[0].shape[0]
    crossing_pixels = n.where(n.array(cnt.values()) > 1)[0]

print 'Averaging measurements within a pixel...'
fluxtrack,wgttrack,sqwgttrack = {},{},{}
S,C,w,wM = [],[],[],[]
for k in srcs:
    fluxtrack[k],wgttrack[k],sqwgttrack[k] = {},{},{}
    for i, meas in enumerate(srcfluxes[k]):
        pcrd,pwgt = beama.crd2px(n.array([x[k][i]]),n.array([y[k][i]]),n.array([z[k][i]]),interpolate=True)
        ncrd,nwgt = beama.crd2px(n.array([-x[k][i]]),n.array([-y[k][i]]),n.array([z[k][i]]),interpolate=True)
        for ind,crd in enumerate(pcrd[0]):
            if crd in crossing_pixels:
                if not fluxtrack[k].has_key(crd):
                    fluxtrack[k][crd] = pwgt[0][ind]*meas
                    wgttrack[k][crd] = pwgt[0][ind]
                    sqwgttrack[k][crd] = (pwgt[0][ind])**2
                else:
                    fluxtrack[k][crd] += pwgt[0][ind]*meas
                    wgttrack[k][crd] += pwgt[0][ind]
                    sqwgttrack[k][crd] += (pwgt[0][ind])**2
        for ind,crd in enumerate(ncrd[0]):
            if crd in crossing_pixels:
                if not fluxtrack[k].has_key(crd):
                    fluxtrack[k][crd] = nwgt[0][ind]*meas
                    wgttrack[k][crd] = nwgt[0][ind]
                    sqwgttrack[k][crd] = (nwgt[0][ind])**2
                else:
                    fluxtrack[k][crd] += nwgt[0][ind]*meas
                    wgttrack[k][crd] += nwgt[0][ind]
                    sqwgttrack[k][crd] += (nwgt[0][ind])**2
    for crd in fluxtrack[k].keys():
        S.append(k)
        C.append(crd)
        fluxtrack[k][crd] /= wgttrack[k][crd]
        weight = (fluxtrack[k][crd]**2)*(wgttrack[k][crd])
        w.append(weight)
        wM.append(n.log10(fluxtrack[k][crd])*(fluxtrack[k][crd]**2)*(wgttrack[k][crd]))
        wgttrack[k][crd] = weight

#print n.where(n.array(w)==0)
print 'Constructing matrices...'
dC,dS = {},{}
for i,c in enumerate(crossing_pixels):
    dC[c] = i
for i,k in enumerate(srcs):
    if not dS.has_key(k): dS[k] = i
neq = len(wM)
npix = len(crossing_pixels)
nsrcs = len(srcs)
A = n.zeros((neq+1,npix+nsrcs),dtype=n.float32)

for k in srcs:
    for ind,wgt in enumerate(w):
        A[ind,dC[C[ind]]] = wgt
        A[ind,npix+dS[S[ind]]] = wgt
    if k == opts.fluxcal:
        A[-1,npix+dS[opts.fluxcal]] = 1e16

wM.append(1e16*n.log10(calflux))
wM = n.array(wM,dtype=n.float32)

measurements = {'wgt':wgttrack,'meas':fluxtrack}

print 'Solving equation...'
#B = n.dot(n.linalg.inv(n.dot(n.transpose(A),A)),n.dot(n.transpose(A),wM))
B = n.linalg.lstsq(A,wM)
print 'Saving source info to', opts.outfile+'.npz'
n.savez(opts.outfile+'.npz',solution=B[0],srcnames=srcs,srcfluxes=B[0][-nsrcs:])
for src,flx in zip(srcs,B[0][-nsrcs:]):
    print src, 10**flx

print 'Making beams...'
#bm = 10**(B[:-nsrcs])
bm = 10**(B[0][:-nsrcs])
beama.add(crossing_pixels,n.ones_like(bm),bm)
outnamea = opts.outfile + 'a.fits'
print 'Saving crossing-points beam to', outnamea
beama.to_fits(outnamea, clobber=True)

pname = opts.outfile + '.pkl'
output = open(pname,'wb')
pickle.dump(measurements,output)
output.close

outnameb = opts.outfile + 'b.fits'
fluxes = 10**(B[0][npix:])
for j,k in enumerate(srcs):
    srcgains = srcfluxes[k]/fluxes[j]
    for i, meas in enumerate(srcgains):
        pcrd,pwgt = beama.crd2px(n.array([x[k][i]]),n.array([y[k][i]]),n.array([z[k][i]]),interpolate=True)
        ncrd,nwgt = beama.crd2px(n.array([-x[k][i]]),n.array([-y[k][i]]),n.array([z[k][i]]),interpolate=True)
        pcrd.shape,ncrd.shape = (4,),(4,)
        pwgt.shape,nwgt.shape = (4,),(4,)
        beamb.add(pcrd,pwgt*(fluxes[j])**2,meas)
        beamb.add(ncrd,nwgt*(fluxes[j])**2,meas)
print 'Saving source-tracks beam to', outnameb
beamb.to_fits(outnameb,clobber=True)
