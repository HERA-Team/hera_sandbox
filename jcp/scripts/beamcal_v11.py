#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse, pickle

o = optparse.OptionParser()
o.set_usage('beamcal.py [options]')
a.scripting.add_standard_options(o, cal=True, pol=True, cmap=True, src=True)
o.add_option('-f', '--freq', dest='freq', type='float', default=.150,
    help='Frequency to plot beam for.  Default .150 GHz')
o.add_option('--nside', dest='nside', type='int', default=32,
    help='NSIDE parameter for HEALPix map of beam.')
o.add_option('-o', '--outfile', dest='outfile',
    help='The name of the fits file to create.')
o.add_option('--nomask', dest='nomask', action='store_true',
    help='Do not use masked npz files.')
o.add_option('--pickle',dest='pickle',action='store_true',
    help='Store measurement dictionaries in a pickle file.')
o.add_option('--fluxcal',dest='fluxcal',default='cyg',
    help='The source to use as a flux calibrator.')
             
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


fluxcal = opts.fluxcal
srctimes = {}
srcfluxes = {}
srcgains = {}
for src, npz in zip(srcs,args):
    print 'Reading:', npz
    try: f = n.load(npz)
    except:
        print 'Load file failed.'
        continue
    if opts.nomask: 
        mask = n.zeros_like(f['times'])
    else: mask = f['mask']
    if not srctimes.has_key(src):
        ftimes = n.ma.array(f['times'],mask=mask)
        srctimes[src] = ftimes.compressed()
    if not srcfluxes.has_key(src):
        ffluxes = f['spec']
        ffluxes = n.ma.array(n.mean(ffluxes.real,axis=1),mask=mask)
        srcfluxes[src] = ffluxes.compressed()
        if src == fluxcal:
            calflux = n.reshape(cat[fluxcal].jys,(1,afreqs.size))
            calflux = n.mean(calflux.real,axis=1)
            #print calflux

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
        weight = (fluxtrack[k][crd]**2)*(wgttrack[k][crd]**2)/(sqwgttrack[k][crd])
        w.append(weight)
        wM.append(n.log10(fluxtrack[k][crd])*(fluxtrack[k][crd]**2)*(wgttrack[k][crd]**2)/(sqwgttrack[k][crd]))
        wgttrack[k][crd] = weight
        #w.append((fluxtrack[k][crd]**2)*(wgttrack[k][crd]))
        #wM.append(n.log10(fluxtrack[k][crd])*(fluxtrack[k][crd]**2)*(wgttrack[k][crd]))
        #w.append((fluxtrack[k][crd]**2)*(wgttrack[k][crd])/100**2)
        #wM.append(n.log10(fluxtrack[k][crd])*(fluxtrack[k][crd]**2)*(wgttrack[k][crd])/100**2)

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
    if k == fluxcal:
        A[-1,npix+dS[fluxcal]] = 1e16

wM.append(1e16*n.log10(calflux))
wM = n.array(wM,dtype=n.float32)

measurements = {'wgt':wgttrack,'meas':fluxtrack}

print 'Solving equation...'
#B = n.dot(n.linalg.inv(n.dot(n.transpose(A),A)),n.dot(n.transpose(A),wM))
B = n.linalg.lstsq(A,wM)
n.savez(opts.outfile+'.npz',solution=B[0],srcnames=srcs,srcfluxes=B[0][-nsrcs:])
print 10**B[0]
#print A, wM, 10**B[0]

#p.imshow(A,aspect='auto',interpolation='nearest',vmax=1.5)
#p.colorbar(shrink=0.5)
#p.show()

print 'Making beams...'
#bm = 10**(B[:-nsrcs])
bm = 10**(B[0][:-nsrcs])
beama.add(crossing_pixels,n.ones_like(bm),bm)
if beama[0,0,1] != 0:
    norm = (beama[0]+beama[1]+beama[2]+beama[3])/4
    print 'Normalizing BeamA...'
    beama.map.map/=norm
basename = opts.outfile.split('.')
outnamea = basename[0]+'a.'+basename[1]
beama.to_fits(outnamea, clobber=True)

pname = basename[0]+'.'+basename[1]+'.pkl'
output = open(pname,'wb')
pickle.dump(measurements,output)
output.close

outnameb = basename[0]+'b.'+basename[1]
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
if beamb[0,0,1] != 0:
    norm = (beamb[0]+beamb[1]+beamb[2]+beamb[3])/4
    print 'Normalizing BeamB...'
    beamb.map.map/=norm
beamb.to_fits(outnameb,clobber=True)
