#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse

o = optparse.OptionParser()
o.set_usage('beamcal.py [options]')
a.scripting.add_standard_options(o, cal=True, pol=True, cmap=True, src=True)
o.add_option('-f', '--freq', dest='freq', type='float', default=.150,
    help='Frequency to plot beam for.  Default .150 GHz')
o.add_option('--nside', dest='nside', type='int', default=32,
    help='NSIDE parameter for HEALPix map of beam.')
o.add_option('--fluxes',dest='fluxes',type='string',
    help='The .npz file with source fluxes.')
o.add_option('--nomask', dest='nomask', action='store_true',
    help='Do not use masked npz files.')
o.add_option('--src1',dest='src1',type='string',
    help='First of two sources.')
o.add_option('--src2',dest='src2',type='string',
    help='Second of two sources.')
             
opts,args = o.parse_args(sys.argv[1:])

cmap = p.get_cmap(opts.cmap)

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
    #track[k] = n.append(n.unique(beama.crd2px(x[k],y[k],z[k],interpolate=True)[0]), n.unique(beama.crd2px(-x[k],-y[k],z[k],interpolate=True)[0]))
    track[k] = n.unique(n.append(beama.crd2px(x[k],y[k],z[k],interpolate=True)[0], beama.crd2px(-x[k],-y[k],z[k],interpolate=True)[0]))
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
        w.append((fluxtrack[k][crd]**2)*(wgttrack[k][crd]**2)/(sqwgttrack[k][crd]))
        wM.append(n.log10(fluxtrack[k][crd])*(fluxtrack[k][crd]**2)*(wgttrack[k][crd]**2)/(sqwgttrack[k][crd]))

npz = n.load(opts.fluxes)
names = npz['srcnames']
fluxes = npz['srcfluxes']
s1pos = n.where(names == opts.src1)
s2pos = n.where(names == opts.src2)

alt,az = a.coord.xyz2thphi(beama.px2crd(crossing_pixels))/a.const.deg
for i,pix in enumerate(crossing_pixels):
    #print opts.src1,pix,fluxtrack[opts.src1][pix]/(10**fluxes[s1pos])
    #print opts.src2,pix,fluxtrack[opts.src2][pix]/(10**fluxes[s2pos])
    #beama.add(pix,1,(fluxtrack[opts.src1][pix]/(10**fluxes[s1pos]))/(fluxtrack[opts.src2][pix]/(10**fluxes[s2pos])))
    print opts.src1, pix, alt[i], az[i]
    #print opts.src2, pix, alt[i], az[i]

    
#beama.to_fits('beamtest.fits',clobber=True)
