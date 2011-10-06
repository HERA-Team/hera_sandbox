#! /usr/bin/env python
import aipy as a, numpy as n, sys, optparse#, pickle
import capo as C

o = optparse.OptionParser()
o.add_option('--nside', dest='nside', type='int', default=32,
    help='NSIDE parameter for HEALPix map of beam.')
o.add_option('--no_interp', dest='no_interp', action='store_true',
    help='Do not use sub-pixel interpolation when gridding sources.')
o.add_option('-o', '--outfile', dest='outfile',
    help='The basename of the output files to create.')
o.add_option('--fluxcal',dest='fluxcal',default='cyg,10622.92',
    help='The source,flux to use for calibration.  Default is "cyg,10622.92"')
opts,args = o.parse_args(sys.argv[1:])

afreqs = n.load(args[0])['afreqs']
srctimes,srctracks,x,y,z = C.jcp.read_srcnpz(args, verbose=True)
for src in srctracks: srctracks[src] = n.mean(srctracks[src], axis=1)
srcs = srctimes.keys()

bm = a.map.Map(opts.nside,interp=True)

calsrc,calflux = opts.fluxcal.split(',')
calflux = float(calflux)

# Gridding source tracks
src_bmtrk = {}
sum_src = 0
for k in srcs:
    xk,yk,zk = n.concatenate([x[k],-x[k]]), n.concatenate([y[k],-y[k]]), n.concatenate([z[k],z[k]])
    flx = n.concatenate([srctracks[k], srctracks[k]])
    if opts.no_interp:
        px,wgt = bm.crd2px(xk,yk,zk, interpolate=False), n.ones_like(xk)
    else:
        px,wgt = bm.crd2px(xk,yk,zk, interpolate=True)
        px,wgt = px.flatten(), wgt.flatten()
        flx = n.array([flx,flx,flx,flx]).transpose().flatten()
        #valid = n.where(wgt < .1, 0, 1)
        #px,wgt,flx = px.compress(valid), wgt.compress(valid), flx.compress(valid)
    bm.add(px, wgt, flx)
    src_bmtrk[k] = (bm.map.map.copy(), bm.wgt.map.copy())
    bm.map.map *= 0; bm.wgt.map *= 0
    sum_src += n.where(src_bmtrk[k][1] > 0, 1., 0)

crossing_pixels = n.where(sum_src > 1)[0]

# Creating matrices for least-squares inversion
# B matrix is parameter values (beam pixels, source fluxes) to solve for
A = []  # matrix implementing the equations
M = []  # weighted log measured answers
for j, src in enumerate(srcs):
    for i, px in enumerate(crossing_pixels):
        i += len(srcs)
        mp,wt = src_bmtrk[src]
        if wt[px] == 0: continue  # skip if crossing point isn't in this source track
        m,w = mp[px], wt[px]
        logm, logw = n.log10(m/w), m**2 * w
        Aline = n.zeros(len(srcs) + crossing_pixels.size, dtype=n.float)
        Aline[j],Aline[i] = logw, logw
        A.append(Aline)
        M.append(logm * logw)
# Add equation fixing flux of calibrator source
j = srcs.index(calsrc)
Aline = n.zeros(len(srcs) + crossing_pixels.size, dtype=n.float)
logm,logw = n.log10(calflux), 1e16
Aline[j] = logw
A.append(Aline)
M.append(logm * logw)
    
print 'Solving equation (M = A*B) for B...'
B = n.linalg.lstsq(n.array(A),n.array(M))[0]
src_fluxes = 10**B[:len(srcs)]
bm_resps = 10**B[len(srcs):]

assert(n.all(bm.wgt.map == 0))  # Make sure beam is fresh and clean
bm.add(crossing_pixels,n.ones_like(bm_resps),bm_resps)
outnamea = opts.outfile + 'a.fits'
print 'Saving crossing-points beam to', outnamea
bm.to_fits(outnamea, clobber=True)

outnameb = opts.outfile + 'b.fits'
outmap,outwgt = 0,0
for src,flx in zip(srcs, src_fluxes):
    print 'flux', src, flx
    mp,wt = src_bmtrk[src]
    outmap += mp * flx
    outwgt += wt * flx**2

bm.map.map, bm.wgt.map = outmap, outwgt
print 'Saving source-tracks beam to', outnameb
bm.to_fits(outnameb,clobber=True)
