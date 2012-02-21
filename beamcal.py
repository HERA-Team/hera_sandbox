#! /usr/bin/env python
import aipy as a, numpy as n, sys, optparse#, pickle
import capo as C

o = optparse.OptionParser()
o.add_option('--nside', dest='nside', type='int', default=64,
    help='NSIDE parameter for HEALPix map of beam.')
o.add_option('--no_interp', dest='no_interp', action='store_true',
    help='Do not use sub-pixel interpolation when gridding sources.')
o.add_option('-o', '--outfile', dest='outfile',
    help='The basename of the output files to create.')
o.add_option('--fluxcal',dest='fluxcal',default='cyg,10622.92',
    help='The source,flux to use for calibration.  Default is "cyg,10622.92"')
opts,args = o.parse_args(sys.argv[1:])

afreqs = n.load(args[0])['freq']
srctimes,srcfluxes,srcwgts,x,y,z = C.jcp.read_srcnpz(args, verbose=True)
srcs = srctimes.keys()
if True:
    bigsrcs = []
    for s in srcs:
        peak = n.abs(n.max(srcfluxes[s]) / n.max(srcwgts[s]))
        if peak > 25: bigsrcs.append(s)
    srcs = bigsrcs
    print srcs

bm = a.map.Map(opts.nside,interp=True)

calsrc,calflux = opts.fluxcal.split(',')
calflux = float(calflux)

# Gridding source tracks
src_bmtrk = {}
sum_src = 0
for k in srcs:
    xk,yk,zk = n.concatenate([x[k],-x[k]]), n.concatenate([y[k],-y[k]]), n.concatenate([z[k],z[k]])
    # XXX Should throw away a hemisphere here for efficiency (must get the same answer top/bottom)
    flx = n.concatenate([srcfluxes[k], srcfluxes[k]])
    wgt = n.concatenate([srcwgts[k], srcwgts[k]])
    if opts.no_interp:
        px = bm.crd2px(xk,yk,zk, interpolate=False)
    else:
        px,_wgt = bm.crd2px(xk,yk,zk, interpolate=True)
        px,_wgt = px.flatten(), _wgt.flatten()
        flx = n.array([flx,flx,flx,flx]).transpose().flatten()
        wgt = n.array([wgt,wgt,wgt,wgt]).transpose().flatten() * _wgt
    bm.add(px, wgt, flx)
    src_bmtrk[k] = (bm.map.map.copy(), bm.wgt.map.copy())
    bm.map.map *= 0; bm.wgt.map *= 0
    # XXX currently, this doesn't find self-crossing points
    sum_src += n.where(src_bmtrk[k][1] > 0, 1., 0)

crossing_pixels = n.where(sum_src > 1)[0]
print 'Solving for %d crossing pixels' % crossing_pixels.size

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
        if m <= 0: continue # skip if source measurement is negative (logrithm explodes).  This will bias answer
        logm, logw = n.log10(m/w), (m/w)**2 * w
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
    
A = n.array(A)
print 'Solving equation (M = A*B) for B: inverting %s matrix' % str(A.shape)
B = n.linalg.lstsq(A,n.array(M))[0]
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
    outmap += mp * (wt * flx)
    outwgt += (wt * flx)**2

bm.map.map, bm.wgt.map = outmap, outwgt
print 'Saving source-tracks beam to', outnameb
bm.to_fits(outnameb,clobber=True)
