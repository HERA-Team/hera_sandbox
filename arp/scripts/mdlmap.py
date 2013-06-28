#! /usr/bin/env python
"""
Models visibilities for various catalog sources and creates a new Miriad UV
file containing either the simulated data, or the residual when the model
is removed from measured data.

Author: Aaron Parsons
"""
import numpy as n, aipy as a, optparse, os, sys, ephem

o = optparse.OptionParser()
o.set_usage('mdlvis.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, cal=True, src=True)
o.add_option('-m','--mode', dest='mode', default='sim',
    help='Operation mode.  Can be "sim" (output simulated data), "sub" (subtract from input data), or "add" (add to input data).  Default is "sim"')
o.add_option('-f', '--flag', dest='flag', action='store_true',
    help='If outputting a simulated data set, mimic the data flagging of the original dataset.')
#o.add_option('-m', '--map', dest='map',
#    help='The Healpix map to use for simulation input.')
#o.add_option('--iepoch', dest='iepoch', default=ephem.J2000, 
#    help='The epoch of coordinates in the map. Default J2000.')
#o.add_option('--freq', dest='freq', default=.150, type='float',
#    help='Frequency of flux data in map.')
o.add_option('-n', '--noiselev', dest='noiselev', default=0., type='float',
    help='RMS amplitude of noise (Jy) added to each UV sample of simulation.')
o.add_option('--nchan', dest='nchan', default=256, type='int',
    help='Number of channels in simulated data if no input data to mimic.  Default is 256')
o.add_option('--sfreq', dest='sfreq', default=.075, type='float',
    help='Start frequency (GHz) in simulated data if no input data to mimic.  Default is 0.075')
o.add_option('--sdf', dest='sdf', default=.150/256, type='float',
    help='Channel spacing (GHz) in simulated data if no input data to mimic.  Default is .150/256')
o.add_option('--inttime', dest='inttime', default=10, type='float',
    help='Integration time (s) in simulated data if no input data to mimic.  Default is 10')
o.add_option('--startjd', dest='startjd', default=2454600., type='float',
    help='Julian Date to start observation if no input data to mimic.  Default is 2454600')
o.add_option('--endjd', dest='endjd', default=2454601., type='float',
    help='Julian Date to end observation if no input data to mimic.  Default is 2454601')
o.add_option('--pol', dest='pol', 
    help='Polarizations to simulate (xx,yy,xy,yx) if starting file from scratch.')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
afreq = aa.get_afreqs()
p,d,f = uv.read(raw=True)
no_flags = n.zeros_like(f)
del(uv)

USE_GAUSS_BEAM = True
PLOT = True

# Generate a model of the sky with point sources and a pixel map

# Initialize pixel map
h = a.healpix.HealpixMap(nside=32)
px = n.arange(h.npix())
eor = {}
for ch in xrange(afreq.size):
    #print ch
    #eor[ch] = n.random.normal(size=h.npix())
    eor[ch] = n.random.randint(2,size=h.npix()).astype(n.float)
#h.put(px, n.ones(h.npix()), n.random.normal(size=h.npix()))
#h.put((0,1,0), 1, 1)
#h.put((1,0,0), 1, 2)
#mflx = h[px]
h.map = eor[100]
x,y,z = h.px2crd(px, ncrd=3)
xyz_eq = n.array((x,y,z))

if PLOT:
    import pylab as P
    P.ion()
    SZ = 250
    im = a.img.Img(size=SZ, res=.5)
    ix,iy,iz = im.get_top((SZ,SZ))
    mask = ix.mask.copy()
    im_xyz_top = n.array([ix,iy,iz])
    im_xyz_top.shape = (im_xyz_top.shape[0], im_xyz_top.shape[1] * im_xyz_top.shape[2])
    plt = None

# A pipe for applying the model
curtime = None
def mdl(uv, p, d, f):
    global plt
    global curtime, eqs
    uvw, t, (i,j) = p
    pol = a.miriad.pol2str[uv['pol']]
    if i == j: return p, d, f
    if curtime != t:
        print t
        curtime = t
        aa.set_jultime(t)
        xyz_top = n.dot(aa.eq2top_m, xyz_eq)
        if USE_GAUSS_BEAM:
            #bm = n.where(xyz_top[2] > 0, n.exp(-(xyz_top[0]**2+xyz_top[1]**2)/(2*.25**2)), 0)
            bm = n.where(xyz_top[2] > 0, aa[0].bm_response(xyz_top, pol='x')[100]**2, 0)
            bm = bm.flatten()
        else:
            bm = aa[0].bm_response(xyz_top, pol='x')**2
    #aa.set_active_pol(pol)
    bl_top = aa.get_baseline(i,j,'z')
    d = n.zeros_like(d)
    for ch in xrange(afreq.size):
        w = n.dot(bl_top, xyz_top) * afreq[ch]
        phs = n.exp(-2j*n.pi*w)
        if USE_GAUSS_BEAM: d[ch] = n.sum(eor[ch] * bm * phs)
        else: d[ch] = n.sum(eor[ch] * bm[ch] * phs)
        if PLOT and ch == 100:
            ix,iy,iz = n.dot(n.linalg.inv(aa.eq2top_m), im_xyz_top)
            #h.set_interpol(True)
            im_eor = h[ix,iy,iz]
            if USE_GAUSS_BEAM: buf,h.map = h.map, bm
            else: buf,h.map = h.map, bm[ch]
            im_eor *= h[ix,iy,iz] # multiply by beam
            h.map = phs
            #im_eor *= h[ix,iy,iz].real # multiply by fringe pattern
            h.map = buf
            im_eor.shape = (2*SZ,2*SZ)
            im_eor = n.where(mask, 0, im_eor)
            #if plt is None: plt = P.imshow(im_eor, vmax=1, vmin=-1, origin='lower')
            if plt is None: plt = P.imshow(im_eor, vmax=im_eor.max(), vmin=0, origin='lower')
            else: plt.set_data(im_eor)
            P.draw()
    if not opts.flag: f = no_flags
    if opts.noiselev != 0:
        # Add on some noise for a more realistic experience
        noise_amp = n.random.random(d.shape) * opts.noiselev
        noise_phs = n.random.random(d.shape) * 2*n.pi * 1j
        noise = noise_amp * n.exp(noise_phs)
        d += noise * aa.passband(i,j)
    return p, n.where(f, 0, d), f

# Run mdl on all files
for filename in args:
    uvofile = filename + 's'
    print filename,'->',uvofile
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    a.scripting.uv_selector(uvi, opts.ant)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mdl, raw=True, append2hist="MDLVIS: " + ' '.join(sys.argv) + '\n')
