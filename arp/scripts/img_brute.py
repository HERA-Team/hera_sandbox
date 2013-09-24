#! /usr/bin/env python

import aipy as a, numpy as n
import sys, optparse
import capo as C

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('--nside', type='int', default=32,
    help='NSIDE of the Healpix Map to use (default 32)')
opts,args = o.parse_args(sys.argv[1:])

CH = 150

uv = a.miriad.UV(args[-1])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
aa.select_chans(n.array([CH]))
aa.set_active_pol(opts.pol)
fq = aa.get_afreqs()[0]
print 'FQ:', fq
del(uv)

h = a.healpix.HealpixMap(nside=opts.nside)
h.set_interpol(False)
px = n.arange(h.npix())
tx,ty,tz = h.px2crd(px)
valid = n.where(tz >= 0, 1, 0)
tx,ty,tz = n.compress(valid,tx), n.compress(valid,ty), n.compress(valid,tz)
bm_a = aa[0].bm_response((tx,ty,tz), pol=opts.pol[0])[0]
bm_b = aa[0].bm_response((tx,ty,tz), pol=opts.pol[-1])[0]
#bm_y = aa[0].bm_response((tx,ty,tz), pol='y')
bm = bm_a * bm_b
print bm.shape

bl_wgts = {}

# M = A * B, B is sky
M = []
A = []

uv = a.miriad.UV(args[-1])
a.scripting.uv_selector(uv, opts.ant, opts.pol)
curtime = None

cat = a.src.get_catalog(srcs=['cen'], catalogs=['misc'])
src = cat['cen']


for (uvw,t,(i,j)),d,f in uv.all(raw=True):
    d,f = d[CH],f[CH]
    if f: continue # bail if no valid data
    if curtime is None:
        curtime = t
        aa.set_jultime(t)
        cat.compute(aa)
        sx,sy,sz = cat.get_crds('top')
        aa.sim_cache(cat.get_crds('eq', ncrd=3), cat.get_jys(), mfreqs=cat.get('mfreq'))
        print sx, sy, sz
    elif t != curtime:
        print 'Inverting'
        print len(A)
        A = n.array(A)
        print A.shape, A.dtype
        B = n.linalg.lstsq(A,n.array(M))[0]
        print 'done'
        h[tx,ty,tz] = B
        h[.8,0,n.sqrt(1-.8**2)] = 100
        import pylab as p
        C.jcp.plot_hpx(h, mode='log', colorbar=True, drng=3)
        p.show()
        sys.exit(0)
        # reset if continuing
    if len(M) > 100: continue
    bl = a.miriad.ij2bl(i,j)
    amp = aa.passband(i,j)
    if not bl_wgts.has_key(bl):
        bx,by,bz = aa.get_baseline(i,j,'z') * fq # in wavelengths
        print 'Caching', (i,j), t, a.miriad.pol2str[uv['pol']], (bx,by,bz)
        # Add the wgts to bl_wgts to cache the answer
        o = aa.get_phs_offset(i,j)
        phs = n.exp(-2j*n.pi * (bx*tx+by*ty+bz*tz + o))
        #phs = n.exp(-2j*n.pi*fq * (bx*tx+by*ty+(bz*tz-1))) # zenith phased
        bl_wgts[bl] = (bm * n.conjugate(phs)).astype(n.complex64)
        #bl_wgts[bl] = (amp * phs).astype(n.float32)
    #img_wgt = n.sqrt(bx**2 + by**2 + bz**2)
    sd = aa.sim(i,j)
    print n.abs(d), n.abs(sd)
    print n.angle(d), n.angle(sd)
    px = n.argmin(n.abs((tx-sx)**2+(ty-sy)**2+(tz-sz)**2))
    print (sx,sy,sz), (tx[px],ty[px],tz[px])
    print 'W', bx*tx[px]+by*ty[px]+bz*tz[px], aa.gen_uvw(i,j, src=src, w_only=True)
    print n.abs(bl_wgts[bl][px])
    print n.angle(bl_wgts[bl][px])
    print '-'*10
    
    M.append(d / amp) # eventually could wgt d by SNR**2 or something
    A.append(bl_wgts[bl])
    
    
    
