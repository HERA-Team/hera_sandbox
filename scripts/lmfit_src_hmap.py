#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse, re
import glob, os#, pymc
from mpl_toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_bm_hmap.py [options] *.hmap')
a.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('--res', dest='res', type='float', default=.01,
    help='Resolution of plot (in radians).  Default .01')
o.add_option('--lmax', dest='lmax', type='int', default=20,
    help='Smooth to the specified maximum spherical harmonic.  Default 20')
o.add_option('--ns', dest='ns', action='store_true',
    help='Assume North/South beam symmetry.')
o.add_option('--ew', dest='ew', action='store_true',
    help='Assume East/West beam symmetry.')
o.add_option('--rot', dest='rot', action='store_true',
    help='Assume 180-degree rotational beam symmetry.')
o.add_option('--hmap', dest='hmap',
    help='Initial beam model.')
o.add_option('--altmin', dest='altmin', type='float', default=25.,
    help='Minimum altitude (in degrees) to use.  Default 25')
opts,args = o.parse_args(sys.argv[1:])

assert(opts.pol in ['xx','yy'])

re_freq = re.compile(r'_(\d+)\.hmap$')
fq = float(re_freq.search(opts.hmap).groups()[0]) / 1e3
hmap = a.map.Map(fromfits=opts.hmap)
zen_resp = hmap[0,0,1][0]
hmap.map.map /= zen_resp

srclist = [f.split('.')[-1][3:-1] for f in args]
cat = a.src.get_catalog(srclist)
uv = a.miriad.UV(args[0])
ch = int(n.round((fq - uv['sfreq']) / uv['sdf']))
aa = a.cal.get_aa(opts.cal, .01, fq, 1)
cat.compute(aa)
del(uv)

print 'Gathering source data'
src_dat = {}
prms = {}
for src_name, filename in zip(srclist, args):
    print '    ', filename
    src = cat[src_name]
    if src.get_jys()[0] < .01: continue
    uv = a.miriad.UV(filename)
    dat = []
    x,y,z = [],[],[]
    for (crd,t,bl),d,f in uv.all(raw=True):
        if f[ch]: continue
        aa.set_jultime(t)
        src.compute(aa)
        if src.alt < opts.altmin * a.img.deg2rad: continue
        xi,yi,zi = src.get_crds('top', ncrd=3)
        x.append(xi); y.append(yi); z.append(zi)
        dat.append(d[ch])
    sdat = n.abs(n.array(dat))
    x,y,z = n.array(x), n.array(y), n.array(z)
    poly = n.polyfit(x, sdat, deg=6)
    sig = n.sqrt(n.average((sdat - n.polyval(poly, x))**2))
    jys = src.get_jys()[0]
    src_dat[src_name] = (x,y,z,sdat,jys,sig)
    #prms[src_name] = pymc.Normal(src_name, jys, 1./sig**2)
    prms[src_name] = jys
    #print src_name, src_dat[src_name][-2:]

#@pymc.deterministic
def infer_hmap(arg=prms):
    alms = a.healpix.Alm(opts.lmax,opts.lmax)
    _hmapr = a.healpix.HealpixMap(hmap.nside(), dtype=n.float)
    _hmapi = a.healpix.HealpixMap(hmap.nside(), dtype=n.float)
    _alm = a.healpix.Alm(opts.lmax,opts.lmax)
    for L,M in zip(*alms.lm_indices()):
        _alm.set_to_zero()
        _alm[L,M] = 1.0
        _hmapr.from_alm(_alm)
        _alm[L,M] = 1.0j
        _hmapi.from_alm(_alm)
        csumr, csumi, cwgtr, cwgti = 0, 0, 0, 0
        for src_name in arg.keys():
            x,y,z,sdat,jys,sig = src_dat[src_name]
            if opts.pol.startswith('x'): _x,_y = x,y
            else: _x,_y = -y,x
            #jys = arg[src_name]
            _hdatr = _hmapr[_x,_y,z]
            _hdati = _hmapi[_x,_y,z]
            _ones = n.ones_like(sdat)
            print src_name, sig, jys,
            jys = n.average(sdat) / n.average(hmap[_x,_y,z])
            print jys
            print n.sum(_hdatr * sdat) / sig**2 / jys,
            print n.sum(_hdati * sdat) / sig**2 / jys,
            print n.sum(_hdatr**2) / sig**2, 
            print n.sum(_hdati**2) / sig**2 
            print n.sum(_hdatr * sdat) / jys / n.sum(_hdatr**2) + \
                1j * n.sum(_hdati * sdat) / jys / n.sum(_hdati**2)
            print '-' * 70
            csumr += n.sum(_hdatr * sdat) / sig**2 / jys
            csumi += n.sum(_hdati * sdat) / sig**2 / jys
            cwgtr += n.sum(_hdatr**2) / sig**2
            cwgti += n.sum(_hdati**2) / sig**2
        if cwgtr == 0: cwgtr = 1
        if cwgti == 0: cwgti = 1
        alms[L,M] = (csumr/cwgtr + 1j*csumi/cwgti) / n.sqrt(4*n.pi)
        print '**', L,M, alms[L,M]
        print '=' * 70
    _hmapr.from_alm(alms)
    return _hmapr
            
RES = .01
SZ = (int(n.pi/2/RES), int(2*n.pi/RES)+2)
alt,az = n.indices(SZ)
alt = n.pi/2 - alt.astype(n.float) * RES
az = az.astype(n.float) * RES
x,y,z = a.coord.azalt2top((az.flatten(), alt.flatten()))
m = Basemap(projection='ortho', lat_0=90, lon_0=180)
cx,cy = m(az * a.img.rad2deg, alt * a.img.rad2deg)

def plot_hmap(hmap, cnt, max=0, min=-2):
    p.subplot(1, 2, cnt)
    hmap.set_interpol(True)
    if opts.pol == 'xx': data = hmap[x,y,z]
    else: data = hmap[-y,x,z]
    data = n.log10(data**2)
    #max = data.max() * 1.1
    data.shape = az.shape
    m.drawmapboundary()
    m.drawmeridians(n.arange(0,360,30))
    m.drawparallels(n.arange(0,90,10))
    step = (max - min) / 35.
    levels = n.arange(min-step, max+step, step)
    m.contourf(cx, cy, data, levels, linewidths=0)
    #p.colorbar(shrink=.5)
    return m

plot_hmap(hmap, 1)
nhmap = infer_hmap()
plot_hmap(nhmap, 2)
p.show()

for i, filename in enumerate(args):
    freq = freqs[i]
    print 'Plotting', filename, '(%f GHz)' % freq
    ch = int(n.round((freq - uv['sfreq']) / uv['sdf']))
    print 'Using channel %d' % ch
    hmap = a.map.Map(fromfits=filename)
    zen_resp = hmap[0,0,1][0]
    hmap.map.map /= zen_resp
    m = plot_hmap(hmap, i+1)
    p.title('%0.3f GHz' % freq)
    for j in range(3):
        hmap.set_interpol(False)
        for scnt, src_name in enumerate(srclist):
            cat[src_name].update_jys(n.array([freq]))
            jys = cat[src_name].get_jys()[0]
            if jys < 1: continue
            sdat = src_dat[src_name][:,ch]
            valid = n.where(sdat == 0, 0, 1)
            sdat = sdat.compress(valid)
            sx = src_crd[src_name][0].compress(valid)
            sy = src_crd[src_name][1].compress(valid)
            sz = src_crd[src_name][2].compress(valid)
            poly = n.polyfit(sx, sdat, deg=6)
            sm_sdat = n.polyval(poly, sx)
            var = n.sum((sdat - sm_sdat)**2)
            # Account for any beam twist when overlaying sources on calc beam
            if opts.pol == 'xx': rot = aa[0].rot_pol_x
            else: rot = aa[0].rot_pol_y
            sx,sy,sz = n.dot(rot, n.array([sx,sy,sz]))
            if False and j == 0:
                _sx,_sy,_sz = sx.copy(), sy.copy(), sz.copy()
                _sdat = sdat.copy()
                p.subplot(pa1, pa2, 4*pa2 + i+1)
                p.semilogy(_sy, hmap[_sx,_sy,_sz], colors[scnt%len(colors)]+':', label='mdl_'+src_name)
            # Set True for EW mirror-symmetry constraint
            if opts.ew:
                sdat = n.concatenate([sdat, sdat], axis=0)
                sx = n.concatenate([sx, -sx])
                sy = n.concatenate([sy, sy])
                sz = n.concatenate([sz, sz])
            # Set True for NS mirror-symmetry constraint
            if opts.ns:
                sdat = n.concatenate([sdat, sdat], axis=0)
                sx = n.concatenate([sx, sx])
                sy = n.concatenate([sy, -sy])
                sz = n.concatenate([sz, sz])
            # Set True for 180-deg rotational symmetry
            if opts.rot:
                sdat = n.concatenate([sdat, sdat], axis=0)
                sx = n.concatenate([sx, -sx])
                sy = n.concatenate([sy, -sy])
                sz = n.concatenate([sz, sz])
            if False and j == 0:
                paz,palt = a.coord.top2azalt((sx,sy,sz))
                psx,psy = m(paz * a.img.rad2deg, palt * a.img.rad2deg)
                if opts.pol != 'xx': psx,psy = psy,psx
                m.plot(psx, psy, 'k.', markersize=3)
            #swgt = n.ones_like(sdat)
            #swgt = n.ones_like(sdat) * n.log(jys).clip(1,5)
            #swgt = n.ones_like(sdat) * n.sqrt(jys).clip(1,1e3)
            swgt = n.ones_like(sdat) / var * 1e3
            print src_name, var, jys
            #if src_name == 'cyg': swgt *= 100
            # Set True to use overlap points from symmetries above to
            # bootstrap relative gains.
            if False:
                overlap = n.where(hmap.wgt[sx,sy,sz] > 1, 1, 0)
                if n.all(overlap == 0):
                    gain = n.sum(sdat) / n.sum(hmap[sx,sy,sz])
                else:
                    osx = sx.compress(overlap)
                    osy = sy.compress(overlap)
                    osz = sz.compress(overlap)
                    osdat = sdat.compress(overlap)
                    gain = n.sum(osdat) / n.sum(hmap[osx,osy,osz])
                print src_name, gain
            else:
                #gain = jys
                gain = n.sum(sdat) / n.sum(hmap[sx,sy,sz])
            if False and j == 0:
                p.subplot(pa1, pa2, 4*pa2 + i+1)
                p.semilogy(_sy, _sdat/gain, colors[scnt%len(colors)]+'-', label='dat_'+src_name)
            #hmap.put((sx,sy,sz), swgt, (sdat / gain))
            hmap.add((sx,sy,sz), swgt, (sdat / gain))
        if j == 0: plot_hmap(hmap, pa2 + i+1)
        hmap.reset_wgt()
        alm = hmap.map.to_alm(opts.lmax, opts.lmax)
        hmap.map.from_alm(alm)
        if False and j == 2:
            p.subplot(pa1, pa2, 4*pa2 + i+1)
            p.semilogy(_sy, hmap[_sx,_sy,_sz], colors[scnt%len(colors)]+'-.', label='mdl2_'+src_name)
    plot_hmap(hmap, 2*pa2 + i+1)
    if not os.path.exists('new_'+filename): hmap.to_fits('new_'+filename)
    hmap_old = a.map.Map(fromfits=filename)
    hmap_old.map.map /= zen_resp
    hmap_old.map.map = n.abs(hmap.map.map - hmap_old.map.map)
    plot_hmap(hmap_old, 3*pa2 + i+1, max=-2, min=-4)
    
p.show()
