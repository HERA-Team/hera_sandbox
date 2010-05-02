#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse, re
import glob, os
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
opts,args = o.parse_args(sys.argv[1:])

assert(opts.pol in ['xx','yy'])

re_freq = re.compile(r'_(\d+)\.hmap$')
freqs = n.array([float(re_freq.search(f).groups()[0]) / 1e3 for f in args])

colors = 'kbrgmc'
srcfiles = glob.glob('*.bm_*m')
srclist = [f.split('.')[-1][3:-1] for f in srcfiles]
cat = a.src.get_catalog(srclist)
#cat = a.cal.get_catalog(opts.cal, srclist)
uv = a.miriad.UV(srcfiles[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

print 'Gathering source data'
src_dat, src_crd = {}, {}
for src_name, filename in zip(srclist, srcfiles):
    print '    ', filename
    uv = a.miriad.UV(filename)
    src = cat[src_name]
    dat = []
    x,y,z = [],[],[]
    for (crd,t,bl),d,f in uv.all(raw=True):
        #aa.set_jultime(t)
        #src.compute(aa)
        a.phs.ArrayLocation.set_jultime(aa, t)
        a.phs.RadioFixedBody.compute(src, aa)
        if src.alt < 25 * a.img.deg2rad: continue
        xi,yi,zi = src.get_crds('top', ncrd=3)
        x.append(xi); y.append(yi); z.append(zi)
        dat.append(n.where(f,0,n.abs(d)))
    src_dat[src_name] = n.array(dat)
    src_crd[src_name] = (n.array(x), n.array(y), n.array(z))

pa1 = 4
pa2 = len(args)
RES = .01
SZ = (int(n.pi/2/RES), int(2*n.pi/RES)+2)
alt,az = n.indices(SZ)
alt = n.pi/2 - alt.astype(n.float) * RES
az = az.astype(n.float) * RES
x,y,z = a.coord.azalt2top((az.flatten(), alt.flatten()))
m = Basemap(projection='ortho', lat_0=90, lon_0=180)
cx,cy = m(az * a.img.rad2deg, alt * a.img.rad2deg)

def plot_hmap(hmap, cnt, max=0, min=-2):
    p.subplot(pa1, pa2, cnt)
    hmap.set_interpol(True)
    if opts.pol == 'xx': data = hmap[x,y,z]
    else: data = hmap[-y,x,z]
    data = n.log10(data**2)
    data.shape = az.shape
    m.drawmapboundary()
    m.drawmeridians(n.arange(0,360,30))
    m.drawparallels(n.arange(0,90,10))
    step = (max - min) / 35.
    levels = n.arange(min-step, max+step, step)
    m.contourf(cx, cy, data, levels)
    #p.colorbar(shrink=.5)
    return m

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
    for j in range(1):
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
                #gain = n.sum(sdat) / n.sum(hmap[sx,sy,sz])
                gwgt = hmap[sx,sy,sz]**2
                gain = n.sum(sdat*gwgt) / n.sum(gwgt**2)
            print src_name, freq, var, jys, gain#, '(%f)' % (n.sum(sdat)/n.sum(gwgt))
            if False and j == 0:
                p.subplot(pa1, pa2, 4*pa2 + i+1)
                p.semilogy(_sy, _sdat/gain, colors[scnt%len(colors)]+'-', label='dat_'+src_name)
            #hmap.put((sx,sy,sz), swgt, (sdat / gain))
            #hmap.add((sx,sy,sz), swgt, (sdat / gain))
            hmap.add((sx,sy,sz), swgt, n.sqrt(sdat / gain))
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
