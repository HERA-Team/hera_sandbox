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
srcfiles = glob.glob('*+4[0-4]*.bm_*m')
srclist = [f.split('.')[-1][3:-1] for f in srcfiles]
cat = a.src.get_catalog(srclist)
#cat = a.cal.get_catalog(opts.cal, srclist)
uv = a.miriad.UV(srcfiles[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
afreqs = aa.get_afreqs()
del(uv)

print 'Gathering source data'
profile_dat = 0
profile_wgt = 0
src_dat, src_crd = {}, {}
for src_name, filename in zip(srclist, srcfiles):
    print '    ', filename
    uv = a.miriad.UV(filename)
    src = cat[src_name]
    dat = []
    x,y,z = [],[],[]
    for (crd,t,bl),d,f in uv.all(raw=True):
        a.phs.ArrayLocation.set_jultime(aa, t)
        a.phs.RadioFixedBody.compute(src, aa)
        if src.alt < 25 * a.img.deg2rad: continue
        xi,yi,zi = src.get_crds('top', ncrd=3)
        x.append(xi); y.append(yi); z.append(zi)
        dat.append(n.where(f,0,n.abs(d)))
    src_dat[src_name] = n.array(dat) 
    src_crd[src_name] = (n.array(x), n.array(y), n.array(z))
    src.update_jys(afreqs)
    CH = 8
    jys = src.get_jys()[CH]
    if jys < 1: continue
    d = n.abs(src_dat[src_name][:,CH])
    valid = n.where(d > 0, 1, 0)
    vx = src_crd[src_name][0].compress(valid)
    vd = d.compress(valid) / jys
    poly = n.polyfit(vx, vd, deg=4)
    sm_d = n.polyval(poly, vx)
    var = n.sum((vd - sm_d)**2)
    print src_name, var
    hdat, bins = n.histogram(vx, weights=vd, range=(-1,1), bins=500)
    hwgt, bins = n.histogram(vx, range=(-1,1), bins=500)
    profile_dat += hdat / var
    profile_wgt += hwgt / var
    hdat /= hwgt
    p.plot(.5*(bins[:-1]+bins[1:]), n.log10(hdat), '.')
    #p.plot(vx, n.log10(sm_d), '^')
profile_dat /= profile_wgt
p.plot(.5*(bins[:-1]+bins[1:]), n.log10(profile_dat), '^')
p.show()
