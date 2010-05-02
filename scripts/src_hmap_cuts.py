#! /usr/bin/env python
import aipy as a, pylab as p, numpy as n, sys, optparse, re
import glob, os

colors = 'kbrgmc'
o = optparse.OptionParser()
o.set_usage('src_hmap_cuts.py [options] *.hmap')
a.scripting.add_standard_options(o, cal=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

assert(opts.pol in ['xx','yy'])

re_freq = re.compile(r'_(\d+)\.hmap$')
freqs = n.array([float(re_freq.search(f).groups()[0]) / 1e3 for f in args])

#srcfiles = [glob.glob('*.bm_%sm' % (s))[0] for s in srclist]
srcfiles = glob.glob('*.bm_J0[689]*m')
srclist = [f.split('.')[-1][3:-1] for f in srcfiles]
src_dat = {}
src_crd = {}
cat = a.src.get_catalog(srclist)
uv = a.miriad.UV(srcfiles[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

print 'Gathering source data'
for src_name, filename in zip(srclist, srcfiles):
    print '    ', filename
    uv = a.miriad.UV(filename)
    src = cat[src_name]
    dat = []
    x,y,z = [],[],[]
    for (crd,t,bl),d,f in uv.all(raw=True):
        aa.set_jultime(t)
        src.compute(aa)
        if src.alt < 25 * a.img.deg2rad: continue
        xi,yi,zi = src.get_crds('top', ncrd=3)
        x.append(xi); y.append(yi); z.append(zi)
        dat.append(n.where(f,0,n.abs(d)))
    src_dat[src_name] = n.sqrt(n.array(dat))
    src_crd[src_name] = (n.array(x), n.array(y), n.array(z))

# For 150 MHz
#src_gains = {
#    'cyg': 1.0,
#    'vir': 1.00236295962,
#    'cas': 0.986401136049,
#    'crab': 0.914414228674,
#    'her': 1.20956302049,
#    'hyd': 0.252997732775,
#}
src_gains = {
    'cyg': 1.0474487692,
    'vir': 1.04110813432,
    'cas': 1.1247132013,
    'crab': 1.09105007454,
}

p1 = int(n.ceil(n.sqrt(len(args))))
p2 = int(n.ceil(len(args)/float(p1)))
for i, filename in enumerate(args):
    freq = freqs[i]
    print 'Plotting', filename, '(%f GHz)' % freq
    ch = int(n.round((freq - uv['sfreq']) / uv['sdf']))
    print 'Using channel %d' % ch
    hmap = a.map.Map(fromfits=filename)
    zen_resp = hmap[0,0,1][0]
    hmap.map.map /= zen_resp
    hmap2 = a.map.Map(fromfits=filename)
    hmap2.map.map /= zen_resp
    hmap.set_interpol(True)
    p.subplot(p1,p2,i+1)
    p.title(filename + ' (%0.3f GHz)' % freq)
    for j, src_name in enumerate(srclist):
        cat[src_name].update_jys(n.array([freq]))
        sdat = src_dat[src_name][:,ch]
        valid = n.where(sdat == 0, 0, 1)
        sdat = sdat.compress(valid)
        sx = src_crd[src_name][0].compress(valid)
        sy = src_crd[src_name][1].compress(valid)
        sz = src_crd[src_name][2].compress(valid)
        _sx = sx.copy()
        # Account for any beam twist when overlaying sources on calc beam
        if opts.pol == 'xx': rot = aa[0].rot_pol_x
        else: rot = aa[0].rot_pol_y
        sx,sy,sz = n.dot(rot, n.array([sx,sy,sz]))
        hdat = hmap[sx,sy,sz]
        gain = n.sqrt(n.average(sdat**2/hdat**2))
        print src_name, gain
        #gain = src_gains[src_name]
        gain = 1
        p.semilogy(_sx,sdat/gain, colors[j]+'-', label=src_name)
        p.semilogy(_sx,hdat, colors[j]+':', label='mdl_'+src_name)
    p.ylim(.05, 2)
#p.legend()
p.show()
