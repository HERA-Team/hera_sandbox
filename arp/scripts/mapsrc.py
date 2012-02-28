#! /usr/bin/env python
import aipy as a, numpy as n, sys, optparse, re
#import pylab as p

re_filename = re.compile(r'.*_c(\d+)_(\d+)_.*')

o = optparse.OptionParser()
o.set_description(__doc__)
a.scripting.add_standard_options(o, src=True)
o.add_option('--sfreq', dest='sfreq', type='float', default=0.1001953125)
o.add_option('--sdf', dest='sdf', type='float', default=0.000390625)
o.add_option('--nside', dest='nside', type='str',
    help="Manually set NSIDE (possibly degrading map) to a power of 2.")
opts,args = o.parse_args(sys.argv[1:])

srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.src.get_catalog(srclist, cutoff, catalogs)

def ch2freq(chan, sfreq, sdf): return (sfreq + chan*sdf)

srccrd = {}
for src in cat: srccrd[src] = a.coord.radec2eq((cat[src]._ra, cat[src]._dec))

if opts.nside is None: opts.nside = [None]
else: opts.nside = map(int, opts.nside.split(','))

fqs,maps = {}, {}
for map in args:
    print 'Reading', map,
    ch1, ch2 = re_filename.match(map).groups()
    ch = .5 * (float(ch1) + float(ch2))
    fq = ch2freq(ch, sfreq=opts.sfreq, sdf=opts.sdf)
    fqs[map] = fq
    print 'at %0.4f GHz' % fq
    m = a.map.Map(fromfits=map)
    m.set_interpol(False)
    maps[map] = m

for nside in opts.nside:
    if nside is None: nside = m.nside()
    print 'NSIDE:', nside
    srcspec,_fqs = {},[]
    for map in maps:
        _fqs.append(fqs[map])
        if nside is None: m = maps[map]
        else:
            m = a.map.Map(nside=int(nside))
            m.set_interpol(False)
            m.from_map(maps[map])
            #m.set_interpol(True)
            #m.map.map *= (maps[map].nside() / m.nside())**2
        for src in cat: srcspec[src] = srcspec.get(src,[]) + [m[tuple(srccrd[src])]]
    _fqs = n.array(_fqs)
    inds = n.argsort(_fqs)
    _fqs = _fqs[inds]
    for src in cat:
        #print "'%s' : {'jys':10**%f, 'index':[%f], }," % (src
        srcspec[src] = n.array(srcspec[src])[inds]
        outfile = '%s_spec_nside%d.npz' % (src,nside)
        print 'Writing spectrum to', outfile
        n.savez(outfile, spec=srcspec[src], freq=_fqs)
#        p.loglog(_fqs, srcspec[src], '-', label=src)
#p.xticks(n.arange(.1,.2,.02), ['100','120','140','160','180'])
#p.xlim(.115,.185)
#p.ylim(3,3e3)
#p.grid()
#p.legend()
#p.show()
