#! /usr/bin/env python
import aipy as a, numpy as n
import optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('-s', '--src', dest='src', type='str', 
    help='Source to use for calibration.')
o.add_option('--cat', dest='cat', type='str', default='helm,misc',
    help='A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aipy._src.  Some catalogs may require a separate data file to be downloaded and installed.')
o.add_option('-d','--deg', dest='deg', type='int', default=8,
    help='Degree of polynomial to fit to bandpass function.')
o.add_option('-q','--quiet', dest='quiet', action='store_true',
    help='Do not plot anything (be visually quiet).')
o.add_option('--batch', dest='batch', action='store_true',
    help='Operate in batch mode.  Saves each plot to "spec_<srcname>.png".')
opts,args = o.parse_args(sys.argv[1:])

if not opts.quiet:
    if opts.batch: import matplotlib; matplotlib.use('Agg')
    import pylab as p

if opts.src != None:
    srclist,cutoff,catalogs, = a.scripting.parse_srcs(opts.src, opts.cat)
    calsrc = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)
    calsrc = calsrc.values()[0]
else:
    calsrc = None
srclist = [f.split('_')[0] for f in args]
if not calsrc is None: assert(calsrc.src_name in srclist)
srclist,cutoff,catalogs, = a.scripting.parse_srcs(','.join(srclist), opts.cat)
print srclist
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)

filetypes = {}
for src in [calsrc] + cat.values():
    if src is None: continue
    for filename in [f for f in args if f.startswith(src.src_name)]:
        print src.src_name, filename, filename[len(src.src_name):]
        filetypes[filename[len(src.src_name):]] = None
filetypes = filetypes.keys()
print filetypes
filetypes.sort()
srccolors = 'krbg'
colors = {}
for cnt,f in enumerate(filetypes): colors[f] = srccolors[cnt % len(srccolors)]
print colors.keys()

bp_cal = {}
for src in [calsrc] + cat.values():
    if src is None: continue
    srcfiles = [f for f in args if f.startswith(src.src_name)]
    for cnt, filename in enumerate(srcfiles):
        filetype = filename[len(src.src_name):]
        color = colors[filetype]
        print color
        print 'Reading', filename
        _f = open(filename)
        f = n.load(_f)
        spec,afreqs = f['spec'].flatten(), f['freq'].flatten()
        _f.close()
        #dspec = spec - a.rfi.remove_spikes(spec, order=8)
        #sig = n.std(dspec)
        if True:
            valid = n.where(n.logical_and(afreqs > .120, afreqs < 0.170), 1, 0)
            spec = spec.compress(valid)
            afreqs = afreqs.compress(valid)
        p.semilogy(afreqs, spec, label='unmodified spec')
        src.update_jys(afreqs)
        p.semilogy(afreqs, src.jys, label='actual src')
#        print src.jys
#        print afreqs
        bp = n.sqrt(spec / src.jys)
        p.semilogy(afreqs, bp, label='bandpass = spec/src.jys')
        bp_poly = n.polyfit(afreqs, bp, deg=opts.deg)
        if not calsrc is None and src.src_name == calsrc.src_name:
            print 'Calibrating to', src.src_name
            bp_cal[filetype] = bp_poly
        if opts.src is None:
            bp_fit = 1.
        else: bp_fit = n.polyval(bp_cal[filetype], afreqs).clip(.1,10)**2
        spec /= bp_fit
        
        src_poly = n.polyfit(n.log10(afreqs/src.mfreq), n.log10(spec), deg=1)
        n.set_printoptions(threshold=n.nan)
        print 'bp =', list(bp_poly)
        print "'%s':" % src.src_name + "{ 'jys':10**%f, 'index':  %f , }," % (src_poly[-1], src_poly[-2])
        print 'RMS residual:', n.sqrt(n.average((spec - 10**n.polyval(src_poly, n.log10(afreqs/src.mfreq)))**2))
        if n.all(spec <= 0): continue

        if not opts.quiet:
            p.semilogy(afreqs, spec, color+'.', label='Measured')
            p.semilogy(afreqs, 10**n.polyval(src_poly, n.log10(afreqs/src.mfreq)), color+'-', 
                label='Fit Power Law')
            p.legend()
            p.show()
    if not opts.quiet:
        #p.loglog(afreqs, src.jys, color+':', label='%f, %s' % (src._jys, str(src.index)))
        p.xticks(n.arange(.1,.2,.02), ['100','120','140','160','180'])
        p.xlim(afreqs[0], afreqs[-1])
        p.ylim(3,3e3)
        p.grid()
        p.title(src.src_name)
        #p.legend()
        if opts.batch:
            fig = p.gcf()
            fig.set_size_inches(3.,2.25)
            fig.subplots_adjust(left=.125, right=.97, bottom=.1, top=.875)
            outfile = 'spec_%s.png' % (src.src_name)
            print 'Saving to', outfile
            outfile = open(outfile,'w')
            p.savefig(outfile)
            outfile.close()
        #else: p.show()
        #p.clf()
p.show()
