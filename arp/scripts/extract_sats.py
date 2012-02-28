#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, src=True)
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
srclist, cutoff, catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)

tx_width = 2.63671875e-05
afreqs = aa.get_afreqs()

args.sort()
jds = [float('.'.join(filename.split('.')[1:3])) for filename in args]
for cnt, filename in enumerate(args):
  for srcname in cat:
    src = cat[srcname]
    if src.mfreq == .1375:
        print 'Satellite freq for %s not initialized, skipping' % srcname
        continue
    chs = n.argwhere(n.abs(afreqs - src.mfreq) < tx_width)
    aa.select_chans(chs)
    print filename,
    try: times = n.arange(jds[cnt],jds[cnt+1], .0002)
    except(IndexError): times = n.arange(jds[cnt],2*jds[cnt]-jds[cnt-1],.0002)
    skip = True
    for t in times:
        a.phs.AntennaArray.set_jultime(aa, t)
        src.compute(aa)
        if src.alt > 0:
            skip = False
            break
    if skip:
        print ': %s not up, skipping' % src.src_name
        continue
    outfile = filename+'.'+src.src_name
    print '->', outfile
    if os.path.exists(outfile):
        print '   Output file exists.  Skipping...'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(outfile, status='new')
    #sfreq = uvi['sfreq']+chs[0]*uvi['sdf']
    override = { #'nchan':3,
        #'sfreq':sfreq, 
        #'freq':sfreq, 
        #'freqs':(4,)+(3, sfreq, uvi['sdf']) * 4,
    }

    curtime = None
    def mfunc(uv, p, d, f):
        global curtime
        crd,t,(i,j) = p
        if t != curtime:
            aa.set_jultime(t)
            src.compute(aa)
            curtime = t
        if src.alt < 0: return p, None, None
        if True:
            d = n.average(d.take(chs))
            f = f.take(chs)
            f = n.any(f)
            return p, n.array([d]), n.array([f])
        else: return p, d, f

    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='EXTRACT_SAT: extracting %s\n' % src.src_name)
    del(uvo)
        

