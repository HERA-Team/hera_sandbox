#!/usr/bin/env python
"""
A script for filtering using a delay/delay-rate transform.  If a source
is specified, will remove/extract that source.  If none is specified,
will filter/extract in absolute terms.
"""

import aipy as a, numpy as n, os, sys, optparse, ephem
import capo as C

o = optparse.OptionParser()
o.set_usage(__file__ + ' [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True)
o.add_option('-s','--sum', dest='sumfiles', action='store_true',
    help='Sum uv files into a single output npz file')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
freqs = aa.get_afreqs()

MIN_CH,MAX_CH = 60,804
aa.select_chans(n.arange(MIN_CH, MAX_CH))
lstbins = n.arange(0,2*n.pi, C.pspec.LST_RES)
if opts.sumfiles:
    ofile = 'summed__cross.npz'
    if os.path.exists(ofile):
        print ofile, 'exists.  Skipping...'
        sys.exit(0)
    dat = {'freqs':freqs[MIN_CH:MAX_CH]}

for filename in args:
    print 'Reading', filename
    if not opts.sumfiles:
        ofile = '%s__cross.npz' % (filename)
        if os.path.exists(ofile):
            print ofile, 'exists.  Skipping...'
            continue
        dat = {'freqs':freqs[MIN_CH:MAX_CH]}
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    times = []
    src = None
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if len(times) == 0 or t != times[-1]:
            aa.set_jultime(t)
            # need to bin in lst to accumulate coherently...
            lst = lstbins[n.argmin(n.abs(lstbins - aa.sidereal_time()))]
            src = a.phs.RadioFixedBody(lst, aa.lat)
            src.compute(aa)
            times.append(t)
            #bm_wgt = aa[0].bm_response(src.get_crds('top'), pol='y')**2
            #bm_wgt = bm_wgt.flatten()
        bl = a.miriad.ij2bl(i,j)
        x,y,z = aa.get_baseline(i,j, src=src)
        u150,v150,w150 = .150*x, .150*y, .150*z
        # mask out overly short baselines adversely affected by xtalk removal
        if u150 < 12.5: continue
        # ignoring w here, which must be carefully considered.  over FoV, baselines w/ same
        # u,v but diff w will have sources away from phase center come in at slightly diff delays
        # Is this ok?
        bin = C.pspec.uv2bin(u150,v150,lst)
        d = d[MIN_CH:MAX_CH]
        w = n.logical_not(f[MIN_CH:MAX_CH]).astype(n.float)
        if n.average(w) < .5: continue
        d = aa.phs2src(d, src, i, j)
        d /= aa.passband(i,j)
        d *= w
        # Weight by square of primary beam response in the direction we're projecting for.
        #d *= bm_wgt
        #w = bm_wgt**2 * val
        sumkey, wgtkey = 'sum_%d' % (bin), 'wgt_%d' % (bin)
        dat[sumkey] = dat.get(sumkey,0) + d
        dat[wgtkey] = dat.get(wgtkey,0) + w

    if not opts.sumfiles:
        dat['times'] = n.array(times)
        n.savez(ofile, **dat)

if opts.sumfiles:
    dat['times'] = n.array(times)
    n.savez(ofile, **dat)

