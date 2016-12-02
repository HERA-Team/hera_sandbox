#! /usr/bin/env python
"""
Removes crosstalk from PAPER data by modeling it over the course of a day
as a static per-channel cross-coupling that rises and falls uniformly across
the band with changing input amplitude.  Steps for crosstalk removal are:
run "xtalk3.py -o" on 1 day of UV files (which should be flagged for RFI, and
if possible, have a sky model removed) to 
generate *.xtalk.npz files.  Then run "xtalk3.py -r" on same UV files to 
reprocess *.xtalk.npz files, separating them into static shape/uniform gain 
terms that are stored in *.xtalk.rep.npz files.  Finally, run "xtalk3.py -i" on 
UV files from the same JD (but which need not have RFI flagged or a sky model
removed) to use the *.xtalk.rep.npz model to remove crosstalk from the visibility 
data.

Author: Aaron Parsons
"""

import aipy as a, numpy as np, sys, os, optparse
import capo.omni as omni

o = optparse.OptionParser()
o.set_usage('xtalk3.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o,cal=True)
o.add_option('-i', '--infile', dest='infile', action='store_true',
    help='Apply xtalk calibrations generated with the -o option.')
o.add_option('-o', '--outfile', dest='outfile', action='store_true',
    help='Rather than apply the calibrations to the data, store them in a file (named by JD) to apply to a different file with the same JD.')
o.add_option('-r', '--reprocess', dest='reprocess', action='store_true',
    help='Reprocess xtalk files for the specified UV files to generate an average crosstalk profile and a time-dependent gain for scaling that profile.  The results will be saved to *.xtalk.rep.npz files which will have precedence over *.xtalk.npz files.')
o.add_option('-c', '--chan', dest='chan', default='160_720',
    help='Channels range to use for normalization when reprocessing.')
o.add_option('--omnimdl', dest='omnimdl', action='store',
    help='Subtract model of visibility given in the omnical model generation')
o.add_option('--verbose', action='store_true',
    help='Be verbose')
opts,args = o.parse_args(sys.argv[1:])

assert(not opts.infile or not opts.outfile)

if opts.reprocess:
    chans = map(int, opts.chan.split('_'))
    xsum, cnt, gain, times = {}, {}, {}, []
    for filename in args:
        uv = a.miriad.UV(filename)
        (uvw,jd,(i,j)),d,f = uv.read(raw=True)
        xfile = '%f.xtalk.npz' % jd
        print xfile
        times.append(jd)
        if not os.path.exists(xfile):
            print xfile, 'does not exist.  Skipping...'
            continue
        xtalk = np.load(xfile)
        for bl in xtalk.files:
            print bl
#            import IPython; IPython.embed()
            dat = np.array(xtalk[bl])
            adat = np.ma.masked_equal(np.abs(dat), 0)
            if not gain.has_key(bl): gain[bl] = []
            gain[bl].append(np.average(adat[chans[0]:chans[1]]))
            dat /= gain[bl][-1]
            xsum[bl] = xsum.get(bl, 0) + dat
            cnt[bl] = cnt.get(bl, 0) + np.logical_not(adat.mask).astype(np.int)
    for bl in xsum: xsum[bl] /= np.where(cnt[bl] == 0, 1, cnt[bl])
    for c, jd in enumerate(times):
        repfile = '%f.xtalk.rep.npz' % jd
        xtalk = {}
        for bl in xsum: xtalk[bl] = gain[bl][c] * xsum[bl]
        print 'Writing', repfile
        np.savez(repfile, **xtalk)
    import sys; sys.exit(0)

guess, cnt, xtalk = {}, {}, {}
for filename in args:
    print filename,'->',filename+'x'
    if not opts.outfile and os.path.exists(filename+'x'):
        print filename+'x', 'exists.  Skipping...'
        continue
    uv = a.miriad.UV(filename)
    uv.select('auto',0, 0, include=False)
    (uvw,jd,(i,j)),d,f = uv.read(raw=True)
    uv.rewind()
    if opts.infile:
        xfile = '%f.xtalk.rep.npz' % jd
        if not os.path.exists(xfile): xfile = '%f.xtalk.npz' % jd
        if not os.path.exists(xfile):
            print xfile, 'does not exist.  Skipping...'
            continue
        print '    using', xfile
        xtalk = np.load(xfile)
    else:
        guess, cnt, xtalk = {}, {}, {}
        if opts.omnimdl:
            #if pass in omnical file, load in the mdl
            omnifile = opts.omnimdl % '.'.join(filename.split('/')[-1].split('.')[0:4])
            m,_,mdl,_ = omni.from_npz(omnifile)
            #make info object so that we can get a mapping of mdl bls to all bls. 
            aa = a.cal.get_aa(opts.cal, np.array([.150]))
            info = omni.aa_to_info(aa)
            redmapping = {}
            for k in mdl[a.miriad.pol2str[uv['pol']]].keys():
                for gp in info.get_reds():
                    if k in gp:
                        for kk in gp : redmapping[kk] = k
        for (uvw,t,(i,j)),d,f in uv.all(raw=True):
            ti = np.where(m['jds']==t)[0]
            bl = str(a.pol.ijp2blp(i,j,uv['pol']))
            if not guess.has_key(bl): guess[bl], cnt[bl] = 0, 0
            if (i,j) in [(97,80), (80,97), (72,96), (96,72), (43, 88), (88,43)]:
                if opts.verbose: print 'No Model Visibility for {0}'.format((i,j))
                ml = np.zeros_like(d)
            else:
                try:
                    ml = mdl[a.miriad.pol2str[uv['pol']]][redmapping[(i,j)]][ti].flatten()
                except KeyError:
                    ml = mdl[a.miriad.pol2str[uv['pol']]][redmapping[(j,i)]][ti].flatten()
            guess[bl] += np.where(f, 0, d)-np.where(f,0,ml)
            cnt[bl] += np.logical_not(f)
        del(uv)
        for bl in guess: xtalk[bl] = guess[bl] / np.clip(cnt[bl], 1, np.Inf)
    if opts.outfile:
        xfile = '%f.xtalk.npz' % jd
        print 'Writing', xfile
        np.savez(xfile, **xtalk)
    else:
        def mfunc(uv, p, d, f):
            uvw,t,(i,j) = p
            bl = str(a.pol.ijp2blp(i,j,uv['pol']))
            try: return p, d - xtalk[bl], f
            except(KeyError): return p, d, f
        uvi = a.miriad.UV(filename)
        uvo = a.miriad.UV(filename+'x', status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc=mfunc, append2hist='XTALK:'+' '.join(sys.argv)+'\n', raw=True)
