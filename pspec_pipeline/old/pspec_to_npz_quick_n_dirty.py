#! /usr/bin/env python
import aipy as a, numpy as n, capo as C
import sys, optparse, os

o = optparse.OptionParser()
o.add_option('--slots', dest='slots', type='int', default=5,
    help='The number of slots within a uv-lst bin allocated for holding spectra prior to squaring by cross-multiplying slots.')
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True)
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
freqs = uv['sfreq'] + n.arange(uv['nchan']) * uv['sdf']
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

for filename in args:
    print 'Reading', filename
    ofile = filename+'.npz'
    if os.path.exists(ofile):
        print ofile, 'exists.  Skipping...'
        continue
    dat, wgt = {}, {}
    cnt, bls = {}, {}
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        valid = n.logical_not(f).astype(n.float)
        bl = a.miriad.ij2bl(i,j)
        # Grid by highest fq = fastest uv movement
        # Will need to rebin by subband to reap full sensitivity...
        bin = int(uv['bin'])
        if not dat.has_key(bin):
            dat[bin], wgt[bin], cnt[bin] = [], [], 0
            for slot in xrange(opts.slots):
                dat[bin].append(n.zeros_like(d))
                wgt[bin].append(n.zeros_like(valid))
        #print i,j, bin, cnt[bin] % opts.slots
        # May want to assign baselines to specific slots to cross-mults don't
        # have crosstalk in them...
        dat[bin][cnt[bin] % opts.slots] += d
        wgt[bin][cnt[bin] % opts.slots] += valid
        cnt[bin] += 1
        if not bls.has_key(bin): bls[bin] = {}
        bls[bin][bl] = None

    # Can save space by throwing out low SNR bins (but be careful
    # about LST bin boundaries...
    bins = [b for b in cnt if cnt[b] >= opts.slots]
    #bins = [b for b in cnt if cnt[b] >= 10*opts.slots]
    print len(bins)
    dat = n.array([dat[bin] for bin in bins])
    wgt = n.array([wgt[bin] for bin in bins])
    # Store list of 5 contributing baselines to each bin
    for bin in bls: bls[bin] = n.array(bls[bin].keys() + [0]*5)[:5]
    bls = n.array([bls[bin] for bin in bins])
    print '    Writing', ofile
    # should stick some more info in here: LSTBIN, UVBIN
    n.savez(ofile, bins=n.array(bins), dat=dat, wgt=wgt, freqs=freqs, bls=bls)
