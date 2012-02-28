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
    curtime = None
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
            aa.set_jultime(t)
            lst = aa.sidereal_time()
            ubin,vbin,lstbin = C.pspec.bin2uv(C.pspec.uv2bin(0,0,lst))
            zen = a.phs.RadioFixedBody(lstbin, aa.lat)
            zen.compute(aa)
            # not bothering with beam weighting for now
            curtime = t
            sys.stdout.write('.'); sys.stdout.flush()
        bl = a.miriad.ij2bl(i,j)
        u,v,w = aa.gen_uvw(i,j, src=zen) # phased to zenith of bin center
        # Grid by highest fq = fastest uv movement
        # Will need to rebin by subband to reap full sensitivity...
        u,v = u.flatten()[-1],v.flatten()[-1]
        umag = n.sqrt(u**2 + v**2).flatten()
        #if umag > 50: continue
        valid = n.logical_not(f).astype(n.float)
        if n.average(valid) < .5: continue
        # Apply one delay correction; will decohere widefield sky for high w's...
        d = aa.phs2src(d, zen, i, j)
        # Make sure this bandpass is smooth all the way to the edge of valid...
        d /= aa.passband(i,j)
        d *= valid
        # Pick a 1/2 of the uv plane
        if u < 0: u,v,d = -u,-v,n.conj(d)
        bin = C.pspec.uv2bin(u.flatten(), v.flatten(), lst)
        #print (i,j), umag, bin
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
