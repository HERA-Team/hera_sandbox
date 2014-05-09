#!/usr/bin/python
import aipy as a, numpy as n
import optparse, sys, os

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
o.set_usage('xrfi_simple.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-c', '--chan', dest='chan',
    help='Manually flag channels before xrfi processing.  Options are "<chan1 #>,..." (a list of channels), or "<chan1 #>_<chan2 #>" (a range of channels).  Default is None.')
o.add_option('-n', '--nsig', dest='nsig', default=2., type='float',
    help='Number of standard deviations above mean to flag if neither --dt nor --df are specified.  Default 2.')
o.add_option('--df', dest='df', type='float',
    help='Number of standard deviations above mean to flag, after taking derivative of frequency axis')
o.add_option('--dt', dest='dt', type='float',
    help='Number of standard deviations above mean to flag, after taking derivative of time axis')
o.add_option('--combine', dest='combine', action='store_true',
    help='Use the same mask for all baselines/pols (and use thresh to decide how many concidences it takes to flag all data.')
o.add_option('--to_npz',
    help='Instead of applying mask to data, store it as npz of this name.  May only be used along with --combine.')
o.add_option('--from_npz',
    help='Apply mask to data from this npz file (generated with --to_npz).  May only be used along with --combine.')
o.add_option('-t', '--thresh', dest='thresh', default=1, type='int',
    help='Number of flagging coincidences (baselines/pols) required to flag a time/chan.')
opts,args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
if not opts.chan is None:
    chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
else:
    chans = []
    opts.chan = 'None'
del(uv)
if opts.to_npz or opts.from_npz: assert(opts.combine)


for uvfile in args:
    uvofile = uvfile+'R'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    if opts.from_npz:
        print '    Reading flags from', opts.from_npz
        m = n.load(opts.from_npz)
        mask = {'xx':{257:{}}} # Just use dummy values here to mimic structure of mask dictionary
        for cnt,t in enumerate(m['times']):
            mask['xx'][257][t] = m[str(cnt)]
    else:
        uvi = a.miriad.UV(uvfile)
        a.scripting.uv_selector(uvi, opts.ant)
        # Gather all data and each time step
        data,mask,times = {}, {}, []
        for (uvw,t,(i,j)), d, f in uvi.all(raw=True):
            if len(times) == 0 or times[-1] != t: times.append(t)
            blp = a.pol.ijp2blp(i,j,uvi['pol'])
            if not blp in data:
                data[blp] = {}
                mask[blp] = {}
            # Manually flag channels
            f[chans] = 1
            mask[blp][t] = f
            data[blp][t] = d

        # Generate statistical mask
        for blp in data:
            i,j,p = a.pol.blp2ijp(blp)
            data_times = data[blp].keys()
            data_times.sort()
            d = n.array([data[blp][t] for t in data_times])
            m = n.array([mask[blp][t] for t in data_times])
            if opts.df != None:
                ddf = d[:,1:-1] - .5 * (d[:,:-2] + d[:,2:])
                ddf2 = n.abs(ddf)**2
                sig = n.sqrt(n.median(ddf2, axis=1))
                sig.shape = (sig.size,1)
                m[:,0] |= 1; m[:,-1] |= 1
                m[:,1:-1] |= n.where(ddf2/sig**2 > opts.df**2, 1, 0)
            if opts.dt != None:
                ddt = d[1:-1,:] - .5 * (d[:-2,:] + d[2:,:])
                ddt2 = n.abs(ddt)**2
                sig = n.sqrt(n.median(ddt2, axis=0))
                sig.shape = (1,sig.size)
                m[0,:] |= 1; m[-1,:] |= 1
                m[1:-1,:] |= n.where(ddt2/sig**2 > opts.dt**2, 1, 0)
            if opts.df == None and opts.dt == None:
                ad = n.abs(d)
                med = n.median(ad)
                sig = n.sqrt(n.median(n.abs(ad-med)**2))
                m |= n.where(ad > med + opts.nsig * sig, 1, 0)
            for i, t in enumerate(data_times):
                mask[blp][t] |= m[i]
        if opts.combine:
            new_mask = {}
            for blp in mask:
                new_mask[t] = new_mask.get(t,0)+mask[blp][t].astype(n.int)
            for t in new_mask:
                m = n.where(new_mask[t] >= opts.thresh, 1, 0)
                for blp in mask:
                    mask[blp][t] = m
        del(uvi)

    if opts.to_npz:
        print '    Writing flags to', opts.to_npz
        m = {}
        for blp in mask.keys():
            times = n.array(mask[blp].keys())
            for i,t in enumerate(mask[blp].keys()):
                m[str(i)] = mask[blp][t]
        m['times'] = times
        #_m = mask.values()[0].values()[0]
        #times = n.array(_m.keys())
        #for cnt,t in enumerate(times): m[str(cnt)] = _m[t]
        #m['times'] = times
        n.savez(opts.to_npz, **m)
    else:
        # Generate a pipe for applying the mask to data as it comes in.
        def rfi_mfunc(uv, preamble, data, flags):
            uvw, t, (i,j) = preamble
            blp = a.pol.ijp2blp(i,j,uv['pol'])
            try:
                if opts.from_npz:
                    m = mask['xx'][257][t]
                else:
                    m = mask[blp][t]
            except(KeyError):
                m = n.ones_like(flags)
            return preamble, n.where(m, 0, data), m

        uvi = a.miriad.UV(uvfile)
        uvo = a.miriad.UV(uvofile, status='new')
        uvo.init_from_uv(uvi)
        uvo.pipe(uvi, mfunc=rfi_mfunc, raw=True, append2hist=' '.join(sys.argv)+'\n')
