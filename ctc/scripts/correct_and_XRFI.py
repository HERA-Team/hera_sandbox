#!/usr/bin/env python
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

########
# XRFI #
########

for uvfile in args:
    uvofile = uvfile+'cR'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue

    uvi = a.miriad.UV(uvfile)
    a.scripting.uv_selector(uvi, opts.ant)
    # Gather all data and each time step
    data,mask,times = {}, {}, []
    for (uvw,t,(i,j)), d, f in uvi.all(raw=True):

        if len(times) == 0 or times[-1] != t:
            times.append(t)

        blp = a.pol.ijp2blp(i, j, uvi['pol'])
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
    new_mask = {}
    for blp in mask.keys():
        for t in mask[blp].keys():
            new_mask[t] = new_mask.get(t,0)+mask[blp][t].astype(n.int)
        del(mask[blp])
    for t in new_mask.keys():
        m = n.where(new_mask[t] >= opts.thresh, 1, 0)
        mask[t] = m
    del(uvi)

    ###########
    # CORRECT #
    ###########

    def same_col(i, j):
        return i/8 == j/8 and i != j

    rewire = {}
    #rewire = {0:1, 1:24, 2:45, 3:100}

    aa = a.phs.ArrayLocation(('-30:43:17.5','21:25:41.9'))
    nints = 0
    curtime = None
    def mfunc(uv, preamble, data, flags):
        global curtime
        global nints
        uvw, t, (i,j) = preamble
        blp = a.pol.ijp2blp(i, j, uv['pol'])
        m = mask[t]
        p1,p2 = a.miriad.pol2str[uv['pol']]

        try:
            i = rewire[i]
        except(KeyError):
            pass
        try:
            j = rewire[j]
        except(KeyError):
            pass

        if i > j:
            i, j = j, i
            d = np.conjugate(d)

        #if i == j and (p1,p2) == ('y','x'):
        #    return None, None, None

        #if same_col(i,j):
        #    return None, None, None

        if t != curtime:
            curtime = t
            aa.set_jultime(t)
            uvo['lst'] = uvo['ra'] = uvo['obsra'] = aa.sidereal_time()
            nints += 1

        preamble = (uvw, t, (i,j))
        return preamble, n.where(m, 0, data), m

    override = {
            'latitud': aa.lat,
            'dec': aa.lat,
            'obsdec': aa.lat,
            'longitu': aa.long,
            'nchan': 1024,
            'nants': 128,
            'ngains': 256,
            'nspect0': 128,
            'telescop': 'PAPER',
            'nints': nints
            }

    uvi = a.miriad.UV(uvfile)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist=' '.join(sys.argv)+'\n')
