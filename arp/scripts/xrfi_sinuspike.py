#! /usr/bin/env python
import aipy as a, numpy as n, capo as C
import optparse, sys, os

def kurtosis(d,w,axis=0,rm_avg=True):
    '''Only works for wgts consisting of 1,0'''
    wgt = n.sum(w, axis=axis).clip(1,n.Inf)
    if rm_avg:
        # XXX this step may fail if axis!=0
        avg = n.sum(d, axis=axis) / wgt
        res = d - avg
    else: res = d
    m4 = n.sum(n.abs(res)**4, axis=axis) / wgt
    m2 = n.sum(n.abs(res)**2, axis=axis) / wgt
    return m4/n.where(m2==0, 1, m2**2) - 3

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True)
o.set_usage('xrfi_sinuspike.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-c', '--chan', dest='chan',
    help='Manually flag channels before xrfi processing.  Options are "<chan1 #>,..." (a list of channels), or "<chan1 #>_<chan2 #>" (a range of channels).  Default is None.')
o.add_option('-n','--nsig', dest='nsig', type='float', 
    help='Number of standard deviations above mean to flag')
o.add_option('--df', dest='df', type='float', 
    help='Kurtosis number above which to flag channels')
o.add_option('--dt', dest='dt', type='float',
    help='Kurtosis number above which to flag integrations')
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


for uvfile in args:
    uvofile = uvfile+'R'
    print uvfile,'->',uvofile
    if os.path.exists(uvofile):
        print uvofile, 'exists, skipping.'
        continue
    uvi = a.miriad.UV(uvfile)
    fqs = a.cal.get_freqs(uvi['sdf'], uvi['sfreq'], uvi['nchan'])
    a.scripting.uv_selector(uvi, opts.ant)
    # Gather all data and each time step
    data,mask,times = {}, {}, []
    for (uvw,t,(i,j)), d, f in uvi.all(raw=True):
        if len(times) == 0 or times[-1] != t: times.append(t)
        bl = a.miriad.ij2bl(i,j)
        pol = uvi['pol']
        if not pol in data:
            data[pol] = {}
            mask[pol] = {}
        if not bl in data[pol]:
            data[pol][bl] = {}
            mask[pol][bl] = {}
        # Manually flag channels
        f[chans] = 1
        mdl,f = C.arp.sinuspike(d,fqs,f,nsig=opts.nsig)
        mask[pol][bl][t] = f
        data[pol][bl][t] = n.where(f, 0, d - mdl)

    # Generate statistical mask
    for pol in data:
      for bl in data[pol]:
        i, j = a.miriad.bl2ij(bl)
        data_times = data[pol][bl].keys()
        data_times.sort()
        d = n.array([data[pol][bl][t] for t in data_times])
        m = n.array([mask[pol][bl][t] for t in data_times])
        w = n.logical_not(m).astype(n.int)
        if opts.df != None:
            k_fq = kurtosis(d,w,axis=0)
            #import pylab; pylab.plot(k_fq); pylab.show()
            m[:,n.where(k_fq>opts.df)] = 1
        if opts.dt != None:
            k_t = kurtosis(d,w,axis=1,rm_avg=False)
            #import pylab; pylab.plot(k_t); pylab.show()
            m[n.where(k_t>opts.dt)] = 1
        for i, t in enumerate(data_times): mask[pol][bl][t] |= m[i]

    new_mask = {}
    for pol in mask:
      for bl in mask[pol]:
        for t in mask[pol][bl]:
            new_mask[t] = new_mask.get(t,0)+mask[pol][bl][t].astype(n.int)
    for t in new_mask:
        m = n.where(new_mask[t] >= opts.thresh, 1, 0)
        for pol in mask:
          for bl in mask[pol]:
            mask[pol][bl][t] = m

    # Generate a pipe for applying the mask to data as it comes in.
    def rfi_mfunc(uv, preamble, data, flags):
        uvw, t, (i,j) = preamble
        bl = a.miriad.ij2bl(i,j)
        m = mask.values()[0].values()[0][t]
        return preamble, n.where(m, 0, data), m

    del(uvi)
    uvi = a.miriad.UV(uvfile)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=rfi_mfunc, raw=True, append2hist=' '.join(sys.argv)+'\n')



