#!/usr/bin/env python
"""
A script for filtering using a delay/delay-rate transform.  If a source
is specified, will remove/extract that source.  If none is specified,
will filter/extract in absolute terms.
"""

import aipy as a, numpy as n, os, sys, optparse, ephem

o = optparse.OptionParser()
o.set_usage('uv2pk.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, cal=True)
o.add_option('--clean', dest='clean', type='float', default=1e-4,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
#src = a.phs.RadioFixedBody(ra='12:30', dec='40:00', name='12:30_40:00')

#MIN_CH,MAX_CH,SUBBAND = 120,720,100
MIN_CH,MAX_CH,SUBBAND = 60,804,100
aa.select_chans(n.arange(MIN_CH, MAX_CH))
UV_RES = 1.5
LST_RES = 2*n.pi/36

def uv2bin(u,v,lst,uv_res=UV_RES, lst_res=LST_RES):
    return (int(n.around(u / uv_res) + 4096) * 8192 + int(n.around(v / uv_res) + 4096)) * 8192 + lst/lst_res
def bin2uv(bin, uv_res=UV_RES, lst_res=LST_RES):
    v = ((bin/8192) % 8192 - 4096) * float(uv_res)
    u = (bin / 8192**2 - 4096) * float(uv_res)
    lst = (bin % 8192) * float(LST_RES)
    return u,v, lst

freqs = aa.get_afreqs()
lstbins = n.arange(0,2*n.pi, LST_RES)

for filename in args:
    print 'Reading', filename
    exists = False
    for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
        ch0 = b*SUBBAND/2
        ch1 = ch0 + SUBBAND
        fq = n.average(freqs[ch0:ch1])
        ofile = '%s__%5.3f.npz' % (filename, fq)
        if os.path.exists(ofile):
            print ofile, 'exists.  Skipping...'
            exists = True
            break
    if exists: continue
    _sum,_wgt = {},{}
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
            bm_wgt = aa[0].bm_response(src.get_crds('top'), pol='y')**2
            bm_wgt = bm_wgt.flatten()
        bl = a.miriad.ij2bl(i,j)
        u,v,w = aa.gen_uvw(i,j, src=src)
        u,v = u.flatten(), v.flatten()
        d = d[MIN_CH:MAX_CH]
        val = n.logical_not(f[MIN_CH:MAX_CH]).astype(n.float)
        # mask out overly short baselines adversely affected by xtalk removal
        val *= n.where(n.abs(u) < 12.5, 0, 1)
        d = aa.phs2src(d, src, i, j)
        #print n.sqrt(n.average(n.abs(d)**2)),
        #print n.sqrt(n.average(n.abs(aa.passband(i,j))**2)),
        d /= aa.passband(i,j)
        #print n.sqrt(n.average(n.abs(d)**2))
        #print '------------------'
        d *= val
        if n.average(val) < .5: continue
        if False:  # Looks like this helps reduce sidelobes of strong sources at some expense in SNR
            gain = n.sqrt(n.average(val**2))
            ker = n.fft.ifft(val)
            _d = n.fft.ifft(d)
            _d, info = a.deconv.clean(_d, ker, tol=opts.clean)
            # Not sure if dividing by gain here is the right thing... once clean components are removed, want
            # residuals to be in original data units.
            #if True: _d = info['res'] / gain
            if True: _d = info['res']
            else: _d += info['res'] / gain
            d = n.fft.fft(_d) * val
        # Weight by square of primary beam response in the direction we're projecting for.
        #wgt = bm_wgt**2 * val
        d *= bm_wgt
        w = bm_wgt**2 * val
        for b in range((MAX_CH-MIN_CH) / (SUBBAND/2)):
            ch0 = b*SUBBAND/2
            ch1 = ch0 + SUBBAND
            _u = n.average(u[ch0:ch1])
            _v = n.average(v[ch0:ch1])
            bin = uv2bin(_u,_v,lst,uv_res=UV_RES,lst_res=LST_RES)
            if not _sum.has_key(b): _sum[b],_wgt[b] = {},{}
            _sum[b][bin] = _sum[b].get(bin,0) + d[ch0:ch1]
            _wgt[b][bin] = _wgt[b].get(bin,0) + w[ch0:ch1]

    for b in _sum:
        ch0 = b*SUBBAND/2
        ch1 = ch0 + SUBBAND
        fq = n.average(freqs[ch0:ch1])
        ofile = '%s__%5.3f.npz' % (filename, fq)
        d = {'freqs':freqs[ch0:ch1], 'uvres':UV_RES, 'lstres':LST_RES}
        for bin in _sum[b]:
            d['sum_%d' % (bin)] = _sum[b][bin]
            d['wgt_%d' % (bin)] = _wgt[b][bin]
        print '   Writing', ofile
        n.savez(ofile, **d)

