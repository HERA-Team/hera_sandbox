#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
    freqs_df = 0.5 * (freqs[1:] + freqs[:-1])
    jy2T = C.pspec.jy2T(freqs)
    jy2T_df = C.pspec.jy2T(freqs_df)
    raw_dat, raw_wgt = 0., 0.
    df_dat, df_wgt = 0., 0.
    dt_dat, dt_wgt = 0., 0.
    prev_d, prev_w = n.zeros(freqs.size), n.zeros(freqs.size)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        w = n.logical_not(f).astype(n.int)
        #d = n.random.normal(size=d.size) * w
        d_dt = (d - prev_d) * n.logical_and(w, prev_w)
        w_dt = 2 * n.logical_and(w, prev_w)
        d_df = (d[1:] - d[:-1]) * n.logical_and(w[1:], w[:-1])
        w_df = 2 * n.logical_and(w[1:], w[:-1])
        raw_dat += n.abs(d)**2
        raw_wgt += w
        df_dat += n.abs(d_df)**2
        df_wgt += w_df
        dt_dat += n.abs(d_dt)**2
        dt_wgt += w_dt
        prev_d, prev_w = d, w
    raw = n.sqrt(n.where(raw_wgt > 0, raw_dat / raw_wgt, 0)) * jy2T
    df = n.sqrt(n.where(df_wgt > 0, df_dat / df_wgt, 0)) * jy2T_df
    dt = n.sqrt(n.where(dt_wgt > 0, dt_dat / dt_wgt, 0)) * jy2T
    p.plot(freqs, raw, 'k')
    p.plot(freqs_df, df, 'r')
    p.plot(freqs, dt, 'g')
p.xlabel('Frequency (GHz)')
p.ylabel('Temperature (mK)')
p.show()
