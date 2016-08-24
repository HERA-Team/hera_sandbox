#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True)
opts,args = o.parse_args(sys.argv[1:])

jy2T = {}
jy2T_df = {}
raw_dat, raw_wgt = {}, {}
df_dat, df_wgt = {}, {}
dt_dat, dt_wgt = {}, {}
#raw_dat, raw_wgt = 0., 0.
#df_dat, df_wgt = 0., 0.
#dt_dat, dt_wgt = 0., 0.
for filename in args:
    print 'Reading', filename
    ftype = filename.split('.')[-1]
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, opts.ant, opts.pol)
    print 'B =', uv['sdf']
    if not jy2T.has_key(ftype):
        freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
        freqs_df = 0.5 * (freqs[1:] + freqs[:-1])
        jy2T[ftype] = C.pspec.jy2T(freqs)
        jy2T_df[ftype] = C.pspec.jy2T(freqs_df)
    #prev_d, prev_w = n.zeros(freqs.size), n.zeros(freqs.size)
    zeros = n.zeros(jy2T[ftype].size)
    prev_d, prev_w = {}, {}
    curtime = None
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        if t != curtime:
            #if curtime != None: print 'dt =', (t - curtime) / a.ephem.second
            curtime = t
        # Re-calibrate to new picA flux scale
        d *= (363.6/424) * (freqs/.160)**(-.76+.95)
        bl = a.miriad.ij2bl(i,j)
        w = n.logical_not(f).astype(n.int)
        #d = n.random.normal(size=d.size) * w
        d_dt = (d - prev_d.get(bl, zeros)) * n.logical_and(w, prev_w.get(bl,zeros))
        w_dt = 2 * n.logical_and(w, prev_w.get(bl,zeros))
        d_df = (d[1:] - d[:-1]) * n.logical_and(w[1:], w[:-1])
        w_df = 2 * n.logical_and(w[1:], w[:-1])
        raw_dat[ftype] = raw_dat.get(ftype,0) + n.abs(d)**2
        raw_wgt[ftype] = raw_wgt.get(ftype,0) + w
        df_dat[ftype] = df_dat.get(ftype,0) + n.abs(d_df)**2
        df_wgt[ftype] = df_wgt.get(ftype,0) + w_df
        dt_dat[ftype] = dt_dat.get(ftype,0) + n.abs(d_dt)**2
        dt_wgt[ftype] = dt_wgt.get(ftype,0) + w_dt
        prev_d[bl] = d
        prev_w[bl] = w
for ftype in jy2T.keys():
    raw = n.sqrt(n.where(raw_wgt[ftype] > 0, raw_dat[ftype] / raw_wgt[ftype], 0)) * jy2T[ftype]
    df = n.sqrt(n.where(df_wgt[ftype] > 0, df_dat[ftype] / df_wgt[ftype], 0)) * jy2T_df[ftype]
    dt = n.sqrt(n.where(dt_wgt[ftype] > 0, dt_dat[ftype] / dt_wgt[ftype], 0)) * jy2T[ftype]
    p.semilogy(freqs, raw, 'k')
    p.plot(freqs_df, df, 'm')
    #p.plot(freqs, dt, 'c')
p.grid(which='both')
#p.xlim(.115,.187)
#p.ylim(1,3e3)
p.ylim([1,100])
p.xlim([0.110,0.190])
p.xlabel('Frequency [GHz]')
p.ylabel('Brightness Temperature [mK]')
p.show()
