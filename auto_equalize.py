#! /usr/bin/env python
import aipy as a, numpy as n
import sys, optparse

o = optparse.OptionParser()
o.add_option('--deg', dest='deg', type='int', default=8,
    help='Degree of polynomial to fit to passband of each antenna. Default 8')
opts,args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

dat,flg = {}, {}
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    uv.select('auto',-1,-1)
    for (crd,t,(i,j)),d, f in uv.all(raw=True):
        if not dat.has_key(i):
            dat[i] = []
            flg[i] = []
        dat[i].append(d.real)
        flg[i].append(f)

for i in dat:
    flg[i] = n.array(flg[i])
    dat[i] = n.where(flg[i], 0, n.array(dat[i]))

model_auto = (freqs / .150)**-2.5
model_auto.shape = (1, model_auto.size)
#ddc_poly = [65123029498.345612, -97672196551.106552, 64042861492.03476, -24111913057.711517, 5761340839.0835114, -912373536.43163371, 97137219.224425375, -6896379.6290631033, 314561.96293171198, -8363.013373925829, 98.57747247358202]
#ddc_bp = n.polyval(ddc_poly, freqs)

bp, poly, gvt = {}, {}, {}

import pylab as p
for i in dat:
    if i != 1: dat[i] /= model_auto
    #dat[i] /= ddc_bp
    d = n.ma.array(dat[i], mask=flg[i])
    avg = n.ma.average(d, axis=1); avg.shape = (avg.size, 1)
    sum = n.sum(dat[i] / avg, axis=0)
    wgt = n.sum(n.logical_not(flg[i]).astype(n.float), axis=0)
    bp[i] = n.sqrt(sum / wgt.clip(1,n.Inf))
    valid = n.where(wgt > .5 * wgt.max(), 1, 0)
    vfqs = n.compress(valid, freqs)
    vbp = n.compress(valid, bp[i])
    poly[i] = n.polyfit(vfqs, vbp, deg=opts.deg)
    sbp = n.polyval(poly[i], freqs).clip(.5,2); sbp.shape = (1, sbp.size)
    d /= sbp**2
    gvt[i] = n.ma.average(d, axis=1)
    #p.plot(freqs, sbp.squeeze(), label=str(i))
#p.legend(); p.show()

amp = {}
for i in bp:
    #amp[i] = n.sqrt(n.average((gvt[i] / gvt[0])).clip(.25,4))
    amp[i] = n.sqrt(n.average(gvt[i]))
    #amp[i] = n.sqrt(n.median(gvt[i] / gvt[0]))
    print i, amp[i]/amp[0], list(poly[i]),','
    p.subplot(211)
    #p.plot(bp[i])
    valid = n.where(bp[i] > 0, 1, 0)
    vbp = n.compress(valid, bp[i])
    vfq = n.compress(valid, freqs)
    p.plot(vfq, vbp / n.polyval(poly[i], vfq))
    p.subplot(212)
    p.plot(gvt[i] / amp[i]**2, label=str(i))
p.legend()
p.show()

