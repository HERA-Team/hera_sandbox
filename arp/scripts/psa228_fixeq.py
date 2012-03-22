#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys
antis = range(4)
antjs = range(4)
ev_ev = []
ev_od = []
od_od = []
data = {}
for i in antis:
  for j in antjs:
    if i >= j: continue
    print i, j
    bl = a.miriad.ij2bl(i,j)
    uv=a.miriad.UV(sys.argv[-1])
    uv.select("antennae",i,j)
    uv.select("polarization",a.miriad.str2pol["yx"],-1)
    d = n.array([_d for pre,_d,f in uv.all(raw=True)])
    dd = d[1:-1] - .5*(d[2:] + d[:-2])
    stdr = n.std(dd.real,axis=0)
    stdi = n.std(dd.imag,axis=0)
    data[bl] = (stdr,stdi)
    #p.semilogy(stdr*8,"black")
    #p.semilogy(stdi,"blue")
    #p.subplot(211)
    #if (j+i) % 2 == 1: p.semilogy(stdi/stdr,label=bl)
    #else: p.semilogy(stdr/stdi,label=bl)
    #p.subplot(212)
    #p.semilogy(stdr/stdi,label=bl)
    if i % 2 == 0 and j % 2 == 0: ev_ev.append(stdr/stdi)
    elif i % 2 == 1 and j % 2 == 1: od_od.append(stdr/stdi)
    else: ev_od.append(stdi/stdr)
ev_ev = n.median(n.array(ev_ev), axis=0)[100:900]
med_ev_ev = n.median(ev_ev)
val_ev_ev = n.where(n.abs(ev_ev - med_ev_ev) > 5, 0, 1)

ev_od = n.median(n.array(ev_od), axis=0)[100:900]
med_ev_od = n.median(ev_od)
val_ev_od = n.where(n.abs(ev_od - med_ev_od) > 5, 0, 1)

od_od = n.median(n.array(od_od), axis=0)[100:900]
med_od_od = n.median(od_od)
val_od_od = n.where(n.abs(od_od - med_od_od) > 5, 0, 1)

x = n.arange(100,900)
p_ev_ev = n.polyfit(x.compress(val_ev_ev), n.log10(ev_ev.compress(val_ev_ev)), deg=6)
p_ev_od = n.polyfit(x.compress(val_ev_od), n.log10(ev_od.compress(val_ev_od)), deg=6)
p_od_od = n.polyfit(x.compress(val_od_od), n.log10(od_od.compress(val_od_od)), deg=6)

p.subplot(311)
p.semilogy(x, ev_ev, 'k.')
p.semilogy(x, ev_od, 'b.')
p.semilogy(x, od_od, 'g.')
x = n.arange(1024)
p.semilogy(10**n.polyval(p_ev_ev, x), 'k-')
p.semilogy(10**n.polyval(p_ev_od, x), 'b-')
p.semilogy(10**n.polyval(p_od_od, x), 'g-')

for bl in data:
    i,j = a.miriad.bl2ij(bl)
    stdr,stdi = data[bl]
    if i % 2 == 0 and j % 2 == 0:
        r = stdr/stdi
        d = r / 10**n.polyval(p_ev_ev, x)
    elif i % 2 == 1 and j % 2 == 1:
        r = stdr/stdi
        d = r / 10**n.polyval(p_ev_od, x)
    else:
        r = stdi/stdr
        d = r / 10**n.polyval(p_od_od, x)
    p.subplot(312); p.semilogy(x, r, '.')
    p.subplot(313); p.semilogy(x, d, '.')

p.show()
