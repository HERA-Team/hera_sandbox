#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse, re, os

args = sys.argv[1:]

pks = []
for filename in args:
    print 'Reading', filename
    f = n.load(filename)
    kpl,pk,err = f['kpl'], f['pk'], f['err']
    pks.append(pk)
pks = n.array(pks).real
C.arp.waterfall(pks); p.colorbar(shrink=.5); p.show()
pk = n.median(pks, axis=0)
err = n.std(pks, axis=0)
#err = n.sqrt(n.median(n.abs(pks-pk)**2, axis=0))

print 'Writing pspec.npz'
n.savez('pspec.npz', kpl=kpl, pk=pk, err=err)

