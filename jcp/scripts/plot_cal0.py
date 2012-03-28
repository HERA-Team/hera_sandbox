#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o)
opts,args = o.parse_args(sys.argv[1:])


f0 = n.load(args[0])
inpnums = f0['inpnums'].flatten()
del(f0)

print args

C,C_gain,times = {},{},[]
for fname in args:
    f = n.load(fname)
    times.append(f['time'])
    for ind, inp in enumerate(inpnums):
        if not C.has_key(inp):
            C[inp] = []
            C_gain[inp] = []
        C[inp].append(f['C'].flatten()[ind])
        C_gain[inp].append(f['C_gain'].flatten()[ind])

times = n.array(times) - 2455000
#print n.isnan(C.values()).sum(), n.isnan(C_gain.values()).sum()
#print C.values()

for inp in inpnums:
    p.subplot(121)
    p.plot(times,C[inp],'.',label=str(inp))
    p.subplot(122)
    p.plot(times,C_gain[inp],'.',label=str(inp))
p.legend()
p.show()
