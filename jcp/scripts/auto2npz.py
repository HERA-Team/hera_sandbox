#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('-c', '--chan', dest='chan', type='int', default=128,
    help='Channel to use.')
o.add_option('-d', '--decimate', dest='decimate', default=1,
    help='Decimate the actual data.')
o.add_option('-o', '--outfile', dest='outfile', default='auto2npz.npz',
    help='The name of the npz file that is created.')
opts,args = o.parse_args(sys.argv[1:])

aa = a.cal.get_aa(opts.cal, n.array([.150]))

auto,lsts = [],[]
my_auto = None
for filename in args:
    print 'Reading', filename
    uv = a.miriad.UV(filename)
    uv.select('antennae', 8, 8)
    uv.select('polarization',a.miriad.str2pol['xx'],a.miriad.str2pol['xx'])
    uv.select('decimate', opts.decimate, 0)
    for (crd,t,(i,j)),d,f in uv.all(raw=True):
        if my_auto is None: my_auto = i
        if i != my_auto: continue
        auto.append(d[opts.chan]/aa.passband(i,j))
        aa.set_jultime(t)
        lsts.append(aa.sidereal_time())
    
auto = n.array(auto).real
lsts = n.array(lsts)

p.plot(lsts,auto)
p.show()

n.savez(opts.outfile,lsts=lsts,auto=auto)
