#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

afreqs = n.linspace(.1,.2,1024)
chans = n.linspace(0,256,1024)

aa = a.cal.get_aa(opts.cal,afreqs)

#p.semilogy(afreqs,aa.passband(0,0))
p.semilogy(chans,aa.passband(0,0))
p.show()
