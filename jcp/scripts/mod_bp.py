#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, optparse, sys

o = optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
opts,args = o.parse_args(sys.argv[1:])

afreqs = n.linspace(.1,.2,2048)
chans = n.linspace(0,2048,2048)

aa = a.cal.get_aa(opts.cal,afreqs)

bp = aa.passband(0,0)
lo,hi = 500,1700
bp[:lo] = bp[lo]
bp[hi:] = bp[hi]

bp_poly = n.polyfit(afreqs,n.sqrt(bp),deg=9)
print [bp_poly / aa[0].amp]
fit = n.polyval(bp_poly,afreqs)**2

#p.semilogy(afreqs,aa.passband(0,0))
p.semilogy(chans,aa.passband(0,0))
p.semilogy(chans,bp)
p.semilogy(chans,fit)
p.show()
