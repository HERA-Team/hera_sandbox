#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import pfits, sys, optparse, ephem

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('-j','--jdbin',dest='jdbin', type='float', default=600,
    help='Time interval for averaging temperatures (seconds).  Default 600')
opts,args = o.parse_args(sys.argv[1:])

# Collect temperature data
jds, temps = [], []
for f in args:
    print f
    hdu = pfits.FITS(f).get_hdus()[1]
    dat = hdu.get_data()
    jds.append(dat['DMJD'] + 2400000.5)
    temps.append(dat['TEMPERATURE'])

jds = n.concatenate(jds)
temps = n.concatenate(temps)
JDBIN = ephem.second * opts.jdbin
nbins = int((jds[-1]-jds[0])/JDBIN)
jwgts,bins = n.histogram(jds, bins=nbins)
jtemps,bins = n.histogram(jds, weights=temps, bins=nbins)
jd_bins = .5*(bins[1:] + bins[:-1])
temp_prof = jtemps/jwgts

if not opts.cal is None:
    aa = a.cal.get_aa(opts.cal, .1, .1, 1)
    lsts = []
    for jd in jd_bins:
        aa.date = a.phs.juldate2ephem(jd)
        lsts.append(aa.sidereal_time())
    p.subplot(211)
    p.plot(lsts, temp_prof, '.')
    p.subplot(212)

p.plot(jd_bins, temp_prof)
p.show()
