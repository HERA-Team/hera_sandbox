#!/opt/local/Library/Frameworks/Python.framework/Versions/2.6/Resources/Python.app/Contents/MacOS/Python
"""
Rotate zenith UV data to a particular source.  Can specify 'zen' to phase data
to zenith, or nothing at all to just remove delay/offset phase components.
"""

import aipy as a, numpy as n, sys, os, optparse

o = optparse.OptionParser()
o.set_usage('dfreq.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
freqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
sync = 500 * (freqs/.150)**-2.5
inttime,sdf = uv['inttime'], uv['sdf']
if True: inttime = 53.5
print 'INTTIME:', inttime
print 'SDF:', sdf
noiselev = sync / n.sqrt(inttime * sdf * 1e9) 
del(uv)

eor_step = n.where(freqs < .140, .080, 0)
signal = (sync + eor_step).astype(n.complex)

# A pipe to use for phasing to a source
def mfunc(uv, p, d, f):
    uvw,t,(i,j) = p
    if i != j: return p, d, f
    noise = n.random.normal(scale=noiselev)
    if True: return p, .03 * (signal + noise), f
    return p, signal + noise, f

# Process data
for filename in args:
    uvofile = filename + 'E'
    print filename,'->',uvofile
    if os.path.exists(uvofile):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(uvofile, status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True)

