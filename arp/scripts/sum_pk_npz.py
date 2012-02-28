#!/usr/bin/env python
"""
A script for filtering using a delay/delay-rate transform.  If a source
is specified, will remove/extract that source.  If none is specified,
will filter/extract in absolute terms.
"""

import numpy as n
import os, sys, optparse

o = optparse.OptionParser()
o.set_usage('sum_pk_npz.py [options] *.uv')
o.set_description(__doc__)

data = {}
for filename in args:
    print 'Loading', filename
    f = n.load(filename)
    fq = n.average(f['freqs'])
    if not data.has_key(fq):
        data[fq] = {}
    for k in f.files:
        if k == 'freqs':
            data[fq]['freqs'] = f['freqs']
        elif k == 'uvres':
            data[fq]['uvres'] = f['uvres']
        elif k == 'lstres':
            data[fq]['lstres'] = f['lstres']
        elif k.startswith('sum'):
            data[fq][k] = data[fq].get(k,0) + f[k]
        elif k.startswith('wgt'):
            data[fq][k] = data[fq].get(k,0) + f[k]
        else: raise ValueError('Unrecognized keyword: %s' % k)

for fq in data:
    ofile = 'sumpk__%5.3f.npz' % (fq)
    print '   Writing', ofile
    n.savez(ofile, **data[fq])
    
