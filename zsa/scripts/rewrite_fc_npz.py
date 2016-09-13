#! /usr/bin/env python
'''Rewrite firstcal npz files to be consistent with omnical.'''
import numpy as np
import sys,optparse

o = optparse.OptionParser()
o.add_option('--pol', '-p',  action='store', help='Polarization string')
opts,args = o.parse_args(sys.argv[1:])

pols = ['xx', 'yy', 'xy', 'yx']

for f in args:
    d = np.load(f)
    pol_test = f.split('.')[3]
    if pol_test in pols: pol = pol_test #test to see if polarization is 
    elif opts.pol in pols: pol = opts.pol
    else:
        raise RuntimeError('No polarization string provided.')
    data = {}
    for k in d.keys(): 
        if k.isdigit(): 
            data[k+'x'] = d[k].reshape(1,d[k].shape[0]) # reshape to be consistent with omni_apply.py
        elif k.endswith('d'): data['d'+k[:-1]] = d[k]
        else: data[k] = d[k]
    np.savez(f, **data)
