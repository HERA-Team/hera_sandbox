#! /usr/bin/env python
'''Write out xrfi files'''
import aipy as a, capo, sys, optparse, numpy as np

o = optparse.OptionParser()
o.add_option('--boxside', type='int', default=8, 
              help='Size of the side of a box in flagger.')
o.add_option('--sig', type='float', default=6.,
             help='Number of sigma to flag above.')
o.add_option('--sigl', type='float', default=2.,
             help='Number of sigma for secondary cut in flagger.')
opts,args = o.parse_args(sys.argv[1:])


for f,filename in enumerate(args):
    print '    Omnical npz:', filename
    m,_,_,_ = capo.omni.from_npz(filename)
    chi_flag=capo.xrfi.omni_chisq_to_flags(m['chisq'], K=opts.boxside, sigma=opts.sig, sigl=opts.sigl)
    flags = dict(zip(map(str,m['jds']),chi_flag))
    newname = '.'.join(filename.split('.')[:4]) + '.chisqflag.npz'
    print 'Writing', newname
    np.savez(newname, **flags)

