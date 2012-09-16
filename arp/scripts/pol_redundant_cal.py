#! /usr/bin/env python
"""
Calculate antenna-based corrections to co-align redundant array data.
"""

import aipy as a, numpy as n, capo as C
import sys, optparse, os, cPickle

ij2bl,bl2ij = a.miriad.ij2bl,a.miriad.bl2ij

o = optparse.OptionParser()
o.set_usage('red_cal.py *.uv')
a.scripting.add_standard_options(o, ant=True, pol=True)
o.set_description(__doc__)
o.add_option('--name', type='str', default='polcal',
    help="Output name of solution npz. [default=polcal]")
o.add_option('--plot', action='store_true',
    help="Plot the gain and phase residuals after removing the parameters solved for in this script.")
o.add_option('--verbose', action='store_true',
    help="Print a lot of stuff.")
o.add_option('--calpol', type='str', default='xx',
    help="Polarization to calibrate to.  Should be xx or yy.  Default xx.")
o.add_option('--maxiter', type='int', default=10,
    help="Maximum number of iterations to run in redundant calibration.  Default 10.")
#o.add_option('--refant',type='str',default="0,0",
#    help="Choose antenna in zero referenced grid matrix location <row>,<col>. Don't pick bottom row or rightmost column.")
opts, args = o.parse_args(sys.argv[1:])

uv = a.miriad.UV(args[-1])
fqs = a.cal.get_freqs(uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)

calpol = opts.calpol

for filename in args:
    times, d, f = C.arp.get_dict_of_uv_data([filename], opts.ant, opts.pol, verbose=True)
    
    m_bl = {}
    for bl in d:
        for pol in d[bl]:
            if pol == calpol: continue
            wgt = n.logical_not(f[bl][pol]).astype(n.float)
            wgt_cal = n.logical_not(f[bl][calpol]).astype(n.float)
            g,tau,info = C.arp.redundant_bl_cal(d[bl][calpol], wgt_cal, d[bl][pol], wgt,
                fqs, use_offset=False, maxiter=opts.maxiter)
            #print a.miriad.bl2ij(bl), pol, tau
            m_bl[bl] = tau
            #phs = n.exp(-2j*n.pi*(fqs*tau)); phs.shape = (1,) + phs.shape
            #dat = d[bl][pol] * phs
            #p.subplot(131); C.arp.waterfall(d[bl][calpol], mode='phs')
            #p.subplot(132); C.arp.waterfall(dat, mode='phs')
            #p.subplot(133); C.arp.waterfall(dat * n.conj(d[bl][calpol]), mode='phs')
            #p.plot(n.angle(n.sum(dat * n.conj(d[bl][calpol]), axis=0)))
    #p.ylim(-n.pi,n.pi)
    #p.show()
    taus = C.arp.selfcal_diff(m_bl, {0:0})
    for i,tau in taus.items(): taus[i] = tau[0]
    print taus
            
    ofile = filename + '_%s.pkl' % opts.name
    print 'Writing', ofile
    ofile = open(ofile, 'w')
    cPickle.dump(taus, ofile)
