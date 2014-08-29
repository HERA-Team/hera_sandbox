#! /usr/bin/env python
import aipy as a, numpy as n
import pfits, sys, optparse, glob, ephem, os

o = optparse.OptionParser()
o.add_option('--H_balun',dest='H_balun', type='float', default=-0.024,
    help='Temperature coefficient (dB/K) of balun gain.  Default -0.024')
o.add_option('--H_cable',dest='H_cable', type='float', default=-0.018,
    help='Temperature coefficient (dB/K) of cable gain.  Default -0.018')
o.add_option('--gom',dest='gom', default='',
    help='Comma-delimited list of gain-o-meter antennas, which will have a load temperature also divided out of their gain.')
o.add_option('--T0',dest='T0', type='float', default=300,
    help='Temperature where gain should be 1, in Kelvin.  Default 300')
opts,args = o.parse_args(sys.argv[1:])

goms = map(int, opts.gom.split(','))

def gain(H, T): return 10**(H*(T - opts.T0)/10)

def mfunc(uv, p, d, f):
    crd,t,(i,j) = p
    g = gain(opts.H_balun, uv['t_balun']) * gain(opts.H_cable, uv['t_cable'])
    if i in goms: g *= n.sqrt(uv['t_load'] / opts.T0)
    if j in goms: g *= n.sqrt(uv['t_load'] / opts.T0)
    return p, d/g, f
    
for filename in args:
    print filename, '->', filename+'t'
    if os.path.exists(filename+'t'):
        print 'File exists: skipping'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'t', status='new')
    uvo.init_from_uv(uvi)
    uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist=' '.join(sys.argv)+'\n')
    del(uvo)
