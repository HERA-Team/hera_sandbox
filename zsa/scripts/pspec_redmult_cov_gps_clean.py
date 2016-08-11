#!/usr/bin/env python
import aipy as a, numpy as n
import capo
import optparse, sys, os, random

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True)
o.add_option('-t', '--taps', type='int', default=1,
    help='Taps to use in the PFB.  Default 1, which instead uses windowed FFT')
o.add_option('-b', '--boot', type='int', default=20,
    help='Number of bootstraps.  Default is 20')
o.add_option('--plot', action='store_true',
    help='Generate plots')
o.add_option('-v', store='verbose',  action='store_true',
    help='When used with plotting, have verbose plotting.')
o.add_option('-vv', store='vverbose', action='store_true',
    help='When used with plotting, have very verbose plotting.')
o.add_option('--window', dest='window', default='blackman-harris',
    help='Windowing function to use in delay transform.  Default is blackman-harris.  Options are: ' + ', '.join(a.dsp.WINDOW_FUNC.keys()))
o.add_option('--msg', '-m', action='store', default='No MSF',
    help='Message for this data run. Used to say number of antennas, type of data, etc')
opts,args = o.parse_args(sys.argv[1:])

#seed random number generator.
random.seed() 

PLOT = opts.plot
if PLOT: import pylab as p
if opts.verbose:
    VPLOT=PLOT
elif opts.vverbose:
    VVPLOT=PLOT

NBOOT = opts.boot
NTAPS = opts.tabs
if NTAPS > 1: PFB = True
else PFB = False
WINDOW = opts.window


