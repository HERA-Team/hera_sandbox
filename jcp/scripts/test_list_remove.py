#! /usr/bin/env python

import aipy as a
import optparse, sys, os

o = optparse.OptionParser()
o.add_option('--stats', default='all',
    help="Statistics to include in npz meta data file.  Options are 'all' (default), 'cnt' (the number of integrations in each lst bin), 'minmax' (the min and max amplitude of the visibilities in each lst bin', 'median' (the median value in each lst bin), and 'var' (the variance in each lst bin).  Multiple values can be chosen with commas e.g. '--stats=minmax,var'.")
opts, args = o.parse_args(sys.argv[1:])

opts.stats = map(str, opts.stats.split(','))
print opts.stats.__class__.__name__

uvi = a.miriad.UV(args[0])
uvo = a.miriad.UV('test_list_remove.uv', status='new')
uvo.init_from_uv(uvi)

print opts.stats
if opts.stats == ['all']: opts.stats = ['cnt','min','max','median','var']
delstats = []
for stat in opts.stats:
    if stat in uvo.vars(): delstats.append(stat)
print opts.stats
print delstats
for stat in delstats:
    opts.stats.remove(stat)
print opts.stats

del(uvo)
os.system('rm -rf test_list_remove.uv')
