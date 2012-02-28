#! /usr/bin/env python

import aipy as a, sys, os, numpy as n, optparse, ephem as e

o = optparse.OptionParser()
o.set_usage('correct_sat952.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-e', '--extract', dest='extract', action='store_true',
    help='Use offset and nchan to extract only some channels.  Default is to echo all channels across.')
o.add_option('-o', '--offset', dest='offset', type='int', default=588,
    help='Offset from channel 0, default is 588')
o.add_option('-n', '--nchan', dest='nchan', type='int', default=37,
    help='Number of channels to slice out, default is 37')

opts, args = o.parse_args(sys.argv[1:])

if not opts.extract:
    uv = a.miriad.UV(args[0])
    p,d,f = uv.read(raw=True)
    opts.nchan = d.size
    opts.offset = 0
chans = n.arange(opts.offset, opts.offset+opts.nchan)
inttime = 2662400/27e6
dt = e.second * inttime

prev_t = None
blcnt = {}

def mfunc(uv, p, d, f):
    global prev_t, blcnt
    uvw, t, (i,j) = p
    bl = a.miriad.ij2bl(i,j)
    if t != prev_t:
        blcnt = {}
        prev_t = t
    else:
        blcnt[bl] = blcnt.get(bl, 0) + 1
        t = t + dt * blcnt[bl]
    return (uvw, t, (i,j)), n.conj(d.take(chans)), f.take(chans)

for filename in args:
    print filename, '->', filename+'c'
    if os.path.exists(filename+'c'):
        print '    File exists... skipping.'
        continue
    uvi = a.miriad.UV(filename)
    uvo = a.miriad.UV(filename+'c', status='new')
    if opts.extract: sfreq = 0.1215 + (opts.offset * uvi['sdf'])
    else: sfreq = uvi['sfreq']
    override = {'sfreq':sfreq,
        'freq':sfreq,
        'freqs':(4,) + (opts.nchan, sfreq, uvi['sdf']) * 4,
        'nchan':opts.nchan,
        'inttime': inttime,
    }
    uvo.init_from_uv(uvi, override=override)
    uvo.pipe(uvi, mfunc=mfunc, raw=True,
        append2hist='CORRECT: select orbcomm chans, fix conj, fix time\n')
    del(uvo)
