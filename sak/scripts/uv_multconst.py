import aipy as a, numpy as n, sys, optparse, os

o = optparse.OptionParser()
o.set_usage('uv_multconst.py file.uv constant [options]')
o.set_description(__doc__)
o.add_option('--div', dest='div', action='store_true',
    help='Divide by the constant instead of multiplying.')
opts,args = o.parse_args(sys.argv[1:])

assert(len(args) == 2)
uv1 = a.miriad.UV(args[0])
const=complex(args[1])

def mfunc(uv, p, d, f):
    if opts.div: d = d/const
    else: d = d*const
    return p, n.where(f, 0, d), f

if opts.div: filename = args[0] + 'D'
else: filename = args[0] + 'M'

print args[0], '->', filename
uvo = a.miriad.UV(filename, status='new')
uvo.init_from_uv(uv1)
uvo.pipe(uv1, mfunc=mfunc, raw=True, append2hist='UVMULTCONST: --sub=%s file=%s\n' % (opts.div, args[1]))
