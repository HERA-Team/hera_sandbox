import aipy as a, numpy as n, sys, optparse, os

o = optparse.OptionParser()
o.set_usage('uv_multconst.py *.uv [options]')
o.set_description(__doc__)
o.add_option('--const',dest='const',help='constant to multiply/divide by',default=1.)
o.add_option('--div', dest='div', action='store_true',
    help='Divide by the constant instead of multiplying.')
opts,args = o.parse_args(sys.argv[1:])

const = float(opts.const)

def mfunc(uv, p, d, f):
    if opts.div: d = d/const
    else: d = d*const
    return p, n.where(f, 0, d), f

for file in args:
	uv = a.miriad.UV(file)
	if opts.div: filename = file + 'D'
	else: filename = file + 'M'
	print file, '->', filename
	uvo = a.miriad.UV(filename, status='new')
	uvo.init_from_uv(uv)
	uvo.pipe(uv, mfunc=mfunc, raw=True, append2hist='UVMULTCONST: --div=%s --const=%s file=%s\n' % (opts.div, opts.const, args.join(' ')))
