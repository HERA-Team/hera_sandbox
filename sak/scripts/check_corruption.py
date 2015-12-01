import aipy, optparse, sys, os

o = optparse.OptionParser()
o.set_usage('check_corruption.py [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

for uvfile in args:
	a = 0
	uv = aipy.miriad.UV(uvfile)
	for preamble, data in uv.all():
		uvw, t, (i,j) = preamble
		a = t
