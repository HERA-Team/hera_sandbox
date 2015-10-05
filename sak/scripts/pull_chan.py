#!~/src/anaconda/envs/PAPER/bin/python

"""Select channel(s) from a uv file and write to a new uv file"""

import aipy as a, sys, optparse, os, numpy as np

o = optparse.OptionParser()
o.set_usage('pull_chan.py [options] *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True)

opts, args = o.parse_args(sys.argv[1:])

# Parse command-line options
uv = a.miriad.UV(args[0])
(j,t,j),j = uv.read()
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
#aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
#aa.select_chans(chans)
#aa.set_active_pol(opts.pol)
del(uv)
#print chans

def mfunc(uv,p,d,f):
	d,f = d.take(chans), f.take(chans)
	return p, np.where(f,0,d), f

for filename in args:
	outname = filename+'f'
	if os.path.exists(outname):
		print '    File exists... skipping.'
		continue
	print filename,'->',outname
	uvi = a.miriad.UV(filename)
	uvo = a.miriad.UV(outname,status='new')
	uvo.init_from_uv(uvi)
	uvo.pipe(uvi,mfunc=mfunc,raw=True,append2hist='CHAN SELECT: '+opts.chan+'\n')
	
	
