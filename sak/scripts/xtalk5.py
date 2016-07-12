#!/usr/bin/env python
"""
Remove cross-talk by subtracting long-time (nightly) average from the data.
This new iteration of xtalk.py is designed to skip the sun in it's average
"""
import aipy as a
import numpy as np
import optparse,sys,os,ephem
from glob import glob

o = optparse.OptionParser()
o.set_usage('xtalk5.py [options] *.uv')
o.set_description(__doc__)
o.add_option('-C', '--cal', action = 'store', default = None, help = 'calfile for processing uv files')
o.add_option('-V', '--verbose', action = 'store_true', default=False, help='toggle verbosity')
o.add_option('-s', '--saveloc',action='store', default=None, help='save location for npz')
opts,args = o.parse_args(sys.argv[1:])

xtalk,cnt = {},{}

print "Reading calfile %s"%opts.cal,
sys.stdout.flush()
aa = a.cal.get_aa(opts.cal, np.array([.15]))
print "Done"

#pyEphem stuff
sun = ephem.Sun()

sa = ephem.Observer()
sa.lon = aa.lon
sa.lat = aa.lat

julDelta = 2415020.

for uvfile in args:
	if uvfile[0:3]=='zen' or uvfile[0:3]=='lst': JD = float(uvfile[4:17])
	elif uvfile[0:3]=='psa' or uvfile[0:3]=='245': JD = float(uvfile[12:25])
	else: print 'JD issue!';sys.exit()
	sa.date = JD-julDelta
	sun.compute(sa)
	sunpos = sun.alt, sun.az
	print '='*(len(uvfile)+8)
	print 'Reading %s'%uvfile
	#DEBUG
	if opts.verbose:
		print sa.date
		print 'Sun alt,az=',sunpos[0],sunpos[1]
	if sunpos[0] > -0.1:
		if opts.verbose: print 'Sun is up, skipping file'
		continue
	else:
		if opts.verbose: print 'Sun is down, proceeding with xtalk modelling'
	#DEBUG
	
	
	uv = a.miriad.UV(uvfile)
	
	for (uvw,t,(i,j)),d,f in uv.all(raw=True):
	    pol = a.miriad.pol2str[uv['pol']]
	    bl = (i,j)
	    if i == j: continue
	    if not pol in xtalk.keys(): 
	        xtalk[pol] = {}
	        cnt[pol] = {}
	    if not bl in xtalk[pol].keys():
	        xtalk[pol][bl] = np.zeros(uv['nchan'],dtype=complex)
	        cnt[pol][bl] = np.zeros(uv['nchan'])
	    xtalk[pol][bl] += np.where(f,0,d)
	    cnt[pol][bl] += np.where(f,0,1)
	del uv
print '='*(len(args[-1])+8)

for pol in xtalk:
    for bl in xtalk[pol]:
        xtalk[pol][bl] /= cnt[pol][bl]

if opts.saveloc is not None: np.savez('%s/xtalk.npz'%opts.saveloc,xtalk=xtalk)

for infile in args:
    outfile = infile+'X'
    print infile,'-->',outfile
    if os.path.exists(outfile):
        print 'File exists, skipping...'
        continue
    uvi = a.miriad.UV(infile)
    uvo = a.miriad.UV(outfile,status='new')
    uvo.init_from_uv(uvi)

    def mfunc(uv,p,d):
        uvw,t,(i,j) = p
        pol = a.miriad.pol2str[uv['pol']]
        bl = (i,j)
        if i != j: d -= xtalk[pol][bl]
        return p,d

    uvo.pipe(uvi,mfunc=mfunc,append2hist='XTALK5 \n')
