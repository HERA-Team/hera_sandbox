import aipy, numpy as np
import sys, os, optparse

o = optparse.OptionParser()
o.add_option('--npzdir', dest='npzdir', help='name of the directory holding outputs of get_chisq4tau.py. Include "/" pls.')
opts,args = o.parse_args(sys.argv[1:])

freqs = np.linspace(.1,.2,num=203)

for uvfile in args:
	uvofile = uvfile+'T'
	print uvfile,'->',uvofile
	if os.path.exists(uvofile):
	    print uvofile, 'exists, skipping.'
	    continue
	#initialize
	uvi = aipy.miriad.UV(uvfile)
	uvo = aipy.miriad.UV(uvofile, status='new')

	#locate npz file
	name = uvfile.split('.')
	name[0]='tau'
	name='.'.join(name)+'.npz'

	npzfile = opts.npzdir+name
	print '    Using: '+npzfile
	
	try: 
		C = np.load(npzfile)
		print '    npz loaded'
	except IOError:
		print 'npz file %s does not exist'%npzfile
		continue
	
	#C = np.load(npzfile)
	X2 = C['chisquared']
	taus= C['tau']

	#construct freq x time array of taus that minimize V across the array
	mins = np.argmin(X2,axis=0)
	final_taus = np.zeros(mins.shape)

	for t in range(mins.shape[0]): final_taus[t,:] = taus[mins[t,:]]

	c_angle = np.exp(-2*np.pi*1.j*freqs*final_taus)
	
	angles = {'cA':c_angle}
	
	#we need to apply these to all yx visibilities each integration
	
	curtime = 0
	t_index = -1

	def mfunc(uv,p,d,f):
		
		if uv['pol']!=-8: return p,d,f
		else:
			global curtime
			global t_index
			global angles
			
			uvw, t, (i,j) = p
			
			if t!=curtime:
				t_index+=1
				curtime = t
			
			angle = angles['cA'][t_index,:]
			d=d*angle
			return p,d,f

	uvo.init_from_uv(uvi)
	uvo.pipe(uvi, mfunc=mfunc, raw=True, append2hist="APPLY_TAUCALC: yx phased to minimize V according to %s \n"%npzfile)
	del(uvo)
	del(curtime)
	del(t_index)
	del(angles)
    
    	
    	
    	
