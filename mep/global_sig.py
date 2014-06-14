#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import useful_functions as uf

def get_coeffs(calfile,na = 32,freqs = n.arange(.1,.2,.01)):
	"""
	This function calculates the coefficients in front of the global signal 
	for the antenna array described in the calfile. The coefficients are 
	determined by the integral of A(l,m)exp[-i2pi(ul+vm)]dldm.
	"""
	aa = a.cal.get_aa(calfile, n.array([.150])) #get antenna array
	im = a.img.Img(size=200, res=.5) #make an image of the sky to get sky coords
	tx,ty,tz = im.get_top(center=(200,200)) #get coords of the zenith?
	valid = n.logical_not(tx.mask)
	tx,ty,tz = tx.flatten(),ty.flatten(),tz.flatten()
	#beam response for an antenna pointing at (tx,ty,tz) with a polarization in x direction
	#amp = A(theta) in notes
	amp = aa[0].bm_response((tx,ty,tz),pol='x')**2 
	#coefficient array: rows baselines; cols frequencies
	coeffs = n.zeros([(na*na-na)/2,len(freqs)],dtype=n.complex)
	# loop through all baselines
	baselines = n.zeros([(na*na-na)/2,3])
	ll=0
	for ii in n.arange(na):
	    for jj in n.arange(ii):
	        bx,by,bz = aa.get_baseline(ii,jj,'z') 
	        #the baseline for antennas ii and jj 
	        print ii,jj,[bx,by,bz]
	        baselines[ll] = [bx,by,bz]
	        kk=0
	        for fq in freqs: #loop over frequencies
	            phs = n.exp(-2j*n.pi*fq * (bx*tx+by*ty+bz*tz)) #fringe pattern
	            phs.shape = amp.shape = im.uv.shape 
	            amp = n.where(valid, amp, 0)
	            phs = n.where(valid, phs, 0)

	            dc_response = n.sum(amp*phs)/n.sum(amp) #this integrates the amplitude * fringe pattern; units mK?
	            jy_response = n.real(dc_response * 100 / C.pspec.jy2T(fq)) # jy2T converts flux density in jansky to temp in mK
	            print '\t',fq, dc_response, jy_response
	            coeffs[ll,kk] = dc_response
	            kk+=1
	        ll+=1
	n.savez('./coeff_data/{0}_data'.format(calfile),baselines=baselines,frequencies=freqs,coeffs=coeffs)
	# n.savetxt('{0}_baselines.txt'.format(calfile),baselines)
	# n.savetxt('{0}_frequencies.txt'.format(calfile),freqs)
	# n.savetxt('{0}_coeffs.txt'.format(calfile),coeffs)
	return baselines,freqs,coeffs

def read_coeffs(calfile,na=32):
	"""
	This function reads the coefficients in front of the global signal 
	from the data file created in get_coeffs run with the smae calfile. 
	"""
	# baselines = n.loadtxt('{0}_baselines.txt'.format('psa898_v002'))
	# freqs = n.loadtxt('{0}_frequencies.txt'.format('psa898_v002'))
	# coeffs = n.loadtxt('{0}_coeffs.txt'.format('psa898_v002'))	
	f = n.load('./coeff_data/{0}_data.npz'.format(calfile))
	freqs = f['frequencies']
	baselines = f['baselines']
	coeffs = f['coeffs']
	baselines = baselines[:na]
	if len(freqs)>1: freqs = freqs[:na]
	coeffs = coeffs[:na]
	return baselines,freqs,coeffs

def test_regress(baselines,coeffs,gs=1,n_sig=5,na=32,readFromFile=False):
	"""
	This function tests the linear regression method of recovering the global 
	signal by generating a visibility VV with random normal noise of amplitude 
	n_sig and then performing a linear regression on VV and the coefficients 
	from get_coeffs(). 
	"""
	print coeffs.shape
	VV = gs*coeffs + n.random.normal(loc=0.0,scale=n_sig,size=[len(coeffs),1])*n.exp(2*n.pi*1j*n.random.rand())
	print VV.shape
	gs_recov,n_recov,RR,redchi = uf.linear_fit(coeffs,VV)
	#gs_recov,n_recov,redchi = uf.linear_fit_new(coeffs,VV)
	print "true gs = {0}\trecovered gs = {1}".format(gs,gs_recov)
	print "true n = {0}\trecovered n = {1}".format(0.0,n_recov)
	print "R = ",RR
	print 'chi = ',redchi
	return gs_recov, redchi

def test_regress_thru_origin(baselines,coeffs,gs=1,n_sig=5,na=32,readFromFile=False):
	"""
	This function tests the linear regression method of recovering the global 
	signal by generating a visibility VV with random normal noise of amplitude 
	n_sig and then performing a linear regression on VV and the coefficients 
	from get_coeffs(). 
	"""
	print coeffs.shape
	VV = gs*coeffs + n.random.normal(loc=0.0,scale=n_sig,size=[len(coeffs),1])
	VV = VV*n.exp(2*n.pi*1j*n.random.rand())
	print VV.shape

	gs_recov,redchi = uf.line_thru_origin_fit(coeffs,VV)

	print "true gs = {0}\trecovered gs = {1}".format(gs,gs_recov)
	print 'chi = ',redchi
	return gs_recov, redchi

def test_regress_vary_n(baselines,coeffs,nants=32,restrictChi=False):
	"""
	This function runs many tests of the linear regression, varying the 
	amount of noise introduced into the data. I hard-coded the global 
	signal strength to be 1 so that it is easy to compare the magnitude
	of the recovered and true global signals and the amplitude of the noise.
	"""
	for jj in n.arange(100):
		gs_diff = n.zeros(20)
		n_sigs = n.logspace(-3,1,num=20)
		for ii,n_sig in enumerate(n_sigs):
			print ii
			#if ii==0 and jj==0:rFF=False 
			#else:rFF=True
			rFF=True
			gs_recov, redchi = test_regress(baselines,coeffs,gs=1,n_sig=n_sig,na=nants,readFromFile=rFF)
			if restrictChi:
				if n.absolute(redchi-1)<0.1:
					gs_diff[ii] = gs_recov - 1.
				else:
					gs_diff[ii] = None
			else:
				gs_diff[ii] = gs_recov - 1.
		p.scatter(n_sigs,n.absolute(gs_diff))
		p.xscale('log')
		p.yscale('log')
		p.xlim(1e-4,1e2)
	p.xlabel('Amplitude of noise relative to global signal\n(I.e. true global signal amplitude is 1)')
	#p.ylabel('Recovered global signal (true gs = 1)')
	p.ylabel('Difference between true and recovered global signal')
	#p.show()
	if restrictChi:
		p.savefig('./figures/gs_diff_vs_n_good_chi.pdf')
	else:
		p.savefig('./figures/gs_diff_vs_n.pdf')
	p.clf()

def test_regress_vary_na(baselines,coeffs,nants=32,restrictChi=False):
	"""
	This function runs many tests of the linear regression, varying the number 
	of antennae in the array. Again, the global signal is hard-coded to be 1.
	"""
	for jj in n.arange(100):
		nas = n.arange(2,nants)
		gs_diff = n.zeros(len(nas))
		for ii,na in enumerate(nas):
			gs_recov, redchi = test_regress(baselines,coeffs,gs=1,n_sig=.1,na=na,readFromFile=True)
			if restrictChi:
				if n.absolute(redchi-1)<0.1:
					gs_diff[ii] = gs_recov - 1.
				else:
					gs_diff[ii] = None
			else:
				gs_diff[ii] = gs_recov - 1.
		p.scatter(nas,gs_diff)
	p.xlabel('Number of antenna')
	p.ylabel('Difference between true and recovered global signal')
	#p.show()
	if restrictChi:
		p.savefig('./figures/gs_diff_vs_na_good_chi.pdf')
	else:
		p.savefig('./figures/gs_diff_vs_na.pdf')
	p.clf()
	

def test_regress_vary_bsln(baselines,coeffs,nants=32,restrictChi=False):
	"""
	This is the exact same function as test_regress_vary_na except that it 
	plots the number of baselines on the x axis instead of the number of 
	antennae.
	"""
	for jj in n.arange(100):
		nas = n.arange(2,nants)
		gs_diff = n.zeros(len(nas))
		for ii,na in enumerate(nas):
			gs_recov, redchi = test_regress(baselines,coeffs,gs=1,n_sig=.1,na=na,readFromFile=True)
			if restrictChi:
				if n.absolute(redchi-1)<0.1:
					gs_diff[ii] = gs_recov - 1.
				else:
					gs_diff[ii] = None
			else:
				gs_diff[ii] = gs_recov - 1.
		p.scatter(nas*(nas-1)/2,gs_diff)
	p.xlabel('Number of baselines')
	p.ylabel('Difference between true and recovered global signal')
	#p.show()
	if restrictChi:
		p.savefig('./figures/gs_diff_vs_bsln_good_chi.pdf')
	else:
		p.savefig('./figures/gs_diff_vs_bsln.pdf')
	p.clf()


if __name__=='__main__': 
	#calfile='psa898_v002'
	calfile='basic_amp_aa'
	na=8
	baselines,freqs,coeffs = get_coeffs(calfile,na=na,freqs=n.array([.1,]))
	#baselines,freqs,coeffs = read_coeffs(calfile,na=na)
	print coeffs
	test_regress_vary_n(baselines,coeffs,nants=na) 
	test_regress_vary_na(baselines,coeffs,nants=na)
	test_regress_vary_bsln(baselines,coeffs,nants=na)


	# if readFromFile:
	# 	baselines,freqs,coeffs = read_coeffs(calfile,na=na)
	# else:
	# 	baselines,freqs,coeffs = get_coeffs(calfile,na=na,freqs=n.array([.1,]))
	# 	coeffs=coeffs[:,0]
	
