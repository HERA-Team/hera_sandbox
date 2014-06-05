#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import homemade_plotting as hp

def get_coeffs(calfile='psa898_v002',na = 32,freqs = n.arange(.1,.2,.01)):
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
	            phs.shape = amp.shape = im.uv.shape #?????
	            amp = n.where(valid, amp, 0)
	            phs = n.where(valid, phs, 0)

	            dc_response = n.sum(amp*phs)/n.sum(amp) #this integrates the amplitude * fringe pattern; units mK?
	            jy_response = n.real(dc_response * 100 / C.pspec.jy2T(fq)) # jy2T converts flux density in jansky to temp in mK
	            print '\t',fq, dc_response, jy_response
	            coeffs[ll,kk] = dc_response
	            kk+=1
	        ll+=1
	#n.savez(filename,baselines=baselines,frequencies=freqs,coeffs=coeffs)
	#f = n.load(filename)
	#freqs = f[frequencies]
	n.savetxt('{0}_baselines.txt'.format(calfile),baselines)
	n.savetxt('{0}_frequencies.txt'.format(calfile),freqs)
	n.savetxt('{0}_coeffs.txt'.format(calfile),coeffs)
	return baselines,freqs,coeffs

def test_regress(gs=10,n_sig=5,na=32,readFromFile=False):
	if readFromFile:
		baselines = n.loadtxt('{0}_baselines.txt'.format('psa898_v002'))
		freqs = n.loadtxt('{0}_frequencies.txt'.format('psa898_v002'))
		coeffs = n.loadtxt('{0}_coeffs.txt'.format('psa898_v002'))
		baselines = baselines[:na]
		#freqs = freqs[:na]
		coeffs = coeffs[:na]
	else:
		baselines,freqs,coeffs = get_coeffs(na=na,freqs=n.array([.1,]))
		coeffs=coeffs[:,0]
	coeffs=coeffs.T
	VV = gs*coeffs + n.random.normal(loc=0.0,scale=n_sig,size=len(coeffs))
	print VV.shape
	VV = VV*n.exp(2*n.pi*1j*n.random.rand())
	print VV.shape
	#AA = n.array([coeffs,n.ones(len(coeffs))])
	# AA = n.vstack([coeffs, n.ones(len(coeffs))]).T
	# print AA
	# print VV.shape
	# print AA.shape
	#gs_recov,n_recov = n.linalg.lstsq(AA,VV)[0]
	gs_recov,n_recov,RR,redchi = hp.linear_fit(coeffs,VV)
	print "true gs = {0}\trecovered gs = {1}".format(gs,gs_recov)
	print "true n = {0}\trecovered n = {1}".format(0.0,n_recov)
	print "R = ",RR
	return gs_recov

def test_regress_vary_n():
	for jj in n.arange(100):
		gs_diff = n.zeros(20)
		n_sigs = n.logspace(-3,1,num=20)
		for ii,n_sig in enumerate(n_sigs):
			print ii
			#if ii==0 and jj==0:rFF=False 
			#else:rFF=True
			rFF=True
			gs_diff[ii] = test_regress(gs=1,n_sig=n_sig,readFromFile=rFF) - 1.
		p.plot(n_sigs,gs_diff)
		p.xscale('log')
	p.xlabel('Amplitude of noise relative to global signal\n(I.e. true global signal amplitude is 1)')
	p.ylabel('Difference between true and recovered global signal')
	#p.show()
	p.savefig('gs_diff_vs_n.pdf')
	p.clf()

def test_regress_vary_na():
	for jj in n.arange(100):
		nas = n.arange(3,32)
		gs_diff = n.zeros(len(nas))
		for ii,na in enumerate(nas):
			gs_diff[ii] = test_regress(gs=1,n_sig=.1,na=na,readFromFile=True) - 1.
		p.plot(nas,gs_diff)
	p.xlabel('Number of antenna')
	p.ylabel('Difference between true and recovered global signal')
	#p.show()
	p.savefig('gs_diff_vs_na.pdf')
	p.clf()
	

def test_regress_vary_bsln():
	for jj in n.arange(100):
		nas = n.arange(3,32)
		gs_diff = n.zeros(len(nas))
		for ii,na in enumerate(nas):
			gs_diff[ii] = test_regress(gs=1,n_sig=.1,na=na,readFromFile=True) - 1.
		p.plot(nas*(nas-1)/2,gs_diff)
	p.xlabel('Number of baselines')
	p.ylabel('Difference between true and recovered global signal')
	#p.show()
	p.savefig('gs_diff_vs_bsln.pdf')
	p.clf()


if __name__=='__main__': 
	#get_coeffs(na=32,freqs=n.array([.1,]))
	test_regress_vary_n()
	test_regress_vary_na()
	test_regress_vary_bsln()


