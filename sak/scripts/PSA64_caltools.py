import pylab
import numpy as np
import omnical as omni
import aipy, glob, sys, os

"""

	

"""

def ij2bl(i,j): 
	if i<=j: return (j+1)*j/2 + i
	if i>j: return (i+1)*i/2 + j 
	
def MLE_all(target_folder,BLINDEX,
reduninfo,
breaklen=100,
plot=True,save=True,show=False):
	
	BLINDEX = int(BLINDEX)
	
	print 'Target folder:',target_folder
	print 'Getting all uvcRREcCPSA128O files'
	filelist = []
	for f in glob.glob(target_folder+'/*.uvcRREcCXPSA128O'): filelist.append(f)
	filelist = sorted(filelist)
	
	assert(len(filelist)>=1)
	print 'Number of files:',len(filelist)
	c=0
	for F in filelist:
		if c>breaklen: break #for testing purposes. Default=100 i.e. more files than I'd probably use this on at any given time.
		print 'File',c
		c+=1
		#print 'File read-in'
		dataxy,txy,timingxy,lstxy,flagsxy = omni.calibration_omni.importuvs([F],{'xy':-7},timingTolerance=1e9,init_mem=.5e9)
		datayx,tyx,timingyx,lstyx,flagsyx = omni.calibration_omni.importuvs([F],{'yx':-8},timingTolerance=1e9,init_mem=.5e9)
	
		info = omni.calibration_omni.read_redundantinfo(reduninfo)
		unit = np.argsort(np.linalg.norm(info['ubl'],axis=-1))[BLINDEX] #<--redun info index
		#I can change the [<number>] @ the end of one_unit for different redundant baseline lengths
		unit_bls = info['subsetant'][info['ublindex'][unit][:,:2].astype('int')] #ij values of all good antennae
		if c==1: print F
		if c==1: print 'Number of this-length baselines:',unit_bls.shape[0]
		if c==1: print 'Baseline type',unit_bls[0]
		
		tau_xy_arr = [] #<--- tau_xy per baseline
		tau_spec_arr=[] #<--- tau_xy per frequency
		
		for bl in unit_bls:
			i = bl[0]
			j = bl[1]
			if i==50 or j==50: continue #the one offending antenna that only gets flagged sometimes over night 2456240
			xy = dataxy[0,:,:,ij2bl(i,j)]
			yx = datayx[0,:,:,ij2bl(i,j)]
			
			numerator = sum(xy*yx.conjugate())
			denominator = sum(np.absolute(yx)**2.)
			"""
			Maximum likelihood estimator for tau_xy is
			
			V_xy V_yx^* / | V_yx |^2 = exp(-2*pi*i*nu*tau_xy)
			
			so taking the real part gives
			
			tau_xy = np.arccos(Re(V_xy V_yx^* / | V_yx |^2))/2*pi*nu
			"""
			
			tau_xy = np.arccos((numerator/denominator).real)/(2.*np.pi*np.linspace(.1,.2,203))
			tau_spec_arr.append(tau_xy)
			
			TAU = np.nanmean(tau_xy)
			tau_xy_arr.append(TAU)
			
			if plot: pylab.plot(tau_xy)
	tau_xy_arr = np.array(tau_xy_arr)
	print 'AVERAGE TAU_XY = ',np.average(tau_xy_arr),'+/-',np.std(tau_xy_arr)
			
	if plot:
		pylab.xlabel('Frequency bin',size=15)
		pylab.ylabel(r'$\tau$ (ns)',size=20)
		if save: pylab.savefig(target_folder+'/all_taus_vs_freq_BL%s.png'%str(BLINDEX))
		if show: pylab.show()
		pylab.close()
		
		try:
			pylab.hist(tau_xy_arr,bins=10)
			pylab.xlabel(r'$\tau$ (ns)',size=20)
			pylab.ylabel(r'Count',size=15)
			if save: pylab.savefig(target_folder+'/tau_hist_BL%s.png'%str(BLINDEX))
			if show: pylab.show()
			pylab.close()
		except: print "Trouble getting a histogram to work"

		pylab.imshow(tau_spec_arr,aspect='auto',interpolation='None')
		pylab.ylabel('Integration',size=15)
		pylab.xlabel('Frequency Bin',size=15)
		pylab.colorbar()
		if save: pylab.savefig(target_folder+'/tau_waterfall_BL%s.png'%str(BLINDEX))
		if show: pylab.show()
		pylab.close()
	
		pylab.plot(np.linspace(.1,.2,203)*1000.,np.nanmean(tau_spec_arr,axis=0))
		pylab.xlabel('Frequency (MHz)',size=15)
		pylab.ylabel(r'$\tau$ (ns)',size=15)
		if save: pylab.savefig(target_folder+'/mean_taus_vs_freq_BL%s.png'%str(BLINDEX))
		if show: pylab.show()
		pylab.close()
		
	print 'Calculations -> '+target_folder+'/%s_BL%s.npz'%(target_folder,str(BLINDEX))
	np.savez(target_folder+'/%s_BL%s.npz'%(target_folder,str(BLINDEX)),tau_waterfall=tau_spec_arr,mean_tau_per_BL=tau_xy_arr)
		
