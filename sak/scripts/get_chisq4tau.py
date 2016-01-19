import pylab, aipy, sys, glob, numpy as np
from omnical import calibration_omni as omni

def ij2bl(i,j): 
	if i<=j: return (j+1)*j/2 + i
	if i>j: return (i+1)*i/2 + j 

nu_arr = np.linspace(.1, .2, num=203) #GHz

def chisq4tau(filename, tau_lim=30., delta_tau=0.01, save_npz=True):
	"""
	define array of possible tau values [always symmetric about zero] 
	and impose freqency coupling.
	
	Use these to calculate the chi squared statistic:
	
	X^2(nu,t) = sum_{baselines b} V^xy_b - V^{yx}_b exp(-2 pi i nu tau)
	
	for a sweep of taus; i.e. the x-to-y delay/angle that minimzes Stokes V
	across the array.
	
	Arguments:
	-- tau_lim: sweep of taus runs (-tau_lim, taulim)
	-- delta_tau: spacing of the sweep
	-- save_npz: saves an npz of the chi square for each tau, time and frequency
	
	Returns: an array of chi squared values for every tau, time and frequency.
	The minimum of the tau axis will give the required minimization for V.
	"""
	
	print '==== '+filename+' ===='
	
	dataxy,txy,timingxy,lstxy,flagsxy = omni.importuvs([filename],{'xy':-7}, timingTolerance=1e9, init_mem=.5e9)
	datayx,tyx,timingyx,lstyx,flagsyx = omni.importuvs([filename],{'yx':-8}, timingTolerance=1e9, init_mem=.5e9)
	
	global tau_arr
	global angle_arr
	
	tau_arr = np.arange(-1*tau_lim, tau_lim, delta_tau) #nanseconds
	angle_arr = np.exp(-2*np.pi*1.j*np.outer(nu_arr,tau_arr))
	
	#calculate the chi for each baseline; sum over the array
	chisquared = np.zeros( (tau_arr.shape[0], dataxy.shape[1], nu_arr.shape[0]) ) #tau x time x freq
	 
	###XXX TODO: VECTORIZE THIS BABY
	###XXX update: I get memory errors when I try to create a baseline axis by
	###### outer-producting with angle_arr. How to fix?
	
	print 'Begin stupid tau/bl loop (vectorizing is giving Mem Errs)' 
	for i in range(angle_arr.shape[1]): # for each tau
		calc = 0
		for bl in range(dataxy.shape[3]): #for each baseline
			XY = dataxy[0,:,:,bl]
			YX = datayx[0,:,:,bl]
			calc += pow(np.absolute(XY - YX*angle_arr[:,i]),2.)
		chisquared[i,:,:]=calc
	
	a=filename.split('.')
	a[0]='tau'
	a='.'.join(a)+'.npz'
	
	if save_npz: np.savez(a, chisquared=chisquared, angle=angle_arr, tau=tau_arr)
	return chisquared

""" EXAMPLE IMPLEMENTATION
fList = []
for f in glob.glob('raw_data/*uvcRRECXR'): fList.append(f)
fList = sorted(fList)

for f in fList: chisq4tau(f, tau_lim=30., delta_tau=0.1)
"""
