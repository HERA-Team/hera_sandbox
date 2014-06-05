import pylab as plt 
import numpy as np

def linear_fit(xdata,ydata):
	if len(xdata.shape)==1: 
		xdata = np.reshape(xdata,(xdata.shape[0],1))
	if len(ydata.shape)==1: 
		ydata = np.reshape(ydata,(ydata.shape[0],1))
	nn = np.concatenate((np.ones_like(xdata),xdata),axis=1)
	MM = np.matrix(nn)
	bb = np.matrix(ydata)
	hh = MM.T*MM
	aa = hh.I*MM.T*bb
	B = aa[0,0] #intercept
	A = aa[1,0] #slope

	# calculate R^2 value
	fdata = A*xdata+B
	ss_tot = np.sum(np.absolute((ydata-np.average(ydata)))**2)
	ss_res = np.sum(np.absolute((ydata-fdata))**2)
	Rsq = 1-ss_res/ss_tot

	# calculate reduced chi^2 value
	sig = ss_tot/len(ydata)
	redchi = np.sum(np.absolute((ydata-fdata))**2/sig)/len(ydata)

	return A,B,Rsq,redchi
