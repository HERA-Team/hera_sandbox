import pylab as plt 
import numpy as np
from scipy import optimize


def general_lstsq_fit_with_err(xdata,ydata,Q,noiseCovar):
    """
    This function takes in x and y data, Q, and a full noise covariance 
    matrix. The function returns <\hat x> and the covariance matrix of 
    the fit.
    """
    xdata = np.array(xdata)
    Q = np.matrix(Q)
    Ninv = np.matrix(noiseCovar).I
    AA = Q.H*Ninv*Q
    AAinv = AA.I
    params = AAinv*Q.H*Ninv*np.matrix(ydata) # params should be 1 by nparam
    return np.array(params), np.array(AAinv)

def line_thru_origin_lstsq_fit_with_err(xdata,ydata,noiseCovar):
    xdata = np.array(xdata)
    gridnum = len(xdata)
    Q = np.matrix(np.zeros([gridnum,1],dtype='complex'))
    for ii in range(gridnum):
        Q[ii] = xdata[ii][0]
    params, variance = general_lstsq_fit_with_err(xdata,ydata,Q,noiseCovar)
    return params, variance

def poly_lstsq_fit_with_err(xdata,ydata,noiseCovar,nparam):
    """
    This function takes in x and y data, a full noise covariance matrix,
    and the degree of the polynomial it will fit (nparam). It fits to a 
    polynomial y = a0 + a1 * x + a2 * x^2 + a3 * x^3 + ... + an * x^n. 
    The function returns the arrays of a's and the covariance matrix of 
    the fit.
    """
    xdata = np.array(xdata)
    gridnum = len(xdata)
    Q = np.matrix(np.zeros([gridnum,nparam],dtype='complex'))
    for ii in range(gridnum):
        for jj in range(nparam):
            Q[ii,jj] = xdata[ii][0]**jj
    params, variance = general_lstsq_fit_with_err(xdata,ydata,Q,noiseCovar)
    return params, variance
    # Ninv = np.matrix(noiseCovar).I
    # AA = Q.H*Ninv*Q
    # AAinv = AA.I
    # params = AAinv*Q.H*Ninv*np.matrix(ydata) # params should be 1 by nparam
    # return np.array(params), np.array(AAinv)

def linear_fit(xdata,ydata):
    if len(xdata.shape)==1: 
        xdata = np.reshape(xdata,(xdata.shape[0],1))
    if len(ydata.shape)==1: 
        ydata = np.reshape(ydata,(ydata.shape[0],1))
    nn = np.concatenate((np.ones_like(xdata),xdata),axis=1)
    #print nn
    MM = np.matrix(nn)
    bb = np.matrix(ydata)
    hh = MM.H*MM
    aa = hh.I*MM.H*bb
    B = aa[0,0] #intercept
    A = aa[1,0] #slope

    # calculate R^2 value
    fdata = A*xdata+B
    ss_tot = np.sum(np.absolute((ydata-np.average(ydata)))**2)
    ss_res = np.sum(np.absolute((ydata-fdata))**2)
    Rsq = 1-ss_res/ss_tot

    # calculate reduced chi^2 value
    sig = ss_tot/len(ydata)
    redchi = np.sum(np.absolute((ydata-fdata))**2/sig)/(len(ydata) - 2 - 1)

    return A,B,Rsq,redchi

def linear_fit_with_err(xdata,ydata,noiseVariance):
    if len(xdata.shape)==1:
        xdata = np.reshape(xdata,(xdata.shape[0],1))
    if len(ydata.shape)==1:
        ydata = np.reshape(ydata,(ydata.shape[0],1))
    nn = np.concatenate((np.ones_like(xdata),xdata),axis=1)
    MM = np.matrix(nn)
    bb = np.matrix(ydata)
    hh = MM.H*MM
    aa = hh.I*MM.H*bb
    B = aa[0,0] #intercept
    A = aa[1,0] #slope
    
    # calculate R^2 value
    fdata = A*xdata+B
    ss_tot = np.sum(np.absolute((ydata-np.average(ydata)))**2)
    ss_res = np.sum(np.absolute((ydata-fdata))**2)
    Rsq = 1-ss_res/ss_tot
    
    # calculate reduced chi^2 value
    sig = ss_tot/len(ydata)
    redchi = np.sum(np.absolute((ydata-fdata))**2/sig)/(len(ydata) - 2 - 1)

    # Calculate the error bars using [A^dag N^-1 A]^-1
    err = np.sqrt(noiseVariance*np.linalg.inv(hh)[1,1])
    return A,B,Rsq,redchi,err

def linear_fit_new(xdata,ydata):
    """
    This doesn't work yet.
    """
    print xdata.shape,ydata.shape
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y: (y - fitfunc(p, x))

    prm0 = np.array([1.0, -1.0]) #initial guess for parameters
    out = optimize.leastsq(errfunc, prm0, args=(xdata, ydata),full_output=1)
    intercept,slope = out[0] 
    covar = out[1] 

    # calculate reduced chi^2 value
    fdata = slope*xdata+intercept
    sig = np.sum(np.absolute((ydata-np.average(ydata)))**2)/len(ydata)
    redchi = np.sum(np.absolute((ydata-fdata))**2/sig)/(len(ydata)-2-1)

    return slope, intercept, redchi

def line_thru_origin_fit(xdata,ydata):
    """
    I'm not sure if this works yet
    """
    # fitfunc = lambda p, x: p*x 
    # errfunc = lambda p, x, y: (y - fitfunc(p, x))
    # p0 = 1.0 #initial guess for parameter
    # out = optimize.leastsq(errfunc, np.array([p0,]), args=(xdata, ydata),full_output=1)
    # pf = out[0]
    pf, _, _, _  = np.linalg.lstsq(xdata,ydata)
    pf = pf[0][0]
    print pf

    # calculate reduced chi^2 value
    fdata = pf*xdata
    sig = np.sum(np.absolute((ydata-np.average(ydata)))**2)/len(ydata)
    redchi = np.sum(np.absolute((ydata-fdata))**2/sig)/(len(ydata)-1-1)

    return pf,redchi

def line_thru_origin_with_err(xdata,ydata,noiseVariance):
    """
    I don't think this works yet
    """
    if len(xdata.shape)==1:
        xdata = np.reshape(xdata,(xdata.shape[0],1))
    if len(ydata.shape)==1:
        ydata = np.reshape(ydata,(ydata.shape[0],1))
    #nn = np.concatenate((np.ones_like(xdata),xdata),axis=1)
    MM = np.matrix(xdata)
    bb = np.matrix(ydata)
    hh = MM.H*MM
    A = hh.I*MM.H*bb
    A = A.item(0)
    print A
    #B = aa[0,0] #intercept
    #A = aa[1,0] #slope
    
    # calculate R^2 value
    fdata = A*xdata
    ss_tot = np.sum(np.absolute((ydata-np.average(ydata)))**2)
    ss_res = np.sum(np.absolute((ydata-fdata))**2)
    Rsq = 1-ss_res/ss_tot
    
    # calculate reduced chi^2 value
    sig = ss_tot/len(ydata)
    redchi = np.sum(np.absolute((ydata-fdata))**2/sig)/(len(ydata) - 2 - 1)

    # Calculate the error bars using [A^dag N^-1 A]^-1
    print noiseVariance*np.linalg.inv(hh)
    err = np.sqrt(noiseVariance*np.linalg.inv(hh).item(0))
    return A,Rsq,redchi,err

def power_law_lstsq_fit(xdata,ydata):
    """
    This power-law fitting is done by first converting
    to a linear equation and then fitting to a straight line.
     y = a * x^b
     log(y) = log(a) + b*log(x) 
    """ 
    logx = np.log10(xdata)
    logy = np.log10(ydata)
    #logyerr = yerr / ydata
 
    # define our (line) fitting function
    fitfunc = lambda p, x: p[0] + p[1] * x
    errfunc = lambda p, x, y: (y - fitfunc(p, x))
    #errfunc = lambda p, x, y, err: (y - fitfunc(p, x)) / err
  
    prm0 = [1.0, -1.0] #initial guess for parameters
    out = optimize.leastsq(errfunc, prm0, args=(logx, logy),full_output=1)#, logyerr), full_output=1)
    prm = out[0] 
    covar = out[1] 
    #print prm
    #print covar

    index = prm[1]
    amp = 10.0**prm[0]
    indexErr = np.sqrt(covar[0][0])
    ampErr = np.sqrt(covar[1][1])*amp

    # calculate reduced chi^2 value
    fdata = amp*xdata**index
    sig = np.sum(np.absolute((ydata-np.average(ydata)))**2)/len(ydata)
    redchi = np.sum(np.absolute((ydata-fdata))**2/sig)/(len(ydata)-2-1)

    return amp, index, ampErr, indexErr, redchi

def legendre_lstsq_fit(xdata,ydata,deg):
    coefs = np.polynomial.legendre.legfit(xdata, ydata, deg)
    fdata = np.polynomial.legendre.legval(xdata,coefs)
    # calculate reduced chi^2 value
    sig = np.sum(np.absolute((ydata-np.average(ydata)))**2)/len(ydata)
    redchi = np.sum(np.absolute((ydata-fdata))**2/sig)/(len(ydata)-deg-1)
    return coefs, redchi

