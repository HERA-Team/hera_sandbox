import pylab as plt 
import numpy as np
from scipy import optimize

#  _____                 _   _                 
# |  ___|   _ _ __   ___| |_(_) ___  _ __  ___ 
# | |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
# |  _|| |_| | | | | (__| |_| | (_) | | | \__ \
# |_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
                                             
def vdot(M,v):
    """
    Dots a matrix and a vector and returns a vector of shape (n,)
    Use because np.dot(M,v) will return a vector (1,n) 
    """
    assert len(v.shape)==1
    vec = np.dot(M,v)
    # print 'uf vec org ',vec.shape
    if len(vec.shape)==2 and vec.shape[1]==1: vec = np.resize(vec,(vec.shape[0],))
    if len(vec.shape)==2 and vec.shape[0]==1: vec = np.resize(vec,(vec.shape[1],))
    # print 'uf M ',M.shape
    # print 'uf v ',v.shape
    # print 'uf vec ',vec.shape
    assert len(vec.shape)==1
    return vec 

def gaussian(sig,xpoints,ypoints,x0=0.0,y0=0.0):
    """
    Returns a gaussian distribution of width sig, centered on x0,y0 
    for the data points in xpoints and ypoints
    """
    gauss = 1/(2*np.pi*sig*sig)*np.exp(-((xpoints-x0)**2+(ypoints-y0)**2)/(2*sig*sig))
    return gauss

def rand_from_covar(C):
    """
    Takes in a covariance matrix C = < x x_dag > and generates a random vector x 
    that is consistent with C.

    We do this by using a Cholesky decomposition of C = L L_dag where L is lower 
    triangular. Then the x we want is x = L z where z has uncorrelated random 
    elements with variance 1. Therefore, < x x_dag > = L < z z_dag > L_dag = L L_dag = C 
    """
    L = np.linalg.cholesky(C)
    z = np.random.randn(C.shape[0])
    x = np.dot(L,z)
    x = np.array(x)
    return x 

def projection_matrix(Z):
    Z = np.matrix(Z)
    #print 'Z = \n',Z
    PP = np.identity(Z.shape[0]) - np.dot(Z,Z.H) #checked by hand
    return PP 

def pseudo_inverse(MM,num_remov=1):
    """
    Computes a matrix pseudo inverse based on equation A4 in Max Tegmark's 
    1997 paper "How to measure the CMB power spectra without losing information"
    """
    eta = np.average(MM)
    #print 'eta = ',eta
    MM = np.matrix(MM)
    #print 'MM = \n',MM
    eig_vals, eig_vecs = np.linalg.eig(MM) # checked
    #print 'eig_vecs.shape = ',eig_vecs.shape
    sorted_ind = np.argsort(np.absolute(eig_vals))
    sorted_vecs = eig_vecs[:,sorted_ind]
    sorted_vals = eig_vals[sorted_ind]
    #print 'eig_vals = \n',sorted_vals
    #print 'eig vecs = \n',sorted_vecs
    if num_remov==None:
        Z = None
        for kk in range(sorted_vals.shape[0]):
            #print np.min(sorted_vals[kk:])/np.max(sorted_vals[kk:])
            if np.min(sorted_vals[kk:])/np.max(sorted_vals[kk:])<1.01*10**-4:
                #print 'removing lambda ',sorted_vals[kk]
                if Z==None:
                    Z = sorted_vecs[:,kk]
                else:
                    Z = np.hstack((Z,sorted_vecs[:,kk]))
            else:
                #print 'broke at lambda ',sorted_vals[kk]
                break       
    else:
        Z = sorted_vecs[:,:num_remov] #checked
    #Z = np.matrix(Z)
    #print 'Z = \n',Z
    PP = projection_matrix(Z) #np.identity(eig_vals.shape[0]) - np.dot(Z,Z.H) #checked by hand
    #print 'PP = \n',PP
    MM_tilde = np.dot(PP,np.dot(MM,PP.H)) #checked
    #print 'MM_tilde = \n',MM_tilde
    AA = MM_tilde + eta*np.dot(Z,Z.H) 
    #print 'AA = \n',AA
    AAinv = np.linalg.inv(AA)
    #np.set_printoptions(threshold='nan')
    #print '\nAAinv = \n',AAinv
    MM_inv = np.dot(PP,np.dot(AAinv,PP.H)) 
    #print 'MM_inv = \n',MM_inv

    #test the inverse
    # v1 = sorted_vecs[:,-1]#+ sorted_vecs[:,-2]
    # v2 = Z[:,0]
    # dotprod = np.dot(v1.H,v2)
    # print dotprod[0,0]
    # v = v1 - np.dot(v1.H,v2)[0,0]*v2
    # print np.dot(v.H,v2)
    # print 'test vector = \n',v
    # vp = np.dot(np.dot(MM_inv,MM_tilde),np.array(v))
    # print 'recovered vector = \n',vp

    return MM_inv


#  _____ _ _   _   _             
# |  ___(_) |_| |_(_)_ __   __ _ 
# | |_  | | __| __| | '_ \ / _` |
# |  _| | | |_| |_| | | | | (_| |
# |_|   |_|\__|\__|_|_| |_|\__, |
#                          |___/ 

def general_lstsq_fit_with_err(xdata,ydata,Q,noiseCovar,pseudo=False):
    """
    This function takes in x and y data, Q, and a full noise covariance 
    matrix. The function returns <\hat x> and the covariance matrix of 
    the fit.
    """
    xdata = np.array(xdata)
    Q = np.matrix(Q)
    Ninv = np.matrix(noiseCovar).I
    AA = Q.H*Ninv*Q
    if pseudo: AAinv = np.linalg.pinv(AA)
    elif isinstance(pseudo, int): AAinv = pseudo_inverse(AA,num_remov=pseudo)
    else: AAinv = AA.I
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

