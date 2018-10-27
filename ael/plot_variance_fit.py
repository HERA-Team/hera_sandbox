import numpy as np
import sys
import pylab as pl
from scipy.stats import linregress
from scipy.optimize import least_squares

for fname in sys.argv[1:]:

    f = np.load(fname)

    variances = f['variances']
    Nsamps = f['avg_lens_nsamp']
    avg_times = f['avg_lens_time']
    stds = f['var_stds']
    
    xvar = f['avg_lens_nsamp']
    xvar = xvar.astype(float)
    
    ## Using linregress:
    #linfit = linregress(np.log10(xvar), np.log10(variances))
    #b0 = linfit.slope
    #a0 = 10**(linfit.intercept)
    #fitted = a0*xvar**b0
    #print('a = {:.3f}, b = {:.3f}'.format(a0,b0))
    
    def fit_func(x, a, b):
        """
            a = Variance * averaging time
            b = power
            x = averaging time
        """
        #return variances[0]/xvar[0]**b * (x.astype(float))**b 
        return a*(x.astype(float))**(b)
    
    def residual(fit, time_arr, var, weight, c):
        a, b = fit
    #    c = int(c)
        return (np.log10(var[c:]) - np.log10(fit_func(time_arr[c:], a, b)))#*np.sqrt(float(c))#*weight
        #return (var - fit_func(time_arr, a,b))#*weight
    
    
    
    #TODO Define weights properly
    weights = f['nsamps'].astype(float)          #The more samples of variance, the more robust.
    
    guess = (np.max(variances)*xvar[0]**(-1), -1)
    #costs = []
    #cutmax = 500
    #for cutoff in range(0, cutmax):
    #    lsq_res = least_squares(residual, guess, args=(xvar.astype(float), variances, weights, cutoff), method='trf')#, xtol=1e-15, ftol=1e-15, gtol=1e-15)
    #    print cutoff, lsq_res.cost
    #    costs.append(lsq_res.cost)
    
    cutoff=0
    lsq_res = least_squares(residual, guess, args=(xvar.astype(float), variances, weights, cutoff), method='trf')
    
    fitted = fit_func(xvar, *lsq_res.x)
    print(np.mean(lsq_res.fun), lsq_res.x)
    print("fname: {}, Percent error with {:.1f} hours of averaging: {:.4f}%".format(fname, avg_times[-1] / 60., np.abs(np.sqrt(variances[-1])/np.mean(f['avgs'][-1]) )*100.0))
    
    fig, ax = pl.subplots(1,1)
    ax.fill_between(avg_times, variances-stds, variances+stds, alpha=0.75, color='0.75')
    ax.plot(avg_times, fitted, label='Lsq fit')
    ax.plot(avg_times, variances, label='Variance')
    ax.set_title(fname)
    ax.set_yscale('log')
    #ax.set_xscale('log')
    pl.legend()

pl.show()
