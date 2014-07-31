#!/usr/bin/env python

"""

NAME: 
      pspec_sim_v2.py 
PURPOSE:
      -Generates a_lm values that are drawn from a Gaussian distribution with variance specified by C_l
      -Inverse spherical harmonic transform to get sky
      -Plots sky onto a Healpix Map for each frequency
EXAMPLE CALL:
      ./pspec_sim_v2.py 
AUTHOR:
      Carina Cheng

"""

import aipy
import numpy
import pylab
import pyfits
import matplotlib.pyplot as plt
import optparse
import os, sys
import healpy
import scipy
from scipy import integrate
from scipy.interpolate import interp1d

###################################
###    OPTIONS & PARAMETERS     ###
###################################

###options

o = optparse.OptionParser()
o.set_usage('vis_simulation_v2.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--nchan', dest='nchan', default=3, type='int',
             help='Number of channels in simulated data. Default is 3.')
o.add_option('--sfreq', dest='sfreq', default=0.1, type='float',
             help='Start frequency (GHz). Default is 0.1.')
o.add_option('--sdf', dest='sdf', default=.001, type='float',
             help='Channel spacing (GHz).  Default is .001')
o.add_option('--lmax', dest='lmax', default=5, type='int',
             help='Maximum l value. Default is 5.')
opts, args = o.parse_args(sys.argv[1:])

###cosmology parameters

nu21 = 1420.*10**6 #21cm frequency [Hz]
c = 3.*10**5 #[km/s]
H0 = 69.7 #Hubble's constant [km/s/Mpc]
omg_m = 0.28 
omg_lambda = 0.72 

###################################
###          FUNCTIONS          ###
###################################

###P_k function

def P_k(kmag, sigma=0.1, k0=0.1):

    return numpy.exp(-(kmag-k0)**2/(2*sigma**2)) 

###C_l function

def C_l(freq1, freq2, Pk_interp, l_vals): #freqs are entered in GHz
    
    f1 = freq1*10**9 #[Hz]
    f2 = freq2*10**9 #[Hz]
    z1 = (nu21/f1)-1
    z2 = (nu21/f2)-1

    one_over_Ez = lambda zprime:1/numpy.sqrt(omg_m*(1+zprime)**3+omg_lambda)
    Dc1 = (c/H0)*integrate.quad(one_over_Ez,0,z1)[0]
    Dc2 = (c/H0)*integrate.quad(one_over_Ez,0,z2)[0]

    C_ls = []

    for i in range(len(l_vals)):

        integral = lambda k:(2/numpy.pi)*Pk_interp(k)*(k**2)*scipy.special.sph_jn(l_vals[i],k*Dc1)[0][i]*scipy.special.sph_jn(l_vals[i],k*Dc2)[0][i]

        ans = integrate.quad(integral,numpy.min(k_data),numpy.max(k_data))[0]
        
        C_ls.append(ans)

    return C_ls

###Generate correlated random variables (a_lms)

def a_lm(cov_T):

    L = numpy.linalg.cholesky(cov_T)
    Lt = numpy.conj(numpy.transpose(L))

    z = numpy.random.normal(scale=1.0,size=len(cov_T))

    return numpy.einsum('ij,j',L,z) #gives correct column vector shape

###Build C_matrix (function of l and frequency combinations)

def C_matrix(freqs,Pk_interp,l_vals):

    matrix = numpy.zeros((len(l_vals),len(freqs),len(freqs)))
    
    for i in range(len(freqs)):      
        for j in range(len(freqs)):
            Cls = C_l(freqs[i],freqs[j],Pk_interp,l_vals)
            for k in range(len(Cls)):
                matrix[k][i][j] = Cls[k]

    return matrix
 
      
###################################
###     SIM CODE STARTS HERE    ###
###################################
"""
k_data = numpy.arange(0.01,0.2,0.01) #actual data points
Pk_data = P_k(k_data)
#plt.plot(k_data,Pk_data, 'k.')
#plt.show()

Pk_interp = interp1d(k_data, Pk_data, kind='linear')

l_vals = numpy.arange(0,opts.lmax,1) 
a_lms = aipy.healpix.Alm(opts.lmax,opts.lmax) #a_lm object

freqs = numpy.linspace(opts.sfreq,opts.sfreq+opts.sdf*opts.nchan,num=opts.nchan, endpoint=False) #array of frequencies

Cmatrix = C_matrix(freqs,Pk_interp,l_vals)
print Cmatrix
#print Cmatrix[0,0,0] #l=0,freq0 with freq0
#print Cmatrix[1,0,1] #l=1, freq0 with freq1

#for each freq, fill all a_lms and make sky map

dirname = 'pspec'+str(opts.nchan)
os.system('mkdir /Users/carinacheng/capo/ctc/images/pspecs/'+dirname)

for f in range(len(freqs)):

    for i in range(len(Cmatrix)):

        Csmall = Cmatrix[i]

        const = 0

        for j in range(l_vals[i]+1):

            alms = a_lm(Csmall) #alm for one lm-mode
            a_lms[l_vals[i],l_vals[i]-const] = float(alms[f])  #fill a_lm
            const += 1     
            sky_map = aipy.map.Map(nside=512)
            sky_map.from_alm(a_lms) #make sky map
            sky_map.to_fits('/Users/carinacheng/capo/ctc/images/pspecs/'+dirname+'/pspec1'+("%03i" % (f+1))+'.fits', clobber=True) 
"""


###################################
###       ADDITIONAL TESTS      ###
###################################

"""
###test a_lm function (recovering T)

cov_T = [[10,5,2],[5,15,5],[2,5,20]] #symmetric
#print numpy.linalg.eig(cov_T)[0] #check for positive eigenvalues

xxt_test = numpy.zeros_like(cov_T)
num = 1000

for i in range(num):

    x = a_lm(cov_T)
    x = numpy.reshape(x,(3,1))
    xt = numpy.conj(numpy.transpose(x))

    xxt = numpy.dot(x,xt)
    xxt_test = xxt_test + xxt

print xxt_test/num (should be T)
"""

"""
###make single frequency map

Cl = C_l(.15,.15,Pk_interp,l_vals)
#plt.plot(l_vals,Cl, 'k.')
#plt.show()

a_lms = aipy.healpix.Alm(lmax,mmax) #a_lm object

for i in range(len(Cl)):

    Csmall = [[Cl[i]]] #one block of C 

    l_val = l_vals[i]

    const = 0
    
    for j in range(l_val+1):

        alm = a_lm(Csmall)
        a_lms[l_val,l_val-const] = float(alm)
        const += 1

sky_map = aipy.map.Map(nside=512)
sky_map.from_alm(a_lms)

sky_map.to_fits('/Users/carinacheng/capo/ctc/images/test.fits', clobber=True)   
"""
"""
###plotting C_l(freq1,freq2)

freq1 = .1
freq2 = numpy.arange(.14,0.19,0.01) #[GHz]

k_data = numpy.arange(0.01,0.2,0.01) #actual data points
Pk_data = P_k(k_data)

Pk_interp = interp1d(k_data, Pk_data, kind='linear')

l_vals = numpy.arange(0,opts.lmax,1) 

for i in range(len(freq2)):

    freqs = [freq1,freq2[i]]

    Cmatrix = C_matrix(freqs,Pk_interp,l_vals)

    Cnu = []

    for j in range(len(l_vals)):

        Csmall = Cmatrix[j]
        Cnu.append(Csmall[0,1])

    plt.plot(l_vals,Cnu)

plt.show()
"""  






