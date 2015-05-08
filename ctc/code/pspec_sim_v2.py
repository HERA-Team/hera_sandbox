#!/usr/bin/env python

"""

NAME: 
      pspec_sim_v2.py 
PURPOSE:
      -Generates a_lm values that are drawn from a Gaussian distribution with variance specified by C_l
      -Inverse spherical harmonic transform to get sky
      -Plots sky onto a Healpix Map for each frequency
EXAMPLE CALL:
      ./pspec_sim_v2.py --nchan 50 --sdf 0.001 --lmax 200
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
#import healpy
import scipy
#import mpmath
from scipy import integrate
from scipy.interpolate import interp1d

###################################
###    OPTIONS & PARAMETERS     ###
###################################

###options

o = optparse.OptionParser()
o.set_usage('pspec_sim_v2.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--nchan', dest='nchan', default=203, type='int',
             help='Number of channels in simulated data. Default is 203.')
o.add_option('--sfreq', dest='sfreq', default=0.1, type='float',
             help='Start frequency (GHz). Default is 0.1.')
o.add_option('--sdf', dest='sdf', default=0.1/203, type='float',
             help='Channel spacing (GHz).  Default is 0.1/203.')
o.add_option('--lmax', dest='lmax', default=5, type='int',
             help='Maximum l value. Default is 5.')
opts, args = o.parse_args(sys.argv[1:])

###cosmology parameters

nu21 = 1420.*10**6 #21cm frequency [Hz]
c = 3.*10**5 #[km/s]
H0 = 100 #69.7 #Hubble's constant [km/s/Mpc] #can use 100 here to include h units
omg_m = 0.28 
omg_lambda = 0.72 

###################################
###          FUNCTIONS          ###
###################################

###P_k function

def P_k(kmag, sigma=0.01, k0=0.02):

    return kmag*0+1. #flat P(k)
    #return numpy.exp(-(kmag-k0)**2/(2*sigma**2)) #Gaussian
    #return numpy.exp(-(kmag-k0)**2/(2*sigma**2))+numpy.exp(-(kmag-k1)**2/(2*sigma**2)) #2 Gaussians


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

        ###integral form 1, k bins
        #integral = lambda k:(2/numpy.pi)*Pk_interp(k)*(k**2)*scipy.special.sph_jn(l_vals[i],k*Dc1)[0][l_vals[i]]*scipy.special.sph_jn(l_vals[i],k*Dc2)[0][l_vals[i]] 
        #ans1 = integrate.romberg(integral, numpy.min(k_data), numpy.max(k_data),tol=10**-18,rtol=10**-18,divmax=15)

        ###summation instead of integral
        ans2=0
        for kk in range(len(k_data)):
            #val = (2/numpy.pi)*Pk_interp(k_data[kk])*(k_data[kk]**2)*scipy.special.sph_jn(l_vals[i],k_data[kk]*Dc1)[0][l_vals[i]]*scipy.special.sph_jn(l_vals[i],k_data[kk]*Dc2)[0][l_vals[i]]
            
            J_sub = (1/2.)*(2*l_vals[i]+1) #coding out spherical bessel functions (sph_jn doesn't work for high k)
            x1 = k_data[kk]*Dc1
            x2 = k_data[kk]*Dc2
            bessel1 = numpy.sqrt(numpy.pi/2)*scipy.special.jn(J_sub,x1)/(numpy.sqrt(x1))
            bessel2 = numpy.sqrt(numpy.pi/2)*scipy.special.jn(J_sub,x2)/(numpy.sqrt(x2))
            val = (2/numpy.pi)*Pk_interp(k_data[kk])*(k_data[kk]**2)*bessel1*bessel2
            
            ans2+=val
        ans2*=(k_data[1]-k_data[0])

        ###integral form 1, ln k bins
        #integral = lambda u:(2/numpy.pi)*Pk_interp(numpy.exp(u))*(numpy.exp(3*u))*scipy.special.sph_jn(l_vals[i],numpy.exp(u)*Dc1)[0][i]*scipy.special.sph_jn(l_vals[i],numpy.exp(u)*Dc2)[0][i] 
        #ans = integrate.quad(integral,numpy.min(numpy.log(k_data)),numpy.max(numpy.log(k_data)))[0]

        ###integral form 2
        #k_perp = l_vals[i]/((Dc1+Dc2)/2)
        #max_int = numpy.sqrt(numpy.max(k_data)**2-k_perp**2)
        #integral = lambda k_par: (1/(numpy.pi*Dc1*Dc2))*Pk_interp(numpy.sqrt(k_perp**2+k_par**2))*numpy.cos(k_par*(Dc1-Dc2)) 
        #ans = integrate.quad(integral, 0, max_int)[0]

        C_ls.append(ans2)

    return C_ls#, Dc1, Dc2 

###Generate correlated random variables (a_lms)

def a_lm(cov_T):

    za = numpy.random.normal(scale=1/numpy.sqrt(2),size=len(cov_T))
    zb = numpy.random.normal(scale=1/numpy.sqrt(2),size=len(cov_T))
    z = za+1j*zb

    #cholesky decomposition

    #L = numpy.linalg.cholesky(cov_T)

    #return numpy.einsum('ij,j',L,z) #gives correct column vector shape
    
    #alternate way: eigendecomposition

    evals, evecs = numpy.linalg.eig(cov_T)
    V = numpy.transpose(evecs)      
    Lambda = numpy.identity(len(cov_T))
    Lambda = Lambda*evals
    Lambda_one_half=numpy.lib.scimath.sqrt(Lambda) #returns complex number if Lambda is negative & real
    #print numpy.dot(numpy.dot(numpy.linalg.inv(V),Lambda),V) #good check (should get T)
    X = numpy.dot(Lambda_one_half,V)
    Xt = numpy.dot(numpy.linalg.inv(V),Lambda_one_half)
    #print numpy.dot(Xt,X) #good check (should get T)

    #return numpy.einsum('ij,j',X,z)
    return numpy.dot(z,X)
    
###Build C_matrix (function of l and frequency combinations)

def C_matrix(freqs,Pk_interp,l_vals):

    matrix = numpy.zeros((len(l_vals),len(freqs),len(freqs)))
    
    for i in range(len(freqs)):      
        for j in range(len(freqs)):
            if j == 0:
                print '   working on freq '+str(i+1)+'/'+str(len(freqs))
            Cls = C_l(freqs[i],freqs[j],Pk_interp,l_vals)
            for k in range(len(Cls)):
                matrix[k][i][j] = Cls[k]

    return matrix
 
      
###################################
###     PARAMETERS TO CHANGE    ###
###################################

###if using P(k) function
"""
k_data = numpy.arange(0.01,0.5,0.01) #actual data points
Pk_data = P_k(k_data)

Pk_interp = interp1d(k_data, Pk_data, kind='linear')
"""
###if using P(k) points

#delta_sq = numpy.array([0.9413099, 2.047375, 5.35123, 9.884925, 19.12472, 24.78529, 28.9412, 29.81639, 29.235, 29.36708, 34.75652, 54.99383, 87.24504])
#ks = numpy.array([0.0281599, 0.04730816, 0.06677343, 0.09664913, 0.144889, 0.21846043, 0.327466, 0.49140543, 0.73756543, 1.10653729, 1.65998857, 2.39712429, 3.13152])
#Pk_data = delta_sq*(2*numpy.pi**2)/(ks**3)

#Pk_interp = interp1d(ks, Pk_data, kind='linear')

k_data = numpy.arange(0.001,0.5,0.01) #edit this
delta_sq = 0.000505 #K^2
Pk_data = delta_sq*(2*numpy.pi**2)/(k_data**3)
Pk_interp = interp1d(k_data,Pk_data,kind='linear')
#print k_data,Pk_data
###################################
###     SIM CODE STARTS HERE    ###
###################################

l_vals = numpy.arange(0,opts.lmax,1) 
a_lms = aipy.healpix.Alm(opts.lmax,opts.lmax) #a_lm object

freqs = numpy.linspace(opts.sfreq,opts.sfreq+opts.sdf*opts.nchan,num=opts.nchan, endpoint=False) #array of frequencies

print 'generating covariance matrix...'

Cmatrix = C_matrix(freqs,Pk_interp,l_vals)
#print Cmatrix[0,0,0] #l=0,freq0 with freq0
#print Cmatrix[1,0,1] #l=1, freq0 with freq1

#for each freq, fill all a_lms and make sky map

dirname = 'pspec'+str(opts.nchan)+'lmax'+str(opts.lmax)
os.system('mkdir ~/capo/ctc/images/pspecs/'+dirname)

alms_all = []

for i in range(len(Cmatrix)):
        
    #print 'getting Csmall for l = ' + str(i+1)

    Csmall = Cmatrix[i] #one block of C
    #noise = numpy.identity(len(Csmall))
    #noise *= 10**-20
    #Csmall = Csmall + noise
    '''
    #plotting eigenvectors
    evals = numpy.linalg.eig(Csmall)[0]
    evecs = numpy.linalg.eig(Csmall)[1]
    plt.figure(figsize = (12,12))
    plt.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.1, top = 0.9, wspace = 0.3, hspace = 0.3)
    for e in range(len(evals)):           
        plt.subplot(5,4,e+1)
        plt.plot(freqs,evecs[:,e],'k-')
        plt.title('eigenvalue = ' + str(evals[e]), fontsize=8)  
    plt.show()
    '''
    const = 0

    for j in range(l_vals[i]+1):

        alms = a_lm(Csmall) #alm for one lm-mode
        if (l_vals[i]-const) == 0: #real value for m=0
            alms = numpy.real(alms)
        alms_all.append(alms)
        const +=1

#print alms_all #alms_all[0] is lm-mode 00, alms_all[1] is lm-mode 10, alms_all[2] is lm-mode 11, etc. Each item in alms_all contains alms for that specific mode for each frequency.

for f in range(len(freqs)):

    print str(f+1) + '/' + str(len(freqs)) + ': making sky map with freq = ' + str(freqs[f]) + ' GHz'    

    lconst = 0
    mconst = 0

    for j in range(len(alms_all)):
        
        a_lms[lconst,mconst] = alms_all[j][f]

        if lconst == 0 and mconst == 0:
            lconst += 1
            mconst += 1
        elif lconst == mconst:
            mconst -= 1
        elif mconst == 0:
            lconst += 1
            mconst = lconst
        else:
            mconst -= 1
    
    sky_map = aipy.map.Map(nside=512)
    sky_map.from_alm(a_lms) #make sky map
    sky_map.to_fits('/home/cacheng/capo/ctc/images/pspecs/'+dirname+'/pspec1'+("%03i" % (f+1))+'.fits', clobber=True) 


###################################
###       ADDITIONAL TESTS      ###
###################################
"""
###test what lmax is needed

test_ks = numpy.arange(min(k_data),max(k_data),0.001)
max_Pk = Pk_interp(test_ks[0])

for i in range(len(test_ks)): #find highest k with highest Pk
    test_Pk = Pk_interp(test_ks[i])
    if test_Pk >= max_Pk:
        max_Pk = test_Pk
        max_k = test_ks[i] 
        index = i

print 'k=',max_k,', Pkmax=',max_Pk

factor = 0.1 #10% below max Pk
test2_ks = test_ks[index:]

for i in range(len(test2_ks)): #find k when Pk drops to some factor of its max
    test2_Pk = Pk_interp(test2_ks[i])
    if test2_Pk < factor*max_Pk:
        print 'kmax=',test2_ks[i],', Pk=',test2_Pk
        break
"""    

"""
###test a_lm function (recovering T)

cov_T = [[10,5,2],[5,15,5],[2,5,20]] #symmetric

print cov_T

xxt_test = numpy.zeros_like(cov_T)
num = 1

for i in range(num):

    x = a_lm(cov_T)
    x = numpy.reshape(x,(3,1))
    xt = numpy.conj(numpy.transpose(x))

    xxt = numpy.dot(x,xt)
    xxt_test = xxt_test + xxt

print xxt_test/num #(should be T)
"""

"""
###make single frequency map

plt.plot(k_data,Pk_data,'k.')
plt.show()

#Pk_interp = interp1d(k_data, Pk_data, kind='linear')

l_vals = numpy.arange(0,opts.lmax,1) 

Cl = C_l(.15,.15,Pk_interp,l_vals)

#print Cl

plt.plot(l_vals,Cl, 'k.') #l vs. Cl plot
plt.show() #see where l peaks if P(k) is delta function

a_lms = aipy.healpix.Alm(opts.lmax,opts.lmax) #a_lm object

for i in range(len(Cl)):

    Csmall = [[Cl[i]]] #one block of C 

    l_val = l_vals[i]

    const = 0
    
    for j in range(l_val+1):

        alm = a_lm(Csmall)
        if (l_val-const) == 0:
            alm = numpy.real(alm)
        a_lms[l_val,l_val-const] = complex(alm[0])
        const += 1

sky_map = aipy.map.Map(nside=512)
sky_map.from_alm(a_lms)
sky_map.to_fits('/Users/carinacheng/capo/ctc/images/test10.fits', clobber=True)   
"""

"""
###plotting C_l as a function of delta_nu for a particular l (opts.max)

freq1 = .1
freq2 = numpy.arange(0.1,0.15,0.001) #[GHz]
freq2all = numpy.arange(0.007,0.2,0.001)
freq3 = numpy.arange(0.025,0.1,0.001)
freq4 = numpy.arange(0.03,0.18,0.01)

k_data = numpy.arange(0.0,0.006,0.0001) #actual data points
Pk_data = P_k(k_data)

Pk_interp = interp1d(k_data, Pk_data, kind='linear')

l_vals = [opts.lmax]

Cl_vals = []
deltanu_vals = []
deltaDc_vals = []
Cl_vals_all = []
deltanu_vals_all = []
Cl_vals3 = []
deltanu_vals3 = []
Cl_vals4 = []
deltanu_vals4 = []

#print C_l(freq1,0.05,Pk_interp,l_vals)

for i in range(len(freq2)):

    Cl = C_l(freq1,freq2[i],Pk_interp,l_vals)

    #delta_Dc = Cl[1]-Cl[2]
    #Cl = Cl[0]
    #deltaDc_vals.append(delta_Dc)

    Cl_vals.append(float(Cl[0]))
    deltanu_vals.append(freq1-freq2[i])

    #plt.plot(numpy.array(freq2[i]),numpy.array(Cl[2]),'r.')
    plt.plot(deltanu_vals,Cl_vals,'r.')
    #plt.plot(deltaDc_vals,Cl_vals,'b.')

for i in range(len(freq2all)):

    Clall = C_l(freq1,freq2all[i],Pk_interp,l_vals)
    Cl_vals_all.append(float(Clall[0]))
    deltanu_vals_all.append(freq1-freq2all[i])
    plt.plot(deltanu_vals_all,Cl_vals_all,'k-')

for i in range(len(freq3)):

    Cl3 = C_l(freq1,freq3[i],Pk_interp,l_vals)
    Cl_vals3.append(float(Cl3[0]))
    deltanu_vals3.append(freq1-freq3[i])
    plt.plot(deltanu_vals3,Cl_vals3,'b.')

for i in range(len(freq4)):

    Cl4 = C_l(freq1,freq4[i],Pk_interp,l_vals)
    Cl_vals4.append(float(Cl4[0]))
    deltanu_vals4.append(freq1-freq4[i])
    plt.plot(deltanu_vals4,Cl_vals4,'go')

plt.xlabel('$\Delta$nu')
plt.ylabel('C$_{l}$')
plt.title('l='+str(opts.lmax))
plt.ylim(-2E-12,3.5E-12)
plt.show()
"""



