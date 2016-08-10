#!/usr/bin/env python

"""

NAME: 
      pspec_sim.py 
PURPOSE:
      -Generates a cube in k-space filled with T_tilde values that are generated from a Gaussian distribution with variance specified by P(k)
      -Fourier transforms the cube to get T(r)
      -Plots a slice of the cube onto a Healpix Map
EXAMPLE CALL:
      ./pspec_sim.py -N 31 --freq 150
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
from scipy import integrate
from collections import defaultdict #for multidimensional dictionary

#options

o = optparse.OptionParser()
o.set_description(__doc__)
o.add_option('-N', dest='N', default = 3, type = 'int',
             help='Number of pixels per side of the T_tilde cube. Must be an odd number. Default is 3.')
o.add_option('--freq', dest='freq', default = 150., type = 'float',
             help='Observed frequency in MHz. Default is 150.')
opts, args = o.parse_args(sys.argv[1:])

#P(k) function

def P_k(kmag, sigma=0.01, k0=0.02):

    return numpy.exp(-(kmag-k0)**2/(2*sigma**2)) #variance

#cosmology calculations

nu21 = 1420.*10**6 #21cm frequency [Hz]
nuobs = opts.freq*10**6 #observed frequency [Hz]
c = 3.*10**5 #[km/s]
H0 = 100 #Hubble's constant [km/s/Mpc]
omg_m = 0.28 
omg_lambda = 0.72 

z = (nu21/nuobs)-1 #redshift

one_over_E_z = lambda zprime:1/numpy.sqrt(omg_m*(1+zprime)**3+omg_lambda)
Dc = (c/H0)*integrate.quad(one_over_E_z,0,z)[0] #comoving distance [Mpc]

print 'redshift z =', z
print 'co-moving distance for', opts.freq, 'MHz [Mpc] =', Dc

#set resolution

wavelength = (3.*10**8)/nuobs #[m]
baseline = 30. #[m]
theta = wavelength/baseline

#delta = theta*Dc #step size in real space [Mpc]
delta = 45 #test case
L = opts.N*delta #size range in real space [Mpc]

print 'real space volume [Mpc] =', L,'x',L,'x',L
print 'real space resolution [Mpc] =', delta
print 'angular scale range on sky [rad] =', L/Dc
print 'angular resolution on sky [rad] =', theta

#make and fill cube

kx = numpy.fft.fftfreq(int(opts.N),delta)*2*numpy.pi #the 2pi corrects for fourier conventions
kx = kx.copy(); kx.shape = (kx.size,1,1)
ky = kx.copy(); ky.shape = (1,ky.size,1)
kz = kx.copy(); kz.shape = (1,1,kz.size)

#print kx[1,0,0],ky[0,1,0],kz[0,0,1]

k_mag = numpy.sqrt(kx**2+ky**2+kz**2) #3D cube
#print numpy.min(k_mag), numpy.max(k_mag)

stdev = numpy.sqrt(L**3*P_k(k_mag)/2)
a_tilde = numpy.random.normal(scale=stdev) 
b_tilde = numpy.random.normal(scale=stdev) #random num with variance P_k/2
#sampled from Gaussian distribution of variance 'stdev**2'
T_tilde = a_tilde+1j*b_tilde

#make cube Hermitian

thing1 = T_tilde[1:1+(kx.size-1)/2,1:,1:] #eliminates 0th rows and cuts cube in half
T_tilde[:(kx.size-1)/2:-1,:0:-1,:0:-1] = numpy.conj(thing1)

thing2 = T_tilde[1:1+(kx.size-1)/2,:1,1:] 
T_tilde[:(kx.size-1)/2:-1,0::-1,:0:-1] = numpy.conj(thing2)

thing3 = T_tilde[:1,1:,1:1+(kz.size-1)/2]
T_tilde[:1,:0:-1,:(kz.size-1)/2:-1] = numpy.conj(thing3)

thing4 = T_tilde[1:1+(kx.size-1)/2,1:,:1]
T_tilde[:(kx.size-1)/2:-1,:0:-1,:1] = numpy.conj(thing4)

thing5 = T_tilde[1:1+(kx.size-1)/2,:1,:1]
T_tilde[:(kx.size-1)/2:-1,:1,:1] = numpy.conj(thing5)

thing6 = T_tilde[:1,1:1+(ky.size-1)/2,:1]
T_tilde[:1,:(ky.size-1)/2:-1,:1] = numpy.conj(thing6)

thing7 = T_tilde[:1,:1,1:1+(kz.size-1)/2]
T_tilde[:1,:1,:(kz.size-1)/2:-1] = numpy.conj(thing7)

T_tilde[0,0,0] = numpy.real(T_tilde[0,0,0])

T_r = numpy.fft.ifftn(T_tilde)/(delta**3)

#plot slice of cube

slice = T_r[opts.N/2]
plt.imshow(numpy.real(slice))
plt.show()
otherslice = []
for i in T_r:
    otherslice.append(i[opts.N/2])
plt.imshow(numpy.real(otherslice))
plt.show()
otherslice2 = numpy.zeros_like(otherslice)
for i in range(len(T_r)):
    slice = T_r[i]
    for j in range(len(slice)):
        otherslice2[i][j] = slice[j][opts.N/2]
plt.imshow(numpy.real(otherslice2))
plt.show()


"""
#print T_r
#numpy.save('Tr_test.npy',T_r)

#make Healpix Map

img = aipy.map.Map(nside=128) #needs to be nside=512
px = numpy.arange(img.npix())
thetas, phis = img.px2crd(px,ncrd=2)

center_index = (opts.N)/2.-0.5

for i in range(len(thetas)):

    delta_x = Dc*numpy.sin(thetas[i])*numpy.cos(phis[i]) #spherical coordinates
    delta_y = Dc*numpy.sin(thetas[i])*numpy.sin(phis[i])
    delta_z = Dc*numpy.cos(thetas[i])

    px_x = delta_x/delta #number of pixels to move by from the center
    px_y = delta_y/delta
    px_z = delta_z/delta

    if (px_x < opts.N/2) and (px_y < opts.N/2) and (px_z < opts.N/2) and (px_x > -opts.N/2) and (px_y > -opts.N/2) and (px_z > -opts.N/2):

        value = T_r[center_index+round(px_x),center_index+round(px_y),center_index+round(px_z)]
        img.put((thetas[i],phis[i]),1.0,value.real)

img.to_fits('/Users/carinacheng/capo/ctc/code/testcube1.fits', clobber=True)
"""






"""
##OLD VERSION OF CODE: succesfully makes Hermitian cubes but is inefficient

k_min = 0
k_max = 3
k_step = 1

kx = numpy.arange(k_min,k_max,k_step) #size of these must be odd #
ky = numpy.arange(k_min,k_max,k_step)
kz = numpy.arange(k_min,k_max,k_step) 

k_map = []

#strformat = "%05.1f"
strformat = "%1.0f" #test case (1 digit)

#0th frequency

indices = []
indices.append('000')
kname = str(strformat % k_min)+str(strformat % k_min)+str(strformat % k_min)
kmag = numpy.sqrt(k_min**2+k_min**2+k_min**2)
P_k = numpy.exp(-(kmag-k0)**2/(2*sigma**2))
stdev = numpy.sqrt(P_k) 
T_tilde = numpy.random.randn()*stdev+k0+0j #no complex part
k_map.append((kname,T_tilde))

#other freqs 

k_names_pos = []
k_names_neg = []
T_tildes_pos = []
T_tildes_neg = []
indices_pos = []
indices_neg = []

for i in range(len(kx)):
    for j in range(len(ky)):
        for m in range(len(kz)):

            kname = str(strformat % kx[i])+str(strformat % ky[j])+str(strformat % kz[m]) #string (ex: '220')

            if kx[i]+ky[j]+kz[m] != 0 and (kname not in k_names_pos) and (kname not in k_names_neg):

                k_names_pos.append(kname)
                kmag = numpy.sqrt(kx[i]**2+ky[j]**2+kz[m]**2)
                P_k = numpy.exp(-(kmag-k0)**2/(2*sigma**2)) #variance
                stdev = numpy.sqrt(P_k/2) #stdev
                a_tilde = numpy.random.randn()*stdev+k0 
                b_tilde = numpy.random.randn()*stdev+k0 #random num with variance P_k/2
                #sampled from Gaussian distribution of mean k0 and variance 'stdev**2'
                T_tilde = a_tilde+1j*b_tilde
                T_tildes_pos.append(T_tilde)
                indices_pos.append(str(i)+str(j)+str(m))

                max = numpy.max(kx)
        
                if kx[i] == 0:
                    term1 = 0
                else:
                    term1 = max-kx[i]+k_step
                if ky[j] == 0:
                    term2 = 0
                else:
                    term2 = max-ky[j]+k_step
                if kz[m] == 0:
                    term3 = 0
                else:
                    term3 = max-kz[m]+k_step

                maxi = len(kx)-1

                if i == 0:
                    term1i = 0
                else:
                    term1i = maxi-i+1
                if j == 0:
                    term2i = 0
                else:
                    term2i = maxi-j+1
                if m == 0:
                    term3i = 0
                else:
                    term3i = maxi-m+1
            
                kname_neg = str(strformat % term1)+str(strformat % term2)+str(strformat % term3)
                k_names_neg.append(kname_neg)
                T_tilde = numpy.conj(T_tilde)
                T_tildes_neg.append(T_tilde)
                indices_neg.append(str(term1i)+str(term2i)+str(term3i))

for i in range(len(indices_pos)):
    indices.append(indices_pos[i])

for i in range(len(indices_neg)):
    indices.append(indices_neg[::-1][i])

for i in range(len(k_names_pos)):
    k_map.append((k_names_pos[i],T_tildes_pos[i]))

for i in range(len(k_names_neg)):
    k_map.append((k_names_neg[::-1][i],T_tildes_neg[::-1][i]))

d = defaultdict(list)
for index, value in k_map:
    d[index].append(value)

###can look at T_tildes here by k-value (ex: d['kxkykz'])

#putting knames and T_tildes in 1D array (0th freq, then pos freqs, then neg freqs)

k_names = []
T_tildes = []
for i in range(len(k_map)):
    k_names.append(k_map[i][0])
    T_tildes.append(k_map[i][1])

#print k_names, T_tildes
#print k_map

#fourier transform 

arr3d = numpy.zeros((len(kx),len(kx),len(kx)),dtype=complex)

for i in range(len(indices)):
    arr3d[indices[i][0],indices[i][1],indices[i][2]] = T_tildes[i]

T_r = numpy.fft.ifftn(arr3d)

#print T_r

pylab.imshow(numpy.real(T_r[0]),interpolation='nearest')
pylab.colorbar(shrink=0.5)
pylab.show()


###can look at T_r here by index number (ex: T_r[0,0,0])

#linking k's and T_r's

T_rs = []

for i in range(len(indices)):  
    T_rs.append(T_r[indices[i][0],indices[i][1],indices[i][2]])

T_map = []

for i in range(len(k_names)):
    T_map.append((k_names[i], T_rs[i]))
d_r = defaultdict(list)
for index, value in T_map:
    d_r[index].append(value)

###can look at T_r here by k-value (ex: d_r['kxkykz'])

#pylab.plot(numpy.real(T_rs))
#pylab.show()

"""





