#!/usr/bin/env python

"""

NAME: 
      pspec.py 
PURPOSE:
      -Generates a cube in k-space filled with T_tilde values that are generated from a Gaussian distribution with variance specified by P(k)
      -Fourier transforms the cube to get T(r)
EXAMPLE CALL:
      ./pspec.py 
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
from collections import defaultdict #for multidimensional dictionary


#make and fill cube

k0 = 0
sigma = 1 #adjustable parameters

k_min = 0
k_max = 2
k_step = 1 #adjustable parameters

kx = numpy.arange(-k_max+1,k_max,k_step)
ky = numpy.arange(-k_max+1,k_max,k_step)
kz = numpy.arange(k_min,k_max,k_step) #only fill half the values

k_map = []

for i in range(len(kx)):
    for j in range(len(ky)):
        for m in range(len(kz)):

            kname = str(kx[i])+str(ky[j])+str(kz[m]) #string (ex: '-2-20')
            kmag = numpy.sqrt(kx[i]**2+ky[j]**2+kz[m]**2)
            P_k = numpy.exp(-(kmag-k0)**2/(2*sigma**2)) #variance
            stdev = numpy.sqrt(P_k/2) #stdev

            if kx[i]==0 and ky[j]==0 and kz[m]==0: #no complex part
                
                stdev = numpy.sqrt(P_k)
                T_tilde = numpy.random.randn()*stdev+k0+0j
                k_map.append((kname, T_tilde))

            else:

                a_tilde = numpy.random.randn()*stdev+k0 
                b_tilde = numpy.random.randn()*stdev+k0 #random num with variance P_k/2
                #sampled from Gaussian distribution of mean k0 and variance 'stdev**2'
                T_tilde = a_tilde+1j*b_tilde
                k_map.append((kname, T_tilde))

d = defaultdict(list)
for index, value in k_map:
    d[index].append(value)

kz_fill = numpy.arange(-k_max+1,k_min,k_step)

for i in range(len(kx)):
    for j in range(len(ky)):
        for m in range(len(kz_fill)):
            kname = str(-1*kx[i])+str(-1*ky[j])+str(kz_fill[m]) 
            kname_pos = str(kx[i])+str(ky[j])+str(-1*kz_fill[m])
            T_tilde = numpy.conj(d[kname_pos])[0] #conjugate value
            k_map.append((kname, T_tilde))

d = defaultdict(list)
for index, value in k_map:
    d[index].append(value)

###can look at T_tildes here by calling d['kxkykz']

#putting knames and T_tildes in array

names = []
T_tildes = []
for i in range(len(k_map)):
    names.append(k_map[i][0])
    T_tildes.append(k_map[i][1])

#reordering knames and T_tildes to place in 3-dim array

cube_size = len(kx)

names_arr = numpy.zeros((cube_size,cube_size,cube_size),dtype=(str,6))
tildes_arr = numpy.zeros((cube_size,cube_size,cube_size),dtype=complex)
index1 = 0
index2 = 0
index3 = 0
count = 0
for i in range(cube_size):
    index2 = 0
    for j in range(cube_size):
        index3 = 0
        for k in range(cube_size):
            names_arr[index1,index2,index3] = str(names[count])
            tildes_arr[index1,index2,index3] = T_tildes[count]
            index3 += 1
            count += 1
        index2 +=1
    index1 +=1

#fourier transform cube

print names_arr
print tildes_arr

T_r = numpy.fft.ifftn(tildes_arr,axes=(0,1,2))

#print T_r

knames = numpy.ndarray.flatten(names_arr)
T_r = numpy.ndarray.flatten(T_r)

T_map = []

for i in range(len(knames)):
    T_map.append((knames[i], T_r[i]))
d = defaultdict(list)
for index, value in T_map:
    d[index].append(value)

###can look at T_r here by calling d['kxkykz']






###positive k values only

"""
#make and fill cube

k0 = 1
sigma = 1 #adjustable parameters

k_min = 0
k_max = 4
k_step = 1 #adjustable parameters

kx = numpy.arange(k_min,k_max,k_step)
ky = numpy.arange(k_min,k_max,k_step)
kz = numpy.arange(k_min,(k_max-k_min)/2,k_step) #only fill half the values
print kz

k_map = []

for i in range(len(kx)):
    for j in range(len(ky)):
        for m in range(len(kz)):
            kname = str(kx[i])+str(ky[j])+str(kz[m]) #string (ex: '-2-20')
            kmag = numpy.sqrt(kx[i]**2+ky[j]**2+kz[m]**2)
            P_k = numpy.exp(-(kmag-k0)**2/(2*sigma**2)) #variance
            stdev = numpy.sqrt(P_k/2) #stdev
            a_tilde = numpy.random.randn()*stdev+k0 
            b_tilde = numpy.random.randn()*stdev+k0 #random num with variance P_k/2
            #sampled from Gaussian distribution of mean k0 and variance 'stdev**2'
            T_tilde = a_tilde+1j*b_tilde
            k_map.append((kname, T_tilde))

d = defaultdict(list)
for index, value in k_map:
    d[index].append(value)

kz_fill = numpy.arange((k_max-k_min)/2,k_max,k_step) #other half of kz values

for i in range(len(kx)):
    for j in range(len(ky)):
        for m in range(len(kz_fill)):
            kname = str(kx[i])+str(ky[j])+str(kz_fill[m]) 
            maxk = numpy.max(kx)
            kname_pos = str(maxk-kx[i])+str(maxk-ky[j])+str(maxk-kz_fill[m])
            T_tilde = numpy.conj(d[kname_pos])[0] #conjugate value
            k_map.append((kname, T_tilde))

d = defaultdict(list)
for index, value in k_map:
    d[index].append(value)

###can look at T_tildes here by calling d['kxkykz']

#putting knames and T_tildes in array

names = []
T_tildes = []
for i in range(len(k_map)):
    names.append(k_map[i][0])
    T_tildes.append(k_map[i][1])

#reordering knames and T_tildes to place in 3-dim array

cube_size = len(kx)
names_reordered = []
T_tildes_reordered = []
start = 0
start2 = start+len(names)/2
for i in range((len(names)/2)/(cube_size/2)):
    for j in range(cube_size/2):
        names_reordered.append(names[start])
        T_tildes_reordered.append(T_tildes[start])
        start += 1
    for j in range(cube_size/2):
        names_reordered.append(names[start2])
        T_tildes_reordered.append(T_tildes[start2])
        start2 += 1

names_arr = numpy.zeros((cube_size,cube_size,cube_size),dtype=(str,3))
tildes_arr = numpy.zeros((cube_size,cube_size,cube_size),dtype=complex)
index1 = 0
index2 = 0
index3 = 0
count = 0
for i in range(cube_size):
    index2 = 0
    for j in range(cube_size):
        index3 = 0
        for k in range(cube_size):
            names_arr[index1,index2,index3] = str(names_reordered[count])
            tildes_arr[index1,index2,index3] = T_tildes_reordered[count]
            index3 += 1
            count += 1
        index2 +=1
    index1 +=1

#fourier transform cube

print tildes_arr

T_r = numpy.fft.ifftn(tildes_arr,axes=(0,1,2))

print T_r

knames = numpy.ndarray.flatten(names_arr)
T_r = numpy.ndarray.flatten(T_r)

T_map = []

for i in range(len(knames)):
    T_map.append((knames[i], T_r[i]))
d = defaultdict(list)
for index, value in T_map:
    d[index].append(value)

###can look at T_r here by calling d['kxkykz']
"""
