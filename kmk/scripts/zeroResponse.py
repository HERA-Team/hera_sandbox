#! /usr/bin/env python
import aipy as a
import capo
import numpy as np
import matplotlib.pylab as pl

data_dir = '/home/kara/capo/kmk/data/'
data_file = 'zen.2456895.51490.xx.uvcRRE'
uv = a.miriad.UV(data_dir + data_file)
N = uv.__getitem__('nchan')

# NOTE: optimal baseline length formula given in Presley et al 2015 eq. 9
# fill h with flat temperature across whole sky (DC signal)
h = a.healpix.HealpixMap(nside=64)
h.map = np.ones(shape = h.map.shape)
h.map = h.map*4*np.pi/h.npix()
temps = np.ones(1)
#temps = np.array([0.5, 0.6, 0.8]) # sky temperature
x, y, z = xyz = h.px2crd(np.arange(h.npix())) #topocentric

freq = uv.__getitem__('restfreq')*np.ones(1)
#freq = np.array([.120, .150, .180]) # in GHz
c = 3e10 # in cm/s
wvlen = c / freq # observed wavelength

# create arrays for visibilities for each baseline
# single timestep, single polarization
# final index is frequency (to be implemented in the future)

# define our observation equation in the Tegmark 1997 form
# Y = AX + N
# Y = observed visibility
# A = theoretical visibility
# X = constant temperature value for global signal
# N = noise
Y = np.zeros((128, 128, len(freq)), dtype=np.complex)
A = np.zeros(shape = Y.shape, dtype=np.complex) 

time = []
f = {} # dictionary of flags
d = {} # dictionary of spectra information
for (uvw, t, (i,j)), data, flag in uv.all(raw='True'):
    # get visibility information for rest frequency
    bl = a.miriad.ij2bl(i,j)
    if not bl in d.keys(): f[bl] = []; d[bl] = []
    f[bl].append(flag)
    d[bl].append(data)
    time.append(t)

for bl in d.keys():
    i,j = a.miriad.bl2ij(bl)
    # use only first time step for simplicity
    Y[i,j,0] = d[bl][-1][N/2-1]

# import array parameters
aa = a.cal.get_aa('psa6622_v001', uv['sdf'], uv['sfreq'], uv['nchan'])
beam = aa[0].bm_response(xyz, pol='x')**2
beam = beam[0]

# simulate visibilities for each baseline
for bl in d.keys():
    k,j = a.miriad.bl2ij(bl)
    bx, by, bz = aa.get_baseline(k, j)
    for i in range(len(temps)):
        # attenuate sky signal by primary beam
        obs = temps[i] * beam * h.map
        phs = np.exp(-2j*np.pi*freq[i]*(bx*x + by*y + bz*z))
        A[k,j,i] = np.sum(np.where(z>0,obs*phs,0))
        if j != k: A[j,k,i] = np.sum(np.where(z>0,obs*phs,0))
        i += 1

# find value of X from cleverly factoring out A
transjugateA = np.conjugate(np.transpose(A))
invA = np.linalg.inv(transjugateA*A)*transjugateA

X = invA*Y

# under the assumption that noise is Gaussian distributed around zero
# sum together all the estimates of X in order to average out the noise
print np.sum(X)

counts = 0
mean = 0

for i in np.arange(128):
    for j in np.arange(128):
        if X[i,j,0] != 0:
            counts += 1
            mean += X[i,j,0]

print mean/counts
