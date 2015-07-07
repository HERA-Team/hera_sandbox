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
temps = 0.5*np.ones(1)
x, y, z = xyz = h.px2crd(np.arange(h.npix())) #topocentric

freq = uv.__getitem__('restfreq')*np.ones(1)
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

ants = '64_10,65_9,72_22,80_20,88_43,96_53,104_31'

time = []
f = {} # dictionary of flags
d = {} # dictionary of spectra information
a.scripting.uv_selector(uv, ants, 'xx')
for (uvw, t, (i,j)), data, flag in uv.all(raw='True'):
    # get visibility information for rest frequency
    bl = a.miriad.ij2bl(i,j)
    if not bl in d.keys(): f[bl] = []; d[bl] = []
    f[bl].append(flag)
    d[bl].append(data)
    time.append(t)

# import array parameters
aa = a.cal.get_aa('psa6622_v001', uv['sdf'], uv['sfreq'], uv['nchan'])
beam = aa[0].bm_response(xyz, pol='x')**2
beam = beam[0]

# simulate visibilities for each baseline
response = {}
for bl in d.keys():
    k, j = a.miriad.bl2ij(bl)
    bx, by, bz = aa.get_baseline(k, j)
    for i in range(len(freq)):
        # attenuate sky signal by primary beam
        obs = beam * h.map
        phs = np.exp(-2j*np.pi*freq[i]*(bx*x + by*y + bz*z))
        if not bl in response.keys(): response[bl] = []
        response[bl].append(np.sum(np.where(z>0, obs*phs, 0)))
        i += 1

A = np.array([response[bl][0] for bl in response.keys()])
A.shape = (A.size,1)
Y = A*temps[0] + np.random.normal(size=A.shape) + 1j*np.random.normal(size=A.shape)
#Y = np.array([d[bl][0][N/2-1] for bl in d.keys()])

# find value of X from cleverly factoring out A in a way which properly weights 
# the measurements
transjugateA = np.conjugate(np.transpose(A))
normalization = np.linalg.inv(np.dot(transjugateA, A))
invA = normalization*transjugateA

X = np.dot(invA,Y)

# print the estimate of the global signal
print X
