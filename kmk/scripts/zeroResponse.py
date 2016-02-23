#! /usr/bin/env python
import aipy as a
import capo
import numpy as np
import matplotlib.pylab as pl
import optparse, sys

data_dir = '/home/kara/capo/kmk/data/'
data_file = 'zen.2456895.51490.xx.uvcRRE'
calfile = 'test'

def makeFlatMap(nside, freq, Tsky=1.0):
    # NOTE: optimal baseline length formula given in Presley et al 2015 eq. 9
    # fill sky map with flat temperature across whole sky (DC signal)
    hpm = a.healpix.HealpixMap(nside=nside)
    hpm.map = Tsky*np.ones(shape = hpm.map.shape)
    print "flat map size = " + str(hpm.map.shape)
    return hpm_TtoJy(hpm, freq)

def hpm_TtoJy(hpm, freq):
    wvlen = a.const.c / freq
    hpm.map *= (4*np.pi/hpm.npix())*2*a.const.k/wvlen**2 # convert to Jy to make summable
    return hpm

# fill sky map with GSM
# create array of GSM files for simulation of multiple frequencies
# MODIFY TO INCLUDE MULTIPLE GSM FILES AT DIFFERENT FREQUENCIES
gsm_dir = '/home/kara/capo/kmk/gsm/gsm_raw/gsm'
gsm_files = np.array([1001, 1002])

# NOTE: the gsm file is associated with a particular frequency map --
# check to ensure filename matches input frequency
def makeGSMMap(array, nside, filename, freq, gsm_dir='/home/kara/capo/kmk/gsm/gsm_raw/'):
    g = a.healpix.HealpixMap(nside=512)
    gsm = a.healpix.HealpixMap(nside=nside)
    print "GSM map size = " + str(gsm.map.shape)
    d = np.loadtxt(gsm_dir + str(filename) + '.dat')
    g.map = d
    g = hpm_TtoJy(g, freq)
    gsm.from_hpm(g) # hack to lower resolution to prevent memory overload
    # convert to topocentric coordinates
    ga2eq = a.coord.convert_m('eq', 'ga', oepoch=array.epoch) #conversion matrix
    eq2top = a.coord.eq2top_m(array.sidereal_time(),array.lat) #conversion matrix
    ga2eq2top = np.dot(eq2top, ga2eq)
    i, j, k = ijk = np.dot(ga2eq2top,gsm.px2crd(np.arange(gsm.npix()))) #topocentric
    return gsm

# create arrays for visibilities for each baseline
# single timestep, single polarization
# final index is frequency (to be implemented in the future)

# define our observation equation in the Tegmark 1997 form
# Y = AX + N
# Y = observed visibility
# A = theoretical visibility
# X = constant temperature value for global signal
# N = noise

def extractData(uv):
    time = []
    f = {} # dictionary of flags
    d = {} # dictionary of spectra information
    for (uvw, t, (m, n)), data, flag in uv.all(raw='True'):
        # get visibility information for rest frequency
        bl = a.miriad.ij2bl(m, n)
        if not bl in d.keys(): f[bl] = []; d[bl] = []
        f[bl].append(flag)
        d[bl].append(data)
        time.append(t)
    print "data collected"
    return d

# simulate visibilities for each baseline
def calcVis(hpm, beam, bl, coord, freq):
    x,y,z = coord
    i, j = a.miriad.bl2ij(bl)
    print i, j
    bx, by, bz = aa.get_baseline(i, j)
    # attenuate sky signal and visibility by primary beam
    obs_sky = hpm.map
    #obs_sky = beam * hpm.map
    phs = np.exp(np.complex128(-2j*np.pi*freq*(bx*x + by*y + bz*z)))
    print 'bx = ' + str(bx)
    print 'phase = ' + str(phs)
    vis = np.sum(np.where(z>0, obs_sky*phs, 0))
    #vis = np.sum(np.where(z>0, obs_sky, 0))
    return vis

# select 150 MHz and 160 MHz for u-mode calibration test
freq = 0.100
#freqs = freq*np.array([2.0, 1.0])
freqs = freq*np.array([15.0/8, 15.0/9, 15.0/10, 15.0/11, 15.0/12, 15.0/13, 15.0/14, 1.0])
#test_freqs = 1e9*np.array([aa_freqs[102], aa_freqs[122]])

flatSky = makeFlatMap(nside=64, Tsky=1.0, freq=freq)
xyz = flatSky.px2crd(np.arange(flatSky.npix())) #topocentric

# import array parameters
aa = a.cal.get_aa(calfile, freqs)
# select freq = 150 MHz
aa_freqs = aa.get_freqs()
beam = aa[0].bm_response(xyz, pol = 'x')**2
print "array parameters imported"

#test_ants = '(64)_(10,49)'
test_ants = '(64)_(29,24,28,55,34,27,51,57)'
parsed_ants = a.scripting.parse_ants(test_ants,len(freqs))
bl = []
for i in xrange(len(parsed_ants)):
    bl.append(parsed_ants[i][0])

for i in xrange(len(bl)):
    #flatMap = makeFlatMap(nside=64, Tsky=1.0, freq=test_freqs[i])
    #gsmMap = makeGSMMap(gsm_dir=gsm_dir, filename=gsm_files[i], nside=64, freq=test_freqs[i], array=aa)
    gsm_vis = calcVis(hpm=flatSky, beam=beam[i], bl=bl[i], coord=xyz, freq=freqs[i])
    print gsm_vis
