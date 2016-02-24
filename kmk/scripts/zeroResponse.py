#! /usr/bin/env python
import aipy as a
import capo
import numpy as np
import matplotlib.pylab as pl
import optparse, sys

data_dir = '/home/kara/capo/kmk/data/'
data_file = 'zen.2456895.51490.xx.uvcRRE'
calfile = 'test'

def bl2delay(array, bl):
    i, j = a.miriad.bl2ij(bl)
    bx, by, bz = aa.get_baseline(i, j)
    return np.sqrt(bx**2 + by**2 + bz**2)

def calcFreq(array, bl, ref_bl, ref_freq, min_freq, max_freq):
    f = ref_freq * (ref_bl / baseline)
    if max_freq >= f >= min_freq:
        return f
    else:
        raise ValueError("invalid frequency in baseline %d,%d" % (i,j))

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
    vis = np.sum(np.where(z>0, obs_sky*phs, 0))
    return vis

# select 150 MHz and 160 MHz for u-mode calibration test
freq = 0.100

flatSky = makeFlatMap(nside=64, Tsky=1.0, freq=freq)
xyz = flatSky.px2crd(np.arange(flatSky.npix())) #topocentric

# create array of baselines
test_ants = '(64)_(29,24,28,55,34,27,51,57)'
parsed_ants = a.scripting.parse_ants(test_ants,8)
bl = []
for i in xrange(len(parsed_ants)):
    bl.append(parsed_ants[i][0])
bl = np.array(bl)

aa = a.cal.get_aa(calfile, np.array([freq]))
ref_bl = bl2delay(aa, bl.max())
freqs = []
for i in xrange(len(bl)):
    baseline = bl2delay(array=aa, bl=bl[i])
    freqs.append(calcFreq(array=aa, bl=baseline, ref_bl=ref_bl, ref_freq=freq, min_freq=0.100, max_freq=0.200))
freqs = np.array(freqs)

# import array parameters
aa = a.cal.get_aa(calfile, freqs)
aa_freqs = aa.get_freqs()
beam = aa[0].bm_response(xyz, pol = 'x')**2
print "array parameters imported"

for i in xrange(len(bl)):
    #flatMap = makeFlatMap(nside=64, Tsky=1.0, freq=test_freqs[i])
    #gsmMap = makeGSMMap(gsm_dir=gsm_dir, filename=gsm_files[i], nside=64, freq=test_freqs[i], array=aa)
    gsm_vis = calcVis(hpm=flatSky, beam=beam[i], bl=bl[i], coord=xyz, freq=freqs[i])
    print gsm_vis
