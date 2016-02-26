#! /usr/bin/env python
import aipy as a
import capo
import numpy as np
import matplotlib.pylab as pl
import optparse, sys
import os

sim_dir = '/home/kara/capo/kmk/scripts/'
data_dir = '/home/kara/capo/kmk/data/'
data_file = 'zen.2456895.51490.xx.uvcRRE'
calfile = 'test'
gsm_dir = '/home/kara/capo/kmk/gsm/gsm_raw/'

def bl2delay(array, bl):
    """ converts unique baseline numbers into nanosecond delays """
    i, j = a.miriad.bl2ij(bl)
    bx, by, bz = aa.get_baseline(i, j)
    return np.sqrt(bx**2 + by**2 + bz**2)

def calcFreq(array, bl, ref_bl, ref_freq, min_freq, max_freq):
    """ calculates nuCal frequencies based on an initial frequency, reference baseline, 
        and input baseline. must include min and max frequencies for array in order to
        ensure that calculated frequencies are legal. """
    f = ref_freq * (ref_bl / baseline)
    if max_freq >= f >= min_freq:
        return f
    else:
        raise ValueError("invalid frequency in baseline %d,%d" % (i,j))

def makeFlatMap(nside, freq, Tsky=1.0):
    """ fill sky map of given size and frequency with flat temperature across whole sky,
        returns map in Janskys. """
    hpm = a.healpix.HealpixMap(nside=nside)
    hpm.map = Tsky*np.ones(shape = hpm.map.shape)
    print "flat map size = " + str(hpm.map.shape)
    return hpm_TtoJy(hpm, freq)

def hpm_TtoJy(hpm, freq):
    """ converts given Healpix map from brightness temperature to Janskys, provided map
        frequency """
    wvlen = a.const.c / freq
    hpm.map *= (4*np.pi/hpm.npix())*2*a.const.k/wvlen**2 # convert to Jy to make summable
    return hpm

def makeGSM(path, filename, sfreq, sdf, num):
    """ runs the MIT multi-frequency global sky model simulation at given frequencies 

        required input variables:
        path: provides pathname to gsmmf.sh -- need to be in this directory to run sim
        filename: desired output filenames for GSM maps
        sfreq: desired starting frequency for your GSM maps
        sdf: desired frequency spacings of maps
        num: desired number of maps produced 
        
        GSM simulator can be downloaded at 
        <<http://space.mit.edu/~angelica/gsm/index.html>> """
    os.chdir(path)
    os.system('rm -rf args.dat')
    f = open('args.dat', 'w')
    s = "%s %f %f %d" % (filename, sfreq, sdf, num)
    print s
    f.write(s)
    f.close()
    os.system(path+'gsmmf.sh')

def makeGSMMap(array, nside, filename, freq, path='/home/kara/capo/kmk/gsm/gsm_raw/'):
    """ create a Healpix map of a given size filled with a simulated global sky model
        at a given frequency """
    makeGSM(path=path, filename=filename, sfreq=freq, sdf=10, num=1)
    g = a.healpix.HealpixMap(nside=512)
    gsm = a.healpix.HealpixMap(nside=nside)
    print "GSM map size = " + str(gsm.map.shape)
    d = np.loadtxt(path + filename + str(1001) + '.dat')
    g.map = d
    g = hpm_TtoJy(g, freq)
    gsm.from_hpm(g) # hack to lower resolution to prevent memory overload
    # convert to topocentric coordinates
    ga2eq = a.coord.convert_m('eq', 'ga', oepoch=array.epoch) #conversion matrix
    eq2top = a.coord.eq2top_m(array.sidereal_time(),array.lat) #conversion matrix
    ga2eq2top = np.dot(eq2top, ga2eq)
    i, j, k = ijk = np.dot(ga2eq2top,gsm.px2crd(np.arange(gsm.npix()))) #topocentric
    return gsm


def extractData(uv):
    """ create arrays for visibilities for each baseline in uv data file for 
        single timestep, single polarization

        define our observation equation in the Tegmark 1997 form
        Y = AX + N
        Y = observed visibility
        A = theoretical visibility
        X = constant temperature value for global signal
        N = noise 

        NOTE: final index in data directory is frequency 
        (to be implemented in the future) """
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

def calcVis(hpm, beam, bl, coord, freq):
    """ simulate sky visibilities for a given baseline and primary beam, 
        provided with a sky map at a known frequency and its coordinate system """
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

sim_data = []
for i in xrange(len(bl)):
    gsmMap = makeGSMMap(path=gsm_dir, filename='gsm', nside=64, freq=freqs[i], array=aa)
    gsm_vis = calcVis(hpm=gsmMap, beam=beam[i], bl=bl[i], coord=xyz, freq=freqs[i])
    vis_data = [a.miriad.bl2ij(bl[i]), freqs[i], gsm_vis]
    sim_data.append(vis_data)

np.array(sim_data)
np.savez(sim_dir+'sim_output',sim_data)
print sim_data
