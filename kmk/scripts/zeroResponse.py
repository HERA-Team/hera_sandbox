#! /usr/bin/env python
import aipy as a
import capo
import numpy as np
import os
from GlobalSkyModel import GlobalSkyModel

def bl2delay(array, ij):
    """ converts unique baseline numbers into nanosecond delays """
    i,j = ij
    bl = array.get_baseline(i, j)
    return np.sqrt(np.dot(bl, bl))

def calcFreq(array, ij, ref_ij, ref_freq, min_freq, max_freq):
    """ calculates nuCal frequencies (in GHz) based on an initial frequency, reference baseline, 
        and input baseline. must include min and max frequencies for array in order to
        ensure that calculated frequencies are legal. """
    bl = bl2delay(array, ij)
    ref_bl = bl2delay(array, ref_ij)
    f = ref_freq * (ref_bl / bl)
    assert(max_freq >= f >= min_freq)
    return f

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

def makeGSM(path, nside, freq):
        gsmMap = a.healpix.HealpixMap(nside=nside)
        g = a.healpix.HealpixMap(nside=512)
        g.map = GlobalSkyModel(freq=1000*freq, GSMlocation=path, GSMNSIDE=512).map # galactic
        g = hpm_TtoJy(g, freq=freq)
        gsmMap.from_hpm(g)
        return gsmMap

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

def calcVis(aa, hpm, ij, freq, ha):
    # TODO: import array + GSMMap, calculate topocentric coordinates on the 
    # fly, generate PB on the fly, include time
    """ simulate sky visibilities for a given baseline and primary beam, 
        provided with a sky map at a known frequency and its coordinate system """
    gxyz = gx,gy,gz, = hpm.px2crd(np.arange(hpm.npix())) # galactic
    ga2eq = a.coord.convert_m(isys='eq', osys='ga', oepoch=aa.epoch) 
    # conversion matrix so the isys and osys are reversed
    exyz = ex,ey,ez = np.dot(ga2eq, hpm.px2crd(np.arange(hpm.npix()))) # equatorial
    txyz = tx,ty,tz = np.dot(a.coord.eq2top_m(ha, aa.lat), exyz) # topocentric
    # generate proper PB
    beam = aa[0].bm_response(txyz, pol='x')**2 # topocentric
    i,j = ij
    print i, j
    bxyz = bx,by,bz = aa.get_baseline(i, j, src='z')
    # attenuate sky signal and visibility by primary beam
    obs_sky = beam * hpm.map
    phs = np.exp(np.complex128(-2j*np.pi*freq*np.dot(bxyz, txyz)))
    vis = np.sum(np.where(tz>0, obs_sky*phs, 0))
    return vis

if __name__ == '__main__':
    import matplotlib.pylab as pl
    import optparse, sys

    o = optparse.OptionParser()
    o.add_option('--sim_dir', default='/home/kara/capo/kmk/scripts/')
    o.add_option('--data_dir', default='/home/kara/capo/kmk/data/')
    o.add_option('--data_file', default='zen.2456895.51490.xx.uvcRRE')
    o.add_option('--calfile', default='test')
    o.add_option('--gsm_dir', default='/home/kara/capo/kmk/gsm/gsm_raw/')
    o.add_option('--fileout', default='sim_results.uv')

    opts,args = o.parse_args(sys.argv[1:])
    sim_dir = opts.sim_dir
    data_dir = opts.data_dir
    data_file = opts.data_file
    calfile = opts.calfile
    gsm_dir = opts.gsm_dir
    fileout = opts.fileout

    # select 150 MHz and 160 MHz for u-mode calibration test
    freq = 0.100
    aa = a.cal.get_aa(calfile, np.array([freq]))

    flatSky = makeFlatMap(nside=64, Tsky=1.0, freq=freq)

    # create array of baselines
    test_ants = '(64)_(29,24,28,55,34,27,51,57)'
    parsed_ants = a.scripting.parse_ants(test_ants,8)
    ij = []
    for i in xrange(len(parsed_ants)):
        ij.append(a.miriad.bl2ij(parsed_ants[i][0]))
    ij = np.array(ij)

    freqs = []
    for i in xrange(len(ij)):
        freqs.append(calcFreq(array=aa, ij=ij[i], ref_ij=ij[len(ij)-1], ref_freq=freq, min_freq=0.100, max_freq=0.200))
    freqs = np.array(freqs)

    # import array parameters
    aa = a.cal.get_aa(calfile, freqs) # equatorial
    aa_freqs = aa.get_freqs()
    print "array parameters imported"

    #t = np.array([-3.0,-2.0,-1.0,0.0,1.0,2.0,3.0])
    t = np.array([0.0])
    os.system('rm -rf '+sim_dir+fileout)
    sim_data = a.miriad.UV(sim_dir + fileout, status='new')
    for k in xrange(len(ij)):
        gsmMap = makeGSM(path=gsm_dir, nside=64, freq=freqs[k])
        print "GSM made at freq %f" % freqs[k]
        for n in xrange(len(t)): 
            gsm_vis = calcVis(aa=aa, hpm=flatSky, ij=ij[k], freq=freqs[k], ha=t[n])
            i,j = ij[k]
            uvw = aa.gen_uvw(i,j,src='z')
            uvw = np.array([uvw[0,0,k], uvw[1,0,k], uvw[2,0,k]])
            vis_data = np.ma.array([gsm_vis], mask=False)
            preamble = (uvw, aa.get_jultime() + t[n], (i,j))
            print preamble
            print gsm_vis
            sim_data.write(preamble, vis_data)
            print "data written"

