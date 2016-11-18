#! /usr/bin/env python
import aipy as a
import capo
import numpy as np
import os
from GlobalSkyModel import GlobalSkyModel

class AntennaArray(a.pol.AntennaArray):
    def geometric_dly(self, i, j):
        """ converts unique baseline numbers into nanosecond delays """
        bx, by, bz = self.get_baseline(i,j)
        return np.sqrt(bx**2 + by**2 + bz**2)        

def bl2delay(array, bl):
    """ converts unique baseline numbers into nanosecond delays """

def calcFreq(array, ij, ref_ij, ref_freq, min_freq, max_freq):
    """ calculates nuCal frequencies (in GHz) based on an initial frequency, reference baseline, 
        and input baseline. must include min and max frequencies for array in order to
        ensure that calculated frequencies are legal. """
<<<<<<< HEAD
    bl = bl2delay(array, ij)
    ref_bl = bl2delay(array, ref_ij)
    f = ref_freq * (ref_bl / bl)
    assert(max_freq >= f >= min_freq)
=======
    i, j = a.miriad.bl2ij(bl)
    bl = aa.get_baseline(i, j)
    bl = np.sqrt(np.dot(bl, bl))
    f = ref_freq * (ref_bl / baseline)
    assert(max_freq >= f >= min_freq) # check for invalid freq
>>>>>>> acc5d41f1b85003fefd93ae51bdd5e0a4d7bff48
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

<<<<<<< HEAD
def makeGSM(path, nside, freq):
        gsmMap = a.healpix.HealpixMap(nside=nside)
        g = a.healpix.HealpixMap(nside=512)
        g.map = GlobalSkyModel(freq=1000*freq, GSMlocation=path, GSMNSIDE=512).map # galactic
        g = hpm_TtoJy(g, freq=freq)
        gsmMap.from_hpm(g)
        return gsmMap
=======
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

class GSMMap(a.healpix.HealpixMap):
    def from_fits(self, filename):
        

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

>>>>>>> acc5d41f1b85003fefd93ae51bdd5e0a4d7bff48

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
<<<<<<< HEAD
    import matplotlib.pylab as pl
    import optparse, sys
=======
    import optparse, sys
    import matplotlib.pylab as pl
>>>>>>> acc5d41f1b85003fefd93ae51bdd5e0a4d7bff48

    o = optparse.OptionParser()
    o.add_option('--sim_dir', default='/home/kara/capo/kmk/scripts/')
    o.add_option('--data_dir', default='/home/kara/capo/kmk/data/')
    o.add_option('--data_file', default='zen.2456895.51490.xx.uvcRRE')
    o.add_option('--calfile', default='test')
    o.add_option('--gsm_dir', default='/home/kara/capo/kmk/gsm/gsm_raw/')
<<<<<<< HEAD
    o.add_option('--fileout', default='sim_results.uv')

=======
    
>>>>>>> acc5d41f1b85003fefd93ae51bdd5e0a4d7bff48
    opts,args = o.parse_args(sys.argv[1:])
    sim_dir = opts.sim_dir
    data_dir = opts.data_dir
    data_file = opts.data_file
    calfile = opts.calfile
    gsm_dir = opts.gsm_dir
<<<<<<< HEAD
    fileout = opts.fileout

    # select 150 MHz and 160 MHz for u-mode calibration test
    freq = 0.100
    aa = a.cal.get_aa(calfile, np.array([freq]))

    flatSky = makeFlatMap(nside=64, Tsky=1.0, freq=freq)
=======

    # select 150 MHz and 160 MHz for u-mode calibration test
    freq = 0.100

    flatSky = makeFlatMap(nside=64, Tsky=1.0, freq=freq)
    xyz = flatSky.px2crd(np.arange(flatSky.npix())) #topocentric
>>>>>>> acc5d41f1b85003fefd93ae51bdd5e0a4d7bff48

    # create array of baselines
    test_ants = '(64)_(29,24,28,55,34,27,51,57)'
    parsed_ants = a.scripting.parse_ants(test_ants,8)
<<<<<<< HEAD
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

=======
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
>>>>>>> acc5d41f1b85003fefd93ae51bdd5e0a4d7bff48
