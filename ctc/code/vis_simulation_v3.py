#!/usr/bin/env python

"""

NAME: 
      vis_simulation_v3.py 
PURPOSE:
      Models visibilities using power spectra (pspecs) from pspec_sim_v2.py and creates a new Miriad UV file
      Differs from vis_simulation.py in that the sky image uses eq. coordinates and the fringe/beam is rotated with time (interpolation happens for fringe)
EXAMPLE CALL:
      ./vis_simulation_v3.py --nchan 2 --inttime 30000 
IMPORTANT NOTE: 
      Be careful when changing sdf and sfreq because they need to match the pspec files!
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

o = optparse.OptionParser()
o.set_usage('vis_simulation_v3.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--map', dest='map', default='/Users/carinacheng/capo/ctc/images/pspecs/pspec40lmax110/',
             help='Directory where PSPEC files are (labeled pspec1001.fits, pspec1002.fits, etc.). Include final / when typing path.')
o.add_option('--filename', dest='filename', default='/Users/carinacheng/capo/ctc/tables/testpspec.uv',
             help='Filename of created Miriad UV file (ex: test.uv).')
o.add_option('--nchan', dest='nchan', default=203, type='int',
             help='Number of channels in simulated data. Default is 203.')
o.add_option('--inttime', dest='inttime', default=10., type='float',
             help='Integration time (s). Default is 10.') 
o.add_option('--sfreq', dest='sfreq', default=0.1, type='float',
             help='Start frequency (GHz). Default is 0.1.')
o.add_option('--sdf', dest='sdf', default=0.1/203, type='float',
             help='Channel spacing (GHz).  Default is .1/203')
o.add_option('--startjd', dest='startjd', default=2454500., type='float',
             help='Julian Date to start observation.  Default is 2454500')
o.add_option('--endjd', dest='endjd', default=2454501., type='float',
             help='Julian Date to end observation.  Default is 2454501')
o.add_option('--psa', dest='psa', default='psa898_v003', 
             help='Name of PSA file.')
opts, args = o.parse_args(sys.argv[1:])

#miriad uv file set-up

print 'setting up miriad UV file...'

uv = aipy.miriad.UV(opts.filename, status='new')
uv._wrhd('obstype','mixed-auto-cross')
uv._wrhd('history','MDLVIS: created file.\nMDLVIS: ' + ' '.join(sys.argv) + '\n')

uv.add_var('telescop' ,'a'); uv['telescop'] = 'AIPY'
uv.add_var('operator' ,'a'); uv['operator'] = 'AIPY'
uv.add_var('version' ,'a'); uv['version'] = '0.0.1'
uv.add_var('epoch' ,'r'); uv['epoch'] = 2000.
uv.add_var('source'  ,'a'); uv['source'] = 'zenith'
uv.add_var('nchan' ,'i'); uv['nchan'] = opts.nchan
uv.add_var('sdf' ,'d'); uv['sdf'] = opts.sdf
uv.add_var('sfreq' ,'d'); uv['sfreq'] = opts.sfreq
uv.add_var('freq' ,'d'); uv['freq'] = opts.sfreq
uv.add_var('restfreq' ,'d'); uv['restfreq'] = opts.sfreq
uv.add_var('nschan' ,'i'); uv['nschan'] = uv['nchan']
uv.add_var('inttime' ,'r'); uv['inttime'] = opts.inttime
uv.add_var('npol' ,'i'); uv['npol'] = 1
uv.add_var('nspect' ,'i'); uv['nspect'] = 1
uv.add_var('nants' ,'i'); uv['nants'] = 32

#variables just set to dummy values

uv.add_var('vsource' ,'r'); uv['vsource'] = 0.
uv.add_var('ischan'  ,'i'); uv['ischan'] = 1
uv.add_var('tscale'  ,'r'); uv['tscale'] = 0.
uv.add_var('veldop'  ,'r'); uv['veldop'] = 0.

#variables to be updated

uv.add_var('coord' ,'d')
uv.add_var('time' ,'d')
uv.add_var('lst' ,'d')
uv.add_var('ra' ,'d')
uv.add_var('obsra' ,'d')
uv.add_var('baseline' ,'r')
uv.add_var('pol' ,'i')

#get antenna array

print 'getting antenna array...'

filename = opts.psa
aa = aipy.cal.get_aa(filename, uv['sdf'],uv['sfreq'],uv['nchan'])
freqs = aa.get_afreqs()
i = 0
j = 16
bl = aa.get_baseline(i,j) #for antennas 0 and 16; array of length 3 in ns
blx,bly,blz = bl[0],bl[1],bl[2]

#more miriad variables

uv.add_var('latitud' ,'d'); uv['latitud'] = aa.lat
uv.add_var('dec' ,'d'); uv['dec'] = aa.lat
uv.add_var('obsdec' ,'d'); uv['obsdec'] = aa.lat
uv.add_var('longitu' ,'d'); uv['longitu'] = aa.long
uv.add_var('antpos' ,'d'); uv['antpos'] = (numpy.array([ant.pos for ant in aa], dtype = numpy.double)).transpose().flatten() #transpose is miriad convention

#get topocentric coordinates and calculate beam response

print 'calculating beam response...'

img1 = aipy.map.Map(fromfits = opts.map + 'pspec1001.fits', interp=True)

px = numpy.arange(img1.npix()) #number of pixels in map
crd = numpy.array(img1.px2crd(px,ncrd=3)) #aipy.healpix.HealpixMap.px2crd?
t3 = numpy.asarray(crd)
#t3 = t3.compress(t3[2]>=0,axis=1) #gets rid of coordinates below horizon
tx,ty,tz = t3[0], t3[1], t3[2] #1D arrays of top coordinates of img1 (can define to be whatever coordinate system)
bmxx = aa[0].bm_response((t3[0],t3[1],t3[2]), pol='x')**2 #beam response (makes code slow)
bmyy = aa[0].bm_response((t3[0],t3[1],t3[2]), pol='y')**2
sum_bmxx = numpy.sum(bmxx,axis=1)
sum_bmyy = numpy.sum(bmyy,axis=1)

#get equatorial coordinates

e3 = numpy.asarray(crd)
ex,ey,ez = e3[0], e3[1], e3[2] #1D arrays of eq coordinates of img

#loop through frequency to get PSPECS and calculate fringe

print 'getting PSPECs and calculating fringes...'

img = {}
fngxx = {}
fngyy = {}
fluxes = {}

for jj, f in enumerate(freqs):
    img[f] = aipy.map.Map(fromfits = opts.map + 'pspec1' + str(jj+1).zfill(3) + '.fits', interp=True)
    fng = numpy.exp(-2j*numpy.pi*(blx*tx+bly*ty+blz*tz)*f) #fringe pattern
    fngxx[f] = fng*bmxx[jj]/sum_bmxx[jj] #factor used later in visibility calculation
    fngyy[f] = fng*bmyy[jj]/sum_bmyy[jj]
    fluxes[f] = img[f][px] #fluxes preserved in equatorial grid
    print ("%.8f" % f) + ' GHz done'

#loop through time to pull out fluxes and fringe pattern
#loop through frequency to calculate visibility

print 'writing miriad uv file...'

times = numpy.arange(opts.startjd, opts.endjd, uv['inttime']/aipy.const.s_per_day)
flags = numpy.zeros(len(freqs), dtype=numpy.int32)

for ii, t in enumerate(times):

    print 'Timestep %d/%d' %(ii+1, len(times))
    aa.set_jultime(t)
    uv['time'] = t
    uv['lst'] = aa.sidereal_time()
    uv['ra'] = aa.sidereal_time()
    uv['obsra'] = aa.sidereal_time()

    eq2top = aipy.coord.eq2top_m(aa.sidereal_time(),aa.lat) #conversion matrix
    t3 = numpy.dot(eq2top,e3) #topocentric coordinates
    #t3 = t3.compress(t3[2]>=0,axis=1) #gets rid of coordinates below horizon
    tx,ty,tz = t3[0], t3[1], t3[2] 

    dataxx = []
    datayy = []

    px, wgts = img[freqs[0]].crd2px(tx,ty,tz, interpolate=1) #converts coordinates to pixels for first PSPEC file (pixel numbers don't change with frequency)
    #NOTE: img and fluxes are still in equatorial coordinates... just getting pixels here

    for jj, f in enumerate(freqs):

        efngxx = numpy.sum(fngxx[f][px]*wgts, axis=1)
        efngyy = numpy.sum(fngyy[f][px]*wgts, axis=1)
        
        visxx = numpy.sum(fluxes[f]*efngxx)
        visyy = numpy.sum(fluxes[f]*efngyy)
        dataxx.append(visxx)
        datayy.append(visyy)

        print ("%.8f" % f) + ' GHz done'
    
    dataxx = numpy.asarray(dataxx)
    datayy = numpy.asarray(datayy)
   
    preamble = (bl, t, (i,j))
    uv['pol'] = aipy.miriad.str2pol['xx']
    uv.write(preamble, dataxx, flags)
    uv['pol'] = aipy.miriad.str2pol['yy']
    uv.write(preamble, datayy, flags)

del(uv)
