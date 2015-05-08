#!/usr/bin/env python

"""

NAME: 
      vis_simulation_v5.py 
PURPOSE:
      Models visibilities of maps of different polarizations, julian dates, and frequencies
EXAMPLE CALL:
      ./vis_simulation_v5.py --sdf 0.001 --sfreq 0.1 --nchan 10 --inttime 20000 --mappath /data4/paper/exchange_sa_upenn/Carina/ --filename test.uv -a 0_16 `python -c "import numpy; import aipy; print ' '.join(map(str,numpy.arange(2454500,2454501,20000/aipy.const.s_per_day)))"` -C psa898_v003
AUTHOR:
      Carina Cheng

"""

import aipy
import numpy
#import pylab
#import pyfits
#import matplotlib.pyplot as plt
import optparse
import os, sys

o = optparse.OptionParser()
o.set_usage('vis_simulation_v5.py [options] *.uv')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,ant=True)
o.add_option('--mappath', dest='mappath', default='/data4/paper/exchange_sa_upenn/Carina/',
             help='Directory where maps are. Include final / when typing path.')
#o.add_option('--map', dest='map', default='gsm',
#            help='Map type (gsm or pspec).')
o.add_option('--filename', dest='filename', default='/Users/carinacheng/capo/ctc/tables/testpspec.uv',
             help='Filename of created Miriad UV file (ex: test.uv).')
#o.add_option('--nchan', dest='nchan', default=203, type='int',
#             help='Number of channels in simulated data. Default is 203.')
#o.add_option('--inttime', dest='inttime', default=10., type='float',
#             help='Integration time (s). Default is 10.') 
#o.add_option('--sfreq', dest='sfreq', default=0.1, type='float',
#             help='Start frequency (GHz). Default is 0.1.')
#o.add_option('--sdf', dest='sdf', default=0.1/203, type='float',
#             help='Channel spacing (GHz).  Default is .1/203')

opts, args = o.parse_args(sys.argv[1:])

i,j = map(int,opts.ant.split('_'))
times = map(float,args) #converts args to floats

assert(not os.path.exists(opts.filename)) #checks if UV file exists already

freq_file = numpy.load(opts.mappath+'frequency.npz')
freqs = freq_file['freq']/(10**9)
sdf = freqs[1]-freqs[0]

JD_file = numpy.load(opts.mappath+'judian_dates.npz')
jds = JD_file['jd']

nchan = len(freqs)
inttime = (jds[1]-jds[0])*24*60*60 #seconds

print freqs, jds, nchan, inttime

print 'getting antenna array...'

aa = aipy.cal.get_aa(opts.cal, sdf, freqs[0], nchan)
bl = aa.get_baseline(i,j) #for antennas 0 and 16; array of length 3 in ns
blx,bly,blz = bl[0],bl[1],bl[2]

#get topocentric coordinates 

print 'getting coordinates...'

img1 = aipy.map.Map(fromfits = opts.mappath + 'healpix_maps-I-f' + ('%.2f' % (float(freqs[0])*10**3)) + '_l' + ('%.5f' % jds[0]) + '.fits', interp=True)

px = numpy.arange(img1.npix()) #number of pixels in map
crd = numpy.array(img1.px2crd(px,ncrd=3)) #aipy.healpix.HealpixMap.px2crd?
t3 = numpy.asarray(crd)
#t3 = t3.compress(t3[2]>=0,axis=1) #gets rid of coordinates below horizon
tx,ty,tz = t3[0], t3[1], t3[2] #1D arrays of top coordinates of img1 (can define to be whatever coordinate system)

#get equatorial coordinates

e3 = numpy.asarray(crd)
ex,ey,ez = e3[0], e3[1], e3[2] #1D arrays of eq coordinates of img

#loop through frequency to get maps and calculate fringe

print 'getting maps and calculating fringes...'

#loop through time to pull out fluxes and fringe pattern
#loop through frequency to calculate visibility

shape = (len(times),len(freqs))
flags = numpy.zeros(shape, dtype=numpy.int32)
uvgridxx = numpy.zeros(shape, dtype=numpy.complex64)
uvgridyy = numpy.zeros(shape, dtype=numpy.complex64)

#NOTE TO SELF: I STOPPED EDITING AT THIS POINT... FITS FILES HAVE ERROR BEING READ BY AIPY
#NEED TO GET RID OF BEAM STUFF BELOW AND SET IT UP FOR 4 POLS

"""
for jj, f in enumerate(freqs):
    img = aipy.map.Map(fromfits = opts.mappath+opts.map + '1' + str(jj+1).zfill(3) + '.fits', interp=True)
    fng = numpy.exp(-2j*numpy.pi*(blx*tx+bly*ty+blz*tz)*f) #fringe pattern
    aa.select_chans([jj]) #selects specific frequency
    bmxx = aa[0].bm_response((t3[0],t3[1],t3[2]), pol='x')**2
    bmyy = aa[0].bm_response((t3[0],t3[1],t3[2]), pol='y')**2
    sum_bmxx = numpy.sum(bmxx,axis=1)
    sum_bmyy = numpy.sum(bmyy,axis=1)
    fngxx = fng*bmxx[0]/sum_bmxx[0] #factor used later in visibility calculation
    fngyy = fng*bmyy[0]/sum_bmyy[0]
    fluxes = img[px] #fluxes preserved in equatorial grid

    print 'Frequency %d/%d' % (jj+1, len(freqs)) 

    for ii, t in enumerate(times):

        print '   Timestep %d/%d' % (ii+1, len(times))
        aa.set_jultime(t)

        eq2top = aipy.coord.eq2top_m(aa.sidereal_time(),aa.lat) #conversion matrix
        t3rot = numpy.dot(eq2top,e3) #topocentric coordinates
        #t3 = t3.compress(t3[2]>=0,axis=1) #gets rid of coordinates below horizon
        txrot,tyrot,tzrot = t3rot[0], t3rot[1], t3rot[2] 

        pxrot, wgts = img.crd2px(txrot,tyrot,tzrot, interpolate=1) #converts coordinates to pixels for first PSPEC file (pixel numbers don't change with frequency)
        #NOTE: img and fluxes are still in equatorial coordinates... just getting pixels here

        efngxx = numpy.sum(fngxx[pxrot]*wgts, axis=1)
        efngyy = numpy.sum(fngyy[pxrot]*wgts, axis=1)
        
        visxx = numpy.sum(fluxes*efngxx)
        visyy = numpy.sum(fluxes*efngyy)

        uvgridxx[ii,jj] = visxx
        uvgridyy[ii,jj] = visyy

    print ("%.8f" % f) + ' GHz done'

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

#more miriad variables

uv.add_var('latitud' ,'d'); uv['latitud'] = aa.lat
uv.add_var('dec' ,'d'); uv['dec'] = aa.lat
uv.add_var('obsdec' ,'d'); uv['obsdec'] = aa.lat
uv.add_var('longitu' ,'d'); uv['longitu'] = aa.long
uv.add_var('antpos' ,'d'); uv['antpos'] = (numpy.array([ant.pos for ant in aa], dtype = numpy.double)).transpose().flatten() #transpose is miriad convention

for ii, t in enumerate(times):

    print '%d/%d' % (ii+1, len(times))+' done'
    aa.set_jultime(t)
    uv['time'] = t
    uv['lst'] = aa.sidereal_time()
    uv['ra'] = aa.sidereal_time()
    uv['obsra'] = aa.sidereal_time()    

    preamble = (bl, t, (i,j))
    uv['pol'] = aipy.miriad.str2pol['xx']
    uv.write(preamble, uvgridxx[ii], flags[ii])
    uv['pol'] = aipy.miriad.str2pol['yy']
    uv.write(preamble, uvgridyy[ii], flags[ii])

del(uv)
"""
