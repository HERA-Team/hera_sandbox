#!/usr/bin/env python

"""

NAME:
 single_baseline_test.py
PURPOSE:
 Visibility simulator for a single baseline for all frequencies for 1 Sidereal day
 Different from single_baseline.py in that fringe/beam is computed outside of the loop

------------------------------------------------------------------------------------------------
"""


import aipy
import numpy
import pylab
import pyfits
import matplotlib.pyplot as plt

root = '/Users/carinacheng/capo/ctc/'

print 'setting up miriad UV file...'

#miriad uv file set-up

uv = aipy.miriad.UV(root + 'tables/single_baseline_test.uv', status='new')

uv.add_var('telescop' ,'a'); uv['telescop'] = 'AIPY'
uv.add_var('operator' ,'a'); uv['operator'] = 'AIPY'
uv.add_var('version' ,'a'); uv['version'] = '0.0.1'
uv.add_var('epoch' ,'r'); uv['epoch'] = 2000.
uv.add_var('nchan' ,'i'); uv['nchan'] = 203
uv.add_var('sdf' ,'d'); uv['sdf'] = 0.1/uv['nchan']
uv.add_var('sfreq' ,'d'); uv['sfreq'] = 0.1
uv.add_var('freq' ,'d'); uv['freq'] = 0.1
uv.add_var('restfreq' ,'d'); uv['restfreq'] = 0.1
uv.add_var('nschan' ,'i'); uv['nschan'] = uv['nchan']
uv.add_var('inttime' ,'r'); uv['inttime'] = 10.0
uv.add_var('npol' ,'i'); uv['npol'] = 1
uv.add_var('nspect' ,'i'); uv['nspect'] = 1
uv.add_var('nants' ,'i'); uv['nants'] = 32

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

#freqs = numpy.asarray([0.1,0.12,0.16,0.2]) #if only doing certain freqs

filename = 'psa898_v003'
#aa = aipy.cal.get_aa(filename, freqs) #if only doing certain freqs
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

img1 = aipy.map.Map(fromfits = root + 'images/gsm/gsm256/gsm1001.fits')

px = numpy.arange(img1.npix()) #number of pixels in map
crd = numpy.array(img1.px2crd(px,ncrd=3)) #aipy.healpix.HealpixMap.px2crd?
t3 = numpy.asarray(crd)
t3 = t3.compress(t3[2]>=0,axis=1)
tx,ty,tz = t3[0], t3[1], t3[2] #1D arrays of top coordinates of 3Dimg (can define to be whatever coordinate system)
bm = aa[0].bm_response((t3[0],t3[1],t3[2])) #beam response (makes code slow)
sum_bm = numpy.sum(bm,axis=1)

#loop through frequency to get GSM and calculate fringe

print 'getting GSMs and calculating fringes...'

img = {}
fng = {}
for jj, f in enumerate(freqs):
    img[f] = aipy.map.Map(fromfits = root + 'images/gsm/gsm256/gsm1' + str(jj+1).zfill(3) + '.fits', interp=True)
    #XXX east-west baseline only
    fng[f] = numpy.exp(-2j*numpy.pi*(blx*tx+bly*ty+blz*tz)*f) #fringe pattern
    fng[f]*= bm[jj]/sum_bm[jj]
    print str(f) + ' GHz done'

#loop through time to pull out fluxes from changing equatorial coordinates
#loop through frequency to calculate visibility

print 'writing miriad uv file...'

times = numpy.arange(2454500., 2454501., 61*uv['inttime']/aipy.const.s_per_day)
flags = numpy.zeros(len(freqs), dtype=numpy.int32)

for ii, t in enumerate(times):

    print 'Timestep %d/%d' %(ii+1, len(times))
    aa.set_jultime(t)
    uv['time'] = t
    uv['lst'] = aa.sidereal_time()
    uv['ra'] = aa.sidereal_time()
    uv['obsra'] = aa.sidereal_time()

    top2eq = aipy.coord.top2eq_m(aa.sidereal_time(),aa.lat)
    ex,ey,ez = numpy.dot(top2eq,t3) #equatorial coordinates

    data = []

    px, wgts = img[freqs[0]].crd2px(ex,ey,ez, interpolate=1)

    for jj, f in enumerate(freqs):

        fluxes = numpy.sum(img[f][px]*wgts, axis=1)

        vis = numpy.sum(fluxes*fng[f])#/sum_bm
        data.append(vis)

        print str(f) + ' GHz done'
    
    data = numpy.asarray(data)
   
    preamble = (bl, t, (i,j))
    uv['pol'] = aipy.miriad.str2pol['xx']
    uv.write(preamble, data, flags)

del(uv)

