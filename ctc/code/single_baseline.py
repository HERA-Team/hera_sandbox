## NAME: 
#         single_baseline.py
## PURPOSE: 
#         Visibility simulator for a single baseline for all frequencies for 1 Sidereal day

#------------------------------------------------------------------------------------------------

import aipy
import numpy
import pylab
import pyfits
import matplotlib.pyplot as plt


#get antenna array

filename = 'psa898_v003'

#get 3D image

img3d = aipy.map.Map(fromfits = '/Users/carinacheng/Desktop/Carina/UCBResearch/images/lambda_haslam408_dsds_eq.fits') #reads in 3D image; default is nside=512 (3145728 pixels)

#rescale frequency (Haslam map is 408GHz)

r = 150.0/408
f = r**-2.5 #synchrotron emission (spectral index -2.5)

img3d.map.map *= f #3D image data rescaled

#get eq coordinates

baseline = 3000/aipy.const.len_ns #ns
px = numpy.arange(img3d.npix()) #number of pixels in map
crd3d = numpy.array(img3d.px2crd(px,ncrd=3)) #aipy.healpix.HealpixMap.px2crd?
x3d,y3d,z3d = crd3d[0], crd3d[1], crd3d[2] #1D arrays of eq coordinates of 3Dimg (can define to be whatever coordinate system, but eq is most useful here)

#miriad uv file set-up

uv = aipy.miriad.UV('/Users/carinacheng/Desktop/Carina/UCBResearch/tables/single_baseline.uv', status='new')

uv.add_var('telescop' ,'a');   uv['telescop'] = 'AIPY'
uv.add_var('operator' ,'a');   uv['operator'] = 'AIPY'
uv.add_var('version'  ,'a');   uv['version'] = '0.0.1'
uv.add_var('epoch'    ,'r');   uv['epoch'] = 2000.
uv.add_var('sdf'      ,'d');   uv['sdf'] = 0.1/256
uv.add_var('sfreq'    ,'d');   uv['sfreq'] = 0.1
uv.add_var('freq'     ,'d');   uv['freq'] = 0.1
uv.add_var('restfreq' ,'d');   uv['restfreq'] = 0.1
uv.add_var('nchan'    ,'i');   uv['nchan'] = 256
uv.add_var('nschan'   ,'i');   uv['nschan'] = 256
uv.add_var('inttime'  ,'r');   uv['inttime'] = 10.0
uv.add_var('npol'     ,'i');   uv['npol'] = 1
uv.add_var('nspect'   ,'i');   uv['nspect'] = 1
uv.add_var('nants'    ,'i');   uv['nants'] = 32

aa = aipy.cal.get_aa(filename, uv['sdf'],uv['sfreq'],uv['nchan']) #don't know if this is correct?

uv.add_var('latitud'  ,'d');   uv['latitud'] = aa.lat
uv.add_var('dec'      ,'d');   uv['dec'] = aa.lat
uv.add_var('obsdec'   ,'d');   uv['obsdec'] = aa.lat
uv.add_var('longitu'  ,'d');   uv['longitu'] = aa.long
uv.add_var('antpos'   ,'d');   uv['antpos'] = (numpy.array([ant.pos for ant in aa], dtype = numpy.double)).transpose().flatten() #transpose is miriad convention

#variables to be updated

uv.add_var('coord'    ,'d')
uv.add_var('time'     ,'d')
uv.add_var('lst'      ,'d')
uv.add_var('ra'       ,'d')
uv.add_var('obsra'    ,'d')
uv.add_var('baseline' ,'r')
uv.add_var('pol'      ,'i')

#parameters used in loops

times = numpy.arange(2454500., 2454500.001, uv['inttime']/aipy.const.s_per_day)
freqs = numpy.arange(0.1,0.2,uv['sdf'])

for ii, t in enumerate(times):

    data = []

    print 'Timestep %d/%d' %(ii+1, len(times))
    aa.set_jultime(t)
    uv['time'] = t
    uv['lst'] = aa.sidereal_time()
    uv['ra'] = aa.sidereal_time()
    uv['obsra'] = aa.sidereal_time()

    t3d = aipy.coord.eq2top_m(aa.sidereal_time(),aa.lat)
    tx3d, ty3d, tz3d = numpy.dot(t3d,crd3d) #topocentric coordinates
    #bm3d = aa[0].bm_response((tx3d,ty3d,tz3d)) #beam response (makes code slow)
    #bm3d = numpy.where(tz3d < 0, 0, bm3d) #gets rid of beam values below horizon
    #sum_bm3d = numpy.sum(bm3d)
    
    #XXX east-west baseline only
    for jj, f in enumerate(freqs):
        
        fringe3d = numpy.exp(-2j*numpy.pi*tx3d*baseline*f) #fringe pattern

        fluxes3d = img3d.map.map

        p13d = fluxes3d*fringe3d#*bm3d
        sump13d = numpy.sum(p13d)#/sum_bm3d
        data.append(sump13d)

    data = numpy.asarray(data)

    for i, ai in enumerate(aa):
        for j, aj in enumerate(aa):
            if j < i: continue
            crd = aa.get_baseline(i,j)
            preamble = (crd, t, (i,j))
            uv['pol'] = aipy.miriad.str2pol['xx']
            flags = numpy.zeros((uv['nchan'],),dtype=numpy.int32)
            uv.write(preamble, data, flags)

del(uv)


