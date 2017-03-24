#!/usr/bin/env python

"""

NAME: 
      vis_simulation_v4.py 
PURPOSE:
      Set-up for grid engine on folio
      Models visibilities using power spectra (pspecs) from pspec_sim_v2.py and creates a new Miriad UV file
      Differs from vis_simulation.py in that the sky image uses eq. coordinates and the fringe/beam is rotated with time (interpolation happens for fringe)
EXAMPLE CALL: 
      ./vis_simulation_v4.py --sdf 0.001 --sfreq 0.1 --nchan 10 --inttime 20000 --map pspec --mappath /Users/carinacheng/capo/ctc/images/pspecs/pspec100lmax100/ --filename test.uv -a 0_16 `python -c "import numpy; import aipy; print ' '.join(map(str,numpy.arange(2454500,2454501,20000/aipy.const.s_per_day)))"` -C psa898_v003
IMPORTANT NOTE: 
      Be careful when changing sdf and sfreq because they need to match the pspec files!
AUTHOR:
      Carina Cheng

"""

import aipy
import numpy
import ephem as e
import optparse
import os, sys

o = optparse.OptionParser()
o.set_usage('vis_simulation_v4.py [options] *.uv')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,ant=True)
o.add_option('--mappath', dest='mappath', default='/Users/carinacheng/capo/ctc/images/pspecs/pspec40lmax110/',
             help='Directory where maps are. Include final / when typing path.')
o.add_option('--map', dest='map', default='gsm',
            help='Map type (gsm or pspec).')
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
opts, args = o.parse_args(sys.argv[1:])

bls = opts.ant.split(',')
times = map(float,args) #converts args to floats

namexx = '.'.join(opts.filename.split('.')[:-1])+'.xx.uv'
nameyy = '.'.join(opts.filename.split('.')[:-1])+'.yy.uv'

if os.path.exists(namexx) or os.path.exists(nameyy):
    if os.path.exists(namexx):
        print '   %s exists. Skipping...' % namexx
    if os.path.exists(nameyy):
        print '   %s exists. Skipping...' % nameyy
    exit()

print 'getting antenna array...'

aa = aipy.cal.get_aa(opts.cal, opts.sdf, opts.sfreq, opts.nchan)
freqs = aa.get_afreqs()

#get topocentric coordinates and calculate beam response

print 'calculating beam response...'

img1 = aipy.map.Map(fromfits = opts.mappath+opts.map + '1001.fits', interp=True)

px = numpy.arange(img1.npix()) #number of pixels in map
crd = numpy.array(img1.px2crd(px,ncrd=3)) #aipy.healpix.HealpixMap.px2crd?
t3 = numpy.asarray(crd)
tx,ty,tz = t3[0], t3[1], t3[2] #1D arrays of top coordinates of img1 (can define to be whatever coordinate system)

g3 = numpy.asarray(crd) #map is in galactic coordinates

print 'getting maps and calculating fringes...'

shape = (len(times),len(freqs))
flags = {}
uvgridxx = {}
uvgridyy = {}
bls_pos = []

for bb,baseline in enumerate(bls):
    print 'Baseline %d/%d' % (bb+1, len(bls))
    i,j = map(int,baseline.split('_'))
    bl = aa.get_baseline(i,j) #[ns]
    bls_pos.append(bl)
    blx,bly,blz = bl[0],bl[1],bl[2]
    flags[(i,j)] = numpy.zeros(shape, dtype=numpy.int32)
    uvgridxx[(i,j)] = numpy.zeros(shape, dtype=numpy.complex64)
    uvgridyy[(i,j)] = numpy.zeros(shape, dtype=numpy.complex64)

    for jj, f in enumerate(freqs):
        print '   Frequency %d/%d' % (jj+1, len(freqs))
        img = aipy.map.Map(fromfits = opts.mappath+opts.map + '1' + str(jj+1).zfill(3) + '.fits', interp=True)
        #img = aipy.map.Map(fromfits = opts.mappath+opts.map + '1001.fits', interp=True) #reading same map over and over again
        fng = numpy.exp(-2j*numpy.pi*(blx*tx+bly*ty+blz*tz)*f) #fringe pattern
        aa.select_chans([jj]) #selects specific frequency
        bmxx = aa[0].bm_response((t3[0],t3[1],t3[2]), pol='x')[0]**2
        bmyy = aa[0].bm_response((t3[0],t3[1],t3[2]), pol='y')[0]**2
        sum_bmxx = numpy.sum(bmxx)
        sum_bmyy = numpy.sum(bmyy)
        fngxx = fng*bmxx/sum_bmxx #factor used later in visibility calculation
        fngyy = fng*bmyy/sum_bmyy
        fluxes = img[px] #fluxes preserved in galactic grid

        for ii, t in enumerate(times):

            print '      Timestep %d/%d' % (ii+1, len(times))
            aa.set_jultime(t)

            ga2eq = aipy.coord.convert_m('eq','ga',iepoch=e.J2000,oepoch=aa.epoch) #conversion matrix
            eq2top = aipy.coord.eq2top_m(-aa.sidereal_time(),aa.lat) #conversion matrix
            ga2eq2top = numpy.dot(eq2top,ga2eq)
            t3rot = numpy.dot(ga2eq2top,g3) #topocentric coordinates
            txrot = numpy.ma.compressed(numpy.ma.masked_where(t3rot[2]<0,t3rot[0]))
            tyrot = numpy.ma.compressed(numpy.ma.masked_where(t3rot[2]<0,t3rot[1]))
            tzrot = numpy.ma.compressed(numpy.ma.masked_where(t3rot[2]<0,t3rot[2])) #mask coordinates below horizon
            fluxes2 = numpy.ma.compressed(numpy.ma.masked_where(t3rot[2]<0,fluxes)) #mask data below horizon

            pxrot, wgts = img.crd2px(txrot,tyrot,tzrot, interpolate=1) 

            efngxx = numpy.sum(fngxx[pxrot]*wgts, axis=1)
            efngyy = numpy.sum(fngyy[pxrot]*wgts, axis=1)
            
            visxx = numpy.sum(fluxes2*efngxx)
            visyy = numpy.sum(fluxes2*efngyy)

            uvgridxx[(i,j)][ii,jj] = visxx
            uvgridyy[(i,j)][ii,jj] = visyy

        print ("      %.8f" % f) + ' GHz done'

#miriad uv file set-up

print 'setting up miriad UV file...'

uvxx = aipy.miriad.UV(namexx, status='new')
uvxx._wrhd('obstype','mixed-auto-cross')
uvxx._wrhd('history','MDLVIS: created file.\nMDLVIS: ' + ' '.join(sys.argv) + '\n')

uvyy = aipy.miriad.UV(nameyy, status='new')
uvyy._wrhd('obstype','mixed-auto-cross')
uvyy._wrhd('history','MDLVIS: created file.\nMDLVIS: ' + ' '.join(sys.argv) + '\n')

uvxx.add_var('telescop' ,'a'); uvxx['telescop'] = 'AIPY'
uvxx.add_var('operator' ,'a'); uvxx['operator'] = 'AIPY'
uvxx.add_var('version' ,'a'); uvxx['version'] = '0.0.1'
uvxx.add_var('epoch' ,'r'); uvxx['epoch'] = 2000.
uvxx.add_var('source'  ,'a'); uvxx['source'] = 'zenith'
uvxx.add_var('nchan' ,'i'); uvxx['nchan'] = opts.nchan
uvxx.add_var('sdf' ,'d'); uvxx['sdf'] = opts.sdf
uvxx.add_var('sfreq' ,'d'); uvxx['sfreq'] = opts.sfreq
uvxx.add_var('freq' ,'d'); uvxx['freq'] = opts.sfreq
uvxx.add_var('restfreq' ,'d'); uvxx['restfreq'] = opts.sfreq
uvxx.add_var('nschan' ,'i'); uvxx['nschan'] = uvxx['nchan']
uvxx.add_var('inttime' ,'r'); uvxx['inttime'] = opts.inttime
uvxx.add_var('npol' ,'i'); uvxx['npol'] = 1
uvxx.add_var('nspect' ,'i'); uvxx['nspect'] = 1
uvxx.add_var('nants' ,'i'); uvxx['nants'] = 128

uvyy.add_var('telescop' ,'a'); uvyy['telescop'] = 'AIPY'
uvyy.add_var('operator' ,'a'); uvyy['operator'] = 'AIPY'
uvyy.add_var('version' ,'a'); uvyy['version'] = '0.0.1'
uvyy.add_var('epoch' ,'r'); uvyy['epoch'] = 2000.
uvyy.add_var('source'  ,'a'); uvyy['source'] = 'zenith'
uvyy.add_var('nchan' ,'i'); uvyy['nchan'] = opts.nchan
uvyy.add_var('sdf' ,'d'); uvyy['sdf'] = opts.sdf
uvyy.add_var('sfreq' ,'d'); uvyy['sfreq'] = opts.sfreq
uvyy.add_var('freq' ,'d'); uvyy['freq'] = opts.sfreq
uvyy.add_var('restfreq' ,'d'); uvyy['restfreq'] = opts.sfreq
uvyy.add_var('nschan' ,'i'); uvyy['nschan'] = uvyy['nchan']
uvyy.add_var('inttime' ,'r'); uvyy['inttime'] = opts.inttime
uvyy.add_var('npol' ,'i'); uvyy['npol'] = 1
uvyy.add_var('nspect' ,'i'); uvyy['nspect'] = 1
uvyy.add_var('nants' ,'i'); uvyy['nants'] = 128

#variables just set to dummy values

uvxx.add_var('vsource' ,'r'); uvxx['vsource'] = 0.
uvxx.add_var('ischan'  ,'i'); uvxx['ischan'] = 1
uvxx.add_var('tscale'  ,'r'); uvxx['tscale'] = 0.
uvxx.add_var('veldop'  ,'r'); uvxx['veldop'] = 0.

uvyy.add_var('vsource' ,'r'); uvyy['vsource'] = 0.
uvyy.add_var('ischan'  ,'i'); uvyy['ischan'] = 1
uvyy.add_var('tscale'  ,'r'); uvyy['tscale'] = 0.
uvyy.add_var('veldop'  ,'r'); uvyy['veldop'] = 0.

#variables to be updated

uvxx.add_var('coord' ,'d')
uvxx.add_var('time' ,'d')
uvxx.add_var('lst' ,'d')
uvxx.add_var('ra' ,'d')
uvxx.add_var('obsra' ,'d')
uvxx.add_var('baseline' ,'r')
uvxx.add_var('pol' ,'i')

uvyy.add_var('coord' ,'d')
uvyy.add_var('time' ,'d')
uvyy.add_var('lst' ,'d')
uvyy.add_var('ra' ,'d')
uvyy.add_var('obsra' ,'d')
uvyy.add_var('baseline' ,'r')
uvyy.add_var('pol' ,'i')

#get antenna array

#more miriad variables

uvxx.add_var('latitud' ,'d'); uvxx['latitud'] = aa.lat
uvxx.add_var('dec' ,'d'); uvxx['dec'] = aa.lat
uvxx.add_var('obsdec' ,'d'); uvxx['obsdec'] = aa.lat
uvxx.add_var('longitu' ,'d'); uvxx['longitu'] = aa.long
uvxx.add_var('antpos' ,'d'); uvxx['antpos'] = (numpy.array([ant.pos for ant in aa], dtype = numpy.double)).transpose().flatten() #transpose is miriad convention

uvyy.add_var('latitud' ,'d'); uvyy['latitud'] = aa.lat
uvyy.add_var('dec' ,'d'); uvyy['dec'] = aa.lat
uvyy.add_var('obsdec' ,'d'); uvyy['obsdec'] = aa.lat
uvyy.add_var('longitu' ,'d'); uvyy['longitu'] = aa.long
uvyy.add_var('antpos' ,'d'); uvyy['antpos'] = (numpy.array([ant.pos for ant in aa], dtype = numpy.double)).transpose().flatten() #transpose is miriad convention

print 'Writing', namexx, ',', nameyy

for ii, t in enumerate(times):

    aa.set_jultime(t)
    uvxx['time'] = t
    uvxx['lst'] = aa.sidereal_time()
    uvxx['ra'] = aa.sidereal_time()
    uvxx['obsra'] = aa.sidereal_time()    
    uvyy['time'] = t
    uvyy['lst'] = aa.sidereal_time()
    uvyy['ra'] = aa.sidereal_time()
    uvyy['obsra'] = aa.sidereal_time()    

    for bb,bl in enumerate(uvgridxx.keys()):
        preamble = (bls_pos[bb], t, bl)
        uvxx['pol'] = aipy.miriad.str2pol['xx']
        uvxx.write(preamble, uvgridxx[bl][ii], flags[bl][ii])
        uvyy['pol'] = aipy.miriad.str2pol['yy']
        uvyy.write(preamble, uvgridyy[bl][ii], flags[bl][ii])

del(uvxx)
del(uvyy)
