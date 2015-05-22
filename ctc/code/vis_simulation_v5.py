#!/usr/bin/env python

"""

NAME: 
      vis_simulation_v5.py 
PURPOSE:
      Models visibilities of maps of different polarizations, julian dates, and frequencies (adapted for Gianni Bernardi)
EXAMPLE CALL:
      vis_simulation_v5.py --mappath /data4/paper/exchange_sa_upenn/carina/ --filename test.uv -a 0_16 -C psa898_v003 `python -c "import numpy; file=numpy.load('julian_dates.npz'); print ' '.join(map(str,file['jd']))"`
AUTHOR:
      Carina Cheng

"""

import aipy
import numpy
import optparse
import os, sys

o = optparse.OptionParser()
o.set_usage('vis_simulation_v5.py [options] *.uv')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,ant=True)
o.add_option('--mappath', dest='mappath', default='/data4/paper/exchange_sa_upenn/Carina/',
             help='Directory where maps are. Include final / when typing path.')
o.add_option('--filename', dest='filename', default='/Users/carinacheng/capo/ctc/tables/testpspec.uv',
             help='Filename of created Miriad UV file (ex: test.uv).')

opts, args = o.parse_args(sys.argv[1:])

i,j = map(int,opts.ant.split('_'))
times = map(float,args) #converts args to floats
assert(not os.path.exists(opts.filename)) #checks if UV file exists already

freq_file = numpy.load(opts.mappath+'frequency.npz')
freqs = freq_file['freq']/(10**9)
sdf = freqs[1]-freqs[0]
sfreq = freqs[0]

JD_file = numpy.load(opts.mappath+'julian_dates.npz')
jds = JD_file['jd']

nchan = len(freqs)
inttime = (jds[1]-jds[0])*24*60*60 #seconds

#print freqs, jds, nchan, inttime

print 'getting antenna array...'

aa = aipy.cal.get_aa(opts.cal, sdf, freqs[0], nchan)
bl = aa.get_baseline(i,j) #for antennas 0 and 16; array of length 3 in ns
blx,bly,blz = bl[0],bl[1],bl[2]

#get topocentric coordinates 

print 'getting coordinates...'

img1 = aipy.map.Map(fromfits = opts.mappath + 'healpix_maps-I-f' + str(int(freqs[0]*10**3)) + '_j' + ('%.5f' % jds[0]) + '.fits', interp=True)

px = numpy.arange(img1.npix()) #number of pixels in map
crd = numpy.array(img1.px2crd(px,ncrd=3)) #aipy.healpix.HealpixMap.px2crd?
t3 = numpy.asarray(crd)
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
uvgridI = numpy.zeros(shape, dtype=numpy.complex64)
uvgridQ = numpy.zeros(shape, dtype=numpy.complex64)
uvgridU = numpy.zeros(shape, dtype=numpy.complex64)
uvgridV = numpy.zeros(shape, dtype=numpy.complex64)

for jj, f in enumerate(freqs):
    fng = numpy.exp(-2j*numpy.pi*(blx*tx+bly*ty+blz*tz)*f) #fringe pattern
    
    print 'Frequency %d/%d' % (jj+1, len(freqs)) 

    for ii, t in enumerate(times):

        print '   Timestep %d/%d' % (ii+1, len(times))
        aa.set_jultime(t)

        imgI = aipy.map.Map(fromfits = opts.mappath + 'healpix_maps-I-f' + str(int(f*10**3)) + '_j' + ('%.5f' % t) + '.fits', interp=True)
        imgQ = aipy.map.Map(fromfits = opts.mappath + 'healpix_maps-Q-f' + str(int(f*10**3)) + '_j' + ('%.5f' % t) + '.fits', interp=True)
        imgU = aipy.map.Map(fromfits = opts.mappath + 'healpix_maps-U-f' + str(int(f*10**3)) + '_j' + ('%.5f' % t) + '.fits', interp=True)
        imgV = aipy.map.Map(fromfits = opts.mappath + 'healpix_maps-V-f' + str(int(f*10**3)) + '_j' + ('%.5f' % t) + '.fits', interp=True)

        visI = numpy.sum(imgI[px]*fng)
        visQ = numpy.sum(imgQ[px]*fng)
        visU = numpy.sum(imgU[px]*fng)
        visV = numpy.sum(imgV[px]*fng)
        
        uvgridI[ii,jj] = visI
        uvgridQ[ii,jj] = visQ
        uvgridU[ii,jj] = visU
        uvgridV[ii,jj] = visV

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
uv.add_var('nchan' ,'i'); uv['nchan'] = nchan
uv.add_var('sdf' ,'d'); uv['sdf'] = sdf
uv.add_var('sfreq' ,'d'); uv['sfreq'] = sfreq
uv.add_var('freq' ,'d'); uv['freq'] = sfreq
uv.add_var('restfreq' ,'d'); uv['restfreq'] = sfreq
uv.add_var('nschan' ,'i'); uv['nschan'] = uv['nchan']
uv.add_var('inttime' ,'r'); uv['inttime'] = inttime
uv.add_var('npol' ,'i'); uv['npol'] = 4
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
    uv['pol'] = aipy.miriad.str2pol['I']
    uv.write(preamble, uvgridI[ii], flags[ii])
    uv['pol'] = aipy.miriad.str2pol['Q']
    uv.write(preamble, uvgridQ[ii], flags[ii])
    uv['pol'] = aipy.miriad.str2pol['U']
    uv.write(preamble, uvgridU[ii], flags[ii])
    uv['pol'] = aipy.miriad.str2pol['V']
    uv.write(preamble, uvgridV[ii], flags[ii])

del(uv)

