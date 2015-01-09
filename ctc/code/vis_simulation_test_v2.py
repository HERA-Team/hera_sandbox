#!/usr/bin/env python

"""

NAME: 
      vis_simulation_test_v2.py 
PURPOSE:
      -Models visibilities using selected simple test cases and creates a new Miriad UV file
      -Doesn't read in a map (one source vector pointing towards flux of 1)
EXAMPLE CALL:
      ./vis_simulation_test_v2.py  
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

#options

o = optparse.OptionParser()
o.set_usage('vis_simulation_test.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--map', dest='map', default='/Users/carinacheng/capo/ctc/images/gsm/gsm203/',
             help='Directory where GSM files are (labeled gsm1001.fits, gsm1002.fits, etc.). Include final / when typing path.')
o.add_option('--filename', dest='filename', default='/Users/carinacheng/capo/ctc/tables/test.uv',
             help='Filename of created Miriad UV file (ex: test.uv).')
o.add_option('--nchan', dest='nchan', default=203, type='int',
             help='Number of channels in simulated data. Default is 203.')
o.add_option('--case', dest='case', default='pole',
             help='Test case. Can be "east" (source rises on the East) or "pole" (source at South Pole) or "zenith" (source at observers zenith). Default is "pole".')
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

filename = 'psa898_v003'
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

#cases

if opts.case == 'pole':

    exi,eyi,ezi = 0,0,-1

if opts.case == 'east':

    aa.set_jultime(opts.startjd)
    txi,tyi,tzi = -1,0,0 #east?
    top2eq = aipy.coord.top2eq_m(aa.sidereal_time(), aa.lat)
    exi,eyi,ezi = numpy.dot(top2eq,(txi,tyi,tzi)) #equatorial coordinates
    exi,eyi,ezi = 0.43171066, -0.90201214, 6.12323400e-17 #correcting to match with vis_simulation_test.py

if opts.case == 'zenith':

    aa.set_jultime(opts.startjd)
    txi,tyi,tzi = 0,0,1 #zenith
    top2eq = aipy.coord.top2eq_m(aa.sidereal_time(), aa.lat)
    exi,eyi,ezi = numpy.dot(top2eq,(txi,tyi,tzi)) #equatorial coordinates

print 'calculating beam...'

#beam (only used to find sum_bm)

img = aipy.map.Map(nside=512)
px = numpy.arange(img.npix())
crd = numpy.array(img.px2crd(px,ncrd=3))
t3 = numpy.asarray(crd)
bmi = aa[0].bm_response(t3,pol='x')
sum_bm = numpy.sum(bmi,axis=1)

#loop through time
#loop through frequency to calculate visibility

print 'writing miriad uv file...'

times = numpy.arange(2454500.2,2454500.3,uv['inttime']/aipy.const.s_per_day)#opts.startjd, opts.endjd, uv['inttime']/aipy.const.s_per_day)
flags = numpy.zeros(len(freqs), dtype=numpy.int32)

#plotting set-up

plt.figure(figsize = (10,8))
plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.3)
pylab.ion()
plt1 = None
toplot = numpy.zeros(times.shape)
toplot2 = numpy.zeros(times.shape)

for ii, t in enumerate(times):

    print 'Timestep %d/%d' %(ii+1, len(times))
    aa.set_jultime(t)
    uv['time'] = t
    uv['lst'] = aa.sidereal_time()
    uv['ra'] = aa.sidereal_time()
    uv['obsra'] = aa.sidereal_time()

    eq2top = aipy.coord.eq2top_m(aa.sidereal_time(),aa.lat) #conversion matrix
    tx,ty,tz = numpy.dot(eq2top,(exi,eyi,ezi)) #top coordinates
    #could get rid of half the values here

    data = []

    for jj, f in enumerate(freqs):

        fng = numpy.exp(-2j*numpy.pi*(blx*tx+bly*ty+blz*tz)*f) #fringe pattern
        bm = aa[0].bm_response((tx,ty,tz), pol='x')
        fng *= bm[jj][0]/sum_bm[jj]

        toplot[ii] = 1.0
        
        vis = toplot[ii]*fng
        toplot2[ii] = numpy.real(vis)

        data.append(vis)

        print ("%.5f" % f) + ' GHz done'

        """
        if plt1 == None:

            plt.subplot(2,1,1)
            plt1 = pylab.plot(times,toplot,'b.')
            plt.xlabel("Time (Julian Date)")
            plt.ylabel("Sum Flux")
            plt.ylim(0,2)
            plt.subplot(2,1,2)
            plt2 = pylab.plot(times,toplot2,'r.')
            plt.xlabel("Time (Julian Date)")
            plt.ylabel("Visibilities")
            plt.ylim(-1e-6,1e-6)
            pylab.show()

        else:
            
            plt1[0].set_ydata(toplot)
            plt2[0].set_ydata(toplot2)
            pylab.draw()
    
        """

    data = numpy.asarray(data)
    #print data
   
    preamble = (bl, t, (i,j))
    uv['pol'] = aipy.miriad.str2pol['xx']
    uv.write(preamble, data, flags)

#plt.savefig('/Users/carinacheng/capo/ctc/test.png')#time_axis_' + opts.case + '.png')
del(uv)



