#!/usr/bin/env python

"""

NAME: 
      vis_simulation_v4_test.py 
PURPOSE:
      Used for testing vis_simulation_v4.py
EXAMPLE CALL: 
      ./vis_simulation_v4_test.py --sdf 0.001 --sfreq 0.1 --nchan 10 --inttime 20000 --filename test.uv -a 0_16 `python -c "import numpy; import aipy; print ' '.join(map(str,numpy.arange(2454500,2454501,20000/aipy.const.s_per_day)))"` -C psa898_v003
AUTHOR:
      Carina Cheng

"""

import aipy
import numpy
import pylab
import pyfits
import healpy
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import ephem as e
import optparse
import os, sys

o = optparse.OptionParser()
o.set_usage('vis_simulation_v4_test.py [options]')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,ant=True)
o.add_option('--filename', dest='filename', default='/Users/carinacheng/capo/ctc/tables/testpspec.uv',
             help='Filename of created Miriad UV file (ex: test.uv).')
o.add_option('--nchan', dest='nchan', default=1, type='int',
             help='Number of channels in simulated data. Default is 1.')
o.add_option('--inttime', dest='inttime', default=10., type='float',
             help='Integration time (s). Default is 10.') 
o.add_option('--sfreq', dest='sfreq', default=0.1, type='float',
             help='Start frequency (GHz). Default is 0.1.')
o.add_option('--sdf', dest='sdf', default=0.1/203, type='float',
             help='Channel spacing (GHz).  Default is .1/203')
opts, args = o.parse_args(sys.argv[1:])

i,j = map(int,opts.ant.split('_'))
times = map(float,args) #converts args to floats

assert(not os.path.exists(opts.filename)) #checks if UV file exists already

print 'getting antenna array...'

aa = aipy.cal.get_aa(opts.cal, opts.sdf, opts.sfreq, opts.nchan)
freqs = aa.get_afreqs()
bl = aa.get_baseline(i,j) #for antennas 0 and 16; array of length 3 in ns
blx,bly,blz = bl[0],bl[1],bl[2]

#get topocentric coordinates and calculate beam response


### MAKE MAP WITH ONE PIXEL ###
img = aipy.map.Map(nside=512)
img.map.map = numpy.zeros_like(img.map.map) #empty map
value = 1.0
wgts = 1.0

px = numpy.arange(img.npix()) #number of pixels in map
crd = numpy.array(img.px2crd(px,ncrd=3)) #aipy.healpix.HealpixMap.px2crd?
t3 = numpy.asarray(crd)
tx,ty,tz = t3[0],t3[1],t3[2]

#get galactic coordinates from map

g3 = numpy.asarray(crd)

aa.set_jultime(times[0])
txi,tyi,tzi = 1,0,0
top2eq = aipy.coord.top2eq_m(-aa.sidereal_time(), aa.lat) #note the minus sign
exi,eyi,ezi = numpy.dot(top2eq,(txi,tyi,tzi)) #equatorial coordinates
#exi,eyi,ezi = 0,0,-1 #south pole
eq2ga = aipy.coord.convert_m('ga','eq') #input/output mixed up??
gxi,gyi,gzi = numpy.dot(eq2ga,(exi,eyi,ezi)) #galactic coordinates
img.put((gxi,gyi,gzi),wgts,value)  


#####

def herabeam(filename,top):
    class HealpixMap(aipy.healpix.HealpixMap):
    ### class update to get HERA fits beams working ###
        def from_fits(self, filename, hdunum=1, colnum=0):
            hdu = pyfits.open(filename)[hdunum]
            data = hdu.data.field(colnum)
            if not data.dtype.isnative:
                data.dtype = data.dtype.newbyteorder()
                data.byteswap(True)
            data = data.flatten()
            scheme = hdu.header['ORDERING'][:4]
            self.set_map(data, scheme=scheme)
    hmap = HealpixMap(nside=128)
    hmap.from_fits(filename)
    return hmap[(top[0],top[1],top[2])]

print 'getting maps and calculating fringes...'

#loop through time to pull out fluxes and fringe pattern
#loop through frequency to calculate visibility

shape = (len(times),len(freqs))
flags = numpy.zeros(shape, dtype=numpy.int32)
uvgridxx = numpy.zeros(shape, dtype=numpy.complex64)
uvgridyy = numpy.zeros(shape, dtype=numpy.complex64)

for jj, f in enumerate(freqs):
    fng = numpy.exp(-2j*numpy.pi*(blx*tx+bly*ty+blz*tz)*f) #fringe pattern
    aa.select_chans([jj]) #selects specific frequency
    bmxx = aa[0].bm_response((t3[0],t3[1],t3[2]), pol='x')[0]**2 ### PAPER beam
    #bmxx = herabeam('/Users/carinacheng/capo/ctc/code/hera-cst/HERA_DISH_paper_feed_cyl36_150mhz_X_healpix.fits',t3) ### HERA beam
    
    """        
    #Plot Beam in Topocentric (looking down on observer)
    im = aipy.img.Img(800,.5)
    size=1600
    h = aipy.healpix.HealpixMap(nside=512)
    h.map = bmxx
    x,y,z = im.get_top(center=(size/2,size/2))
    v = numpy.logical_not(x.mask)
    d = h[x.flatten(),y.flatten(),z.flatten()]
    d.shape = (size,size)
    d = numpy.where(v,d,numpy.NaN)
    m = Basemap(projection='ortho',lat_0=aa.lat,lon_0=aa.long,rsphere=1.)
    m.imshow(d.real,interpolation='bicubic',origin='lower',cmap='jet')
    plt.show()
    """ 
    
    bmyy = aa[0].bm_response((t3[0],t3[1],t3[2]), pol='y')[0]**2
    sum_bmxx = numpy.sum(bmxx)
    sum_bmyy = numpy.sum(bmyy)
    fngxx = fng*bmxx/sum_bmxx #factor used later in visibility calculation
    fngyy = fng*bmyy/sum_bmyy
    fluxes = img[px] #fluxes preserved in galactic grid

    print 'Frequency %d/%d' % (jj+1, len(freqs)) 

    toplot1 = numpy.zeros(len(times))
    toplot2 = numpy.zeros(len(times))
    plt.figure(figsize=(10,8))
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9, wspace=0.4, hspace=0.3)
    pylab.ion()
    plt1 = None
    
    for ii, t in enumerate(times):

        print '   Timestep %d/%d' % (ii+1, len(times))
        aa.set_jultime(t)
        
        ga2eq = aipy.coord.convert_m('eq','ga',iepoch=aa.epoch,oepoch=e.J2000) #conversion matrix
        eq2top = aipy.coord.eq2top_m(-aa.sidereal_time(),aa.lat) #conversion matrix, note the minus sign
        ga2eq2top = numpy.dot(eq2top,ga2eq) #topocentric coordinates
        t3rot = numpy.dot(ga2eq2top,g3)
        txrot = numpy.ma.compressed(numpy.ma.masked_where(t3rot[2]<0,t3rot[0]))
        tyrot = numpy.ma.compressed(numpy.ma.masked_where(t3rot[2]<0,t3rot[1]))
        tzrot = numpy.ma.compressed(numpy.ma.masked_where(t3rot[2]<0,t3rot[2])) #mask coordinates below horizon
        fluxes2 = numpy.ma.compressed(numpy.ma.masked_where(t3rot[2]<0,fluxes)) #mask data below horizon

        pxrot, wgts = img.crd2px(txrot,tyrot,tzrot, interpolate=1) 

        efngxx = numpy.sum(fngxx[pxrot]*wgts, axis=1)
        efngyy = numpy.sum(fngyy[pxrot]*wgts, axis=1)
        visxx = numpy.sum(fluxes2*efngxx)
        visyy = numpy.sum(fluxes2*efngyy)
        toplot1[ii] = numpy.sum(fluxes2)
        toplot2[ii] = numpy.real(visxx)

        uvgridxx[ii,jj] = visxx
        uvgridyy[ii,jj] = visyy

        if plt1 == None:
            plt.subplot(2,1,1)
            plt1 = plt.plot(times,toplot1,'b-')
            plt.xlabel('Time')
            plt.ylabel('Sum(Fluxes)')
            plt.ylim(0,2)
            plt.subplot(2,1,2)
            plt2 = plt.plot(times,toplot2,'r-')
            #plt3 = plt.plot(times,toplot3,'b-')
            plt.xlabel('Time')
            plt.ylabel('Real(Vis)')
            plt.ylim(-1e-6,1e-6)
            pylab.show()
        else:
            plt1[0].set_ydata(toplot1)
            plt2[0].set_ydata(toplot2)
            #plt3[0].set_ydata(toplot3)
            pylab.draw()

    plt.savefig('/Users/carinacheng/capo/ctc/code/time_axis.png')
    print ("%.8f" % f) + ' GHz done'


"""
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
