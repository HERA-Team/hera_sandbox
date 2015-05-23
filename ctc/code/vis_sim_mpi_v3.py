#!/usr/bin/env python

"""

NAME: 
   vis_sim_mpi.py
PURPOSE:
   Models visibilities of Healpix maps and creates a new Miriad UV file
EXAMPLE CALL:
   mpiexec -np 3 vis_sim_mpi.py --nchan 100 --inttime 10 --sfreq 0.1 --sdf 0.001 -C psa898_v003 --map gsm --mappath /Users/carinacheng/capo/ctc/images/gsm/gsm100/
IMPORTANT NOTE:
   Make sure sdf and sfreq options match those of the maps!
AUTHOR:
   Carina Cheng

"""

import aipy
import numpy
import pyfits
import optparse
import os, sys
from mpi4py import MPI

#user options

o = optparse.OptionParser()
o.set_usage('vis_sim_mpi.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--mappath', dest='mappath', default='/Users/carinacheng/capo/ctc/images/gsm/gsm203/',
             help='Directory where maps are. Include final / when typing path.')
o.add_option('--map', dest='map', default='gsm',
             help='Map type (gsm or pspec).')
o.add_option('--filename', dest='filename', default='/Users/carinacheng/capo/ctc/tables/test.uv',
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
o.add_option('-C', dest='psa', default='psa898_v003', 
             help='Name of calfile.')
o.add_option('--bli', dest='bli', default=0,
             help='Baseline i. Default is 0.')
o.add_option('--blj', dest='blj', default=16,
             help='Baseline j. Default is 16.')
opts, args = o.parse_args(sys.argv[1:])

#MPI set-up

comm = MPI.COMM_WORLD #get MPI communicator object
size = comm.size      #total number of processors
rank = comm.rank      #rank of a process
status = MPI.Status() #MPI status object (contains source and tag)

if rank == 0:

    #miriad uv file set-up

    print '***Master is setting up miriad UV file...'

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

    print '***Master is getting antenna array...'

    calfile = opts.psa
    aa = aipy.cal.get_aa(calfile,uv['sdf'],uv['sfreq'],uv['nchan'])
    freqs = aa.get_afreqs()
    i = opts.bli
    j = opts.blj
    bl = aa.get_baseline(i,j) #array of length 3 in ns

    #more miriad variables
    
    uv.add_var('latitud','d'); uv['latitud'] = aa.lat
    uv.add_var('dec','d'); uv['dec'] = aa.lat
    uv.add_var('obsdec','d'); uv['obsdec'] = aa.lat
    uv.add_var('longitu','d'); uv['longitu'] = aa.long
    uv.add_var('antpos','d'); uv['antpos'] = (numpy.array([ant.pos for ant in aa],dtype=numpy.double)).transpose().flatten() #transpose is miriad convention

    #setting up tasks

    task_index = 0
    num_workers = size-1
    closed_workers = 0
    count = 1

    times = numpy.arange(opts.startjd, opts.endjd, uv['inttime']/aipy.const.s_per_day)
    flags = numpy.zeros(len(freqs), dtype=numpy.int32)

    all_dataxx = numpy.zeros((len(times),len(freqs)),dtype=numpy.complex)
    #all_datayy = numpy.zeros((len(times),len(freqs)),dtype=numpy.complex)

    print '***Master is starting with %d workers and %d time integrations...' % (num_workers, len(times))

    while closed_workers < num_workers:
        worker_data = comm.recv(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG,status=status)
        tag = status.Get_tag()
        source = status.Get_source()
        if tag == 4: #worker ready to get beam/fringe
            comm.send(None,dest=source,tag=5)
        elif tag == 0: #worker ready to rotate sky
            if task_index < len(times):
                print 'Sending timestep %d/%d to processor %d' % (task_index+1,len(times),source)
                t = times[task_index]
                aa.set_jultime(t)
                uv['time'] = t
                uv['lst'] = aa.sidereal_time()
                uv['ra'] = aa.sidereal_time()
                uv['obsra'] = aa.sidereal_time()
                sid_time = float(aa.sidereal_time())
                lat = float(aa.lat)
                comm.send((task_index,sid_time,lat,freqs),dest=source,tag=3)
                task_index += 1
            else:
                comm.send(None,dest=source,tag=2) #exit tag
        elif tag == 1: #done tag
            index = worker_data[0]
            all_dataxx[index] = worker_data[1]
            #all_datayy[index] = worker_data[?]
            print 'Timestep %d/%d complete' % (count,len(times))
            count += 1
        elif tag == 2: #no more workers needed tag
            print 'Processor %d done' % source
            closed_workers += 1

    print '***Master finishing by writing UV file...'
    for kk, t in enumerate(times): #write UV file sequentially
        preamble = (bl,t,(i,j))
        uv['pol'] = aipy.miriad.str2pol['xx']
        uv.write(preamble,all_dataxx[kk],flags)
        #uv['pol'] = aipy.miriad.str2pol['yy']
        #uv.write(preamble,all_datayy[kk],flags)

else: #other processors do this

    print 'Processor %d is ready' % rank

    comm.send(None,dest=0,tag=4)
    task = comm.recv(source=0,tag=MPI.ANY_TAG,status=status)
    tag = status.Get_tag()
    if tag == 5: #ready to get beam/fringe

        print 'Processor %d getting beam/fringe' % rank

        calfile = opts.psa
        aa = aipy.cal.get_aa(calfile,opts.sdf,opts.sfreq,opts.nchan)
        freqs = aa.get_afreqs()
        i = opts.bli
        j = opts.blj
        bl = aa.get_baseline(i,j) #array of length 3 in ns
        blx,bly,blz = bl[0],bl[1],bl[2] #single numbers

        img1 = aipy.map.Map(fromfits = opts.mappath + opts.map + '1001.fits', interp=True)
        px = numpy.arange(img1.npix()) #number of pixels in map
        crd = numpy.array(img1.px2crd(px,ncrd=3))
        t3 = numpy.asarray(crd)
        tx,ty,tz = t3[0],t3[1],t3[2] #1D topocentric coordinates
        bmxx = aa[0].bm_response((t3[0],t3[1],t3[2]),pol='x')**2
        #bmyy = aa[0].bm_response((t3[0],t3[1],t3[2]),pol='y')**2
        sum_bmxx = numpy.sum(bmxx,axis=1)
        #sum_bmyy = numpy.sum(bmyy,axis=1)
        e3 = numpy.asarray(crd)

        """

        fngxx = {}
        #fngyy = {}
        fluxes = {}

        for jj, f in enumerate(freqs):
            img = aipy.map.Map(fromfits = opts.mappath+opts.map+'1'+str(jj+1).zfill(3)+'.fits', interp=True)
            fng = numpy.exp(-2j*numpy.pi*(blx*tx+bly*ty+blz*tz)*f) #fringe pattern
            fngxx[f] = fng*bmxx[jj]/sum_bmxx[jj]
            #fngyy[f] = fng*bmyy[jj]/sum_bmyy[jj]
            fluxes[f] = img[px] #fluxes preserved in equatorial grid

        """

    while True:
        comm.send(None,dest=0,tag=0) #ready tag
        task = comm.recv(source=0,tag=MPI.ANY_TAG,status=status)
        tag = status.Get_tag()
        if tag == 3: #start tag

            index = task[0]
            sid_time,lat = task[1],task[2]
            freqs = task[3]

            eq2top = aipy.coord.eq2top_m(sid_time,lat) #conversion matrix
            t3rot = numpy.dot(eq2top,e3) #topocentric coordinates
            txrot,tyrot,tzrot = t3rot[0],t3rot[1],t3rot[2]

            dataxx = []
            #datayy = []

            img1 = aipy.map.Map(fromfits = opts.mappath+opts.map+'1001.fits',interp=True)
            pxrot,wgts = img1.crd2px(txrot,tyrot,tzrot,interpolate=1)

            for jj, f in enumerate(freqs):
                #"""
                img = aipy.map.Map(fromfits = opts.mappath+opts.map+'1'+str(jj+1).zfill(3)+'.fits', interp=True)
                fng = numpy.exp(-2j*numpy.pi*(blx*tx+bly*ty+blz*tz)*f)
                fngxx = fng*bmxx[jj]/sum_bmxx[jj]
                #fngyy = fng*bmyy[jj]/sum_bmyy[jj]
                fluxes = img[px]

                efngxx = numpy.sum(fngxx[pxrot]*wgts, axis=1)
                #efngyy = numpy.sum(fngyy[pxrot]*wgts, axis=1)
                visxx = numpy.sum(fluxes*efngxx)
                #visyy = numpy.sum(fluxes*efngyy)
                #"""
                """
                efngxx = numpy.sum(fngxx[f][pxrot]*wgts, axis=1)
                #print fngxx[f].shape,fluxes[f].shape,efngxx.shape
                #efngyy = numpy.sum(fngyy[f][pxrot]*wgts, axis=1)
                visxx = numpy.sum(fluxes[f]*efngxx)
                #visyy = numpy.sum(fluxes[f]*efngyy)
                """
                dataxx.append(visxx)
                #datayy.append(visyy)

            dataxx = numpy.asarray(dataxx)
            #datayy = numpy.asarray(datayy)

            comm.send((index,dataxx),dest=0,tag=1)

        elif tag == 2: #exit tag
            break
    comm.send(None,dest=0,tag=2)


