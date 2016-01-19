#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys
import matplotlib.pyplot as plt


###
# This script runs logcal and lincal from Omnical to find a single complex bandpass solution
# Instead of REDS being baseline pairs, REDS is frequency channel pairs
# Instead of solving for gains per antenna, we solve for one number per frequency
# DATA is a dictionary indexed by frequency channel pairs, containing one complex64 (!) number for each pair: This number is the time-averaged correlation between two visibilities at two different frequencies that measure the same k-mode: <v1 v2*>
###

o = optparse.OptionParser()
o.set_usage('omni_run_passband.py [options] *npz')
o.set_description(__doc__)
o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
            help='Path and name of calfile.')
#o.add_option('--omnipath',dest='omnipath',default='',type='string',
#            help='Path to save .npz files. Include final / in path.')
o.add_option('--plot',dest='plot',default=False,action="store_true",
            help='Plot bandpass.')
o.add_option('--save',dest='save',default=False,action="store_true",
            help='Save bandpass to npz file.')
opts,args = o.parse_args(sys.argv[1:])

for fl in args:
    #Read NPZ file
    filename = fl.split('.npz')[0]
    print 'Reading', fl
    npz = numpy.load(fl)

    antpos = numpy.array(numpy.linspace(0,202,203),dtype=int)

    #reds = [[(1, 2),(0, 1),(2, 3)]] #freq pair channels
    #antpos = numpy.array([0,1,2,3])
    #data[(0,1)] = numpy.array([[2.+2.j]],dtype=numpy.complex64)
    #data[(1,2)] = numpy.array([[3.+3.j]],dtype=numpy.complex64)
    #data[(2,3)] = numpy.array([[1.+1.j]],dtype=numpy.complex64)

    data = {}
    datakeys = [] #tuples
    ugroups = [] 
    for k,key in enumerate(npz.keys()): #key must be tuple for data, but string for npz
        if key != 'files':
            ant1,ant2=key.split('(')[1].split(')')[0].split(',')
            ant2 = ant2.split(' ')[1]
            ant1 = int(ant1)
            ant2 = int(ant2)
            datakey = tuple((ant1,ant2))
            datakeys.append(datakey)
            data[datakey] = numpy.array([[npz[key][-1]]])
            ugroups.append(npz[key][1])

    print '   Creating reds'
    numgroups = numpy.max(ugroups)+1 #number of redundant groups
    reds =[[] for i in range(numgroups)]
    for k,key in enumerate(datakeys):
        reds[npz[str(key)][1]].append(key) #fill reds 

    """
    #get rid of flagged channels
    badones = []
    for k,key in enumerate(datakeys): #key must be tuple for data, but it's a string for npz
        if npz[str(key)][-1] == 0: #if flagged channel, don't save data
            badones.append(k)
        else:
            data[key] = numpy.array([[npz[str(key)][-1]]])  
    reds[0] = [i for j, i in enumerate(reds[0]) if j not in badones]
    """

    print '   Creating info'
    info = omnical.info.RedundantInfo()
    info.init_from_reds(reds,antpos)

    """
    #plot u's
    freqs_ghz = []
    bls = []
    aa = aipy.cal.get_aa('psa6622_v003',0.001,0.1,203)
    bl_str,bl_conj,bl2sep_str = capo.zsa.grid2ij(aa.ant_layout)
    freqs = numpy.linspace(0.1,0.2,203,endpoint=False)
    for k,key in enumerate(datakeys):
        chan = key[0]
        freqs_ghz.append(freqs[chan])
        a1,a2 = npz[str(key)][0][0]
        len_ns =  int(bl2sep_str[aipy.miriad.ij2bl(a1,a2)].split(',')[1])*15*10**9
        bls.append(len_ns)
    us = numpy.array(freqs_ghz)*numpy.array(bls)/(3*10**8)
    plt.hist(us,200)
    #plt.plot(us,'k.')
    plt.show()
    """

    #uCal
    print '   Logcal-ing'
    m1,g1,v1 = omnical.calib.redcal(data,info)
    print '   Lincal-ing'
    m2,g2,v2 = omnical.calib.redcal(data,info,gains=g1,vis=v1,uselogcal=False,removedegen=True)

    #import IPython;IPython.embed()

    #Plot bandpass
    xs = []
    ys = []
    phs = []
    for k,key in enumerate(g2.keys()):
        xs.append(key)    
        ys.append(numpy.abs(g2[key][0][0]))
        phs.append(numpy.angle(g2[key][0][0]))
    fit = numpy.polyval(numpy.polyfit(xs,ys,15),xs) #fit
    if opts.plot:
        plt.subplot(1,2,1)
        plt.plot(xs,ys,'k-')
        #plt.plot(xs,fit,'r-')
        plt.xlabel('Frequency Channel')
        plt.title('Amplitude')
        plt.subplot(1,2,2)
        plt.plot(xs,phs,'k-')
        plt.xlabel('Frequency Channel')
        plt.title('Phase')
        plt.show()

    #Save bandpass
    if opts.save:
        print '   Saving', filename+'bandpass.npz'
        d_npz = {}
        d_npz['chans'] = numpy.array(numpy.linspace(0,202,203),dtype=int)
        bp = numpy.zeros(len(d_npz['chans']))
        for chan in d_npz['chans']:
            try: bp[chan] = numpy.abs(g2[chan][0][0])
            except: bp[chan] = 1.0 #if no solution, set bandpass to 1
        d_npz['bandpass'] = bp
        d_npz['files'] = npz['files']
        numpy.savez(filename+'bandpass.npz',**d_npz)

