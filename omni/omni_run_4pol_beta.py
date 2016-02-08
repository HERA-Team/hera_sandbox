#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys

o = optparse.OptionParser()
o.set_usage('omni_run.py [options] *uvcRRE')
o.set_description(__doc__)
o.add_option('--calparXX',dest='calparXX',type='string',
            help='Path and name of .p calpar file for xx pol.')
o.add_option('--calparYY',dest='calparYY',type='string',
            help='Path and name of .p calpar file for yy pol.')
o.add_option('--redinfo',dest='redinfo',type='string',default='',
            help='Path and name of .bin redundant info file.')
o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
            help='Path and name of calfile.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .npz files. Include final / in path.')
o.add_option('--ba',dest='ba',default=None,
            help='Antennas to exclude, separated by commas.')
opts,args = o.parse_args(sys.argv[1:])

print 'Reading',opts.calpar
calparXX = pickle.load(open(opts.calparXX,'rb')) #firstcal initial gains xx
calparYY = pickle.load(open(opts.calparYY,'rb')) #firstcal initial gains yy



if opts.redinfo != '': #reading redinfo file
    print 'Reading',opts.redinfo
    info = omnical.info.RedundantInfoLegacy()
    print '   Getting reds from redundantinfo'
    info.fromfile(opts.redinfo)
    
else: #generate reds from calfile
    aa = aipy.cal.get_aa(opts.cal,numpy.array([.15]))
    print 'Getting reds from calfile'
    if opts.ba:
        ex_ants = []
        for a in opts.ba.split(','):
            ex_ants.append(int(a))
        print '   Excluding antennas:',ex_ants
    else: ex_ants = []
    #Can manually list them below if wanted
    
    info = capo.omni.aa_to_info(aa, ex_ants=ex_ants)

reds = info.get_reds()
reds_02 = omnical.arrayinfo.filter_reds(reds, ubls=[(64,10)]) #30m E/W baselines
#print ','.join(['%d_%d' % (i,j) for i,j in reds_02[0]])


print 'Getting gains from calpar'
#XXX Eventually, it would be nice if calpar had gains with keys 'x' or 'y', instead of 'xx'/'yy'


g0_x, g0_y = {}, {} 
for k in calparXX.keys(): #loop over pol
    g0_x[k] = {}
    for i in xrange(calparXX[k].shape[1]): #loop over antennas
        g0_x[k][i] = calparXX[k][:,i] / numpy.abs(calparXX[k][:,i]) #gains have len(nfreq)

for k in calparYY.keys():
    g0_y[k] = {}
    for i in xrange(calparYY[k].shape[1]):
        g0_y[k][i] = calparYY[k][:,i] / numpy.abs(calparYY[k][:,i])







### Loop Through Compressed Files ###
for f,filename in enumerate(args):
    print 'Reading', filename
    pol = filename.split('.')[-2] #XXX assumes 1 pol per file
    t,d,f = capo.arp.get_dict_of_uv_data([filename],antstr='cross',polstr=pol) #XXX could try reading only bls in reds
    SH = d.values()[0].values()[0].shape #shape of file data (ex: (19,203))
    data,wgts,xtalk = {}, {}, {}
    m2,g2,v2 = {}, {}, {}
    for i in g0[pol]: g0[pol][i] = numpy.resize(g0[pol][i],SH) #resize gains like data
    for bl in d: 
        i,j = aipy.miriad.bl2ij(bl)
        data[(i,j)] = d[bl][pol]
        wgts[(j,i)] = wgts[(i,j)] = numpy.logical_not(f[bl][pol]).astype(numpy.int)
    print '   Logcal-ing' 
    m1,g1,v1 = omnical.calib.redcal(data,info,gains=g0[pol])
    print '   Lincal-ing'
    m2[pol],g2[pol[0]],v2[pol] = omnical.calib.redcal(data, info, gains=g1, vis=v1, uselogcal=False, removedegen=True)
    xtalk[pol] = capo.omni.compute_xtalk(m2[pol]['res'], wgts) #xtalk is time-average of residual
    print '   Saving '+opts.omnipath+filename.split('/')[-1]+'.npz'
    capo.omni.to_npz(opts.omnipath+filename.split('/')[-1]+'.npz', m2, g2, v2, xtalk)
    
