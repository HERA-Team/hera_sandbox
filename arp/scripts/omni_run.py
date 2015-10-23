#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys

o = optparse.OptionParser()
o.set_usage('omni_run.py [options] *uvcRRE')
o.set_description(__doc__)
o.add_option('--calpar',dest='calpar',type='string',
            help='Path and name of .p calpar file.')
#o.add_option('--redinfo',dest='redinfo',type='string',
#            help='Path and name of .bin redundant info file.')
o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
            help='Path and name of calfile.')
#o.add_option('--ubls',dest='ubls',default=None,
#            help='Unique baselines to include. Ex: [(64,49),(64,10)]')
#o.add_option('--ex_ubls',dest='ex_ubls',default=None,
#            help='Unique baselines to exclude. Ex: [(64,49),(64,10)]')
#o.add_option('--ants',dest='ants',default=None,
#            help='Antennas to include. Ex: [64,49,10]')
#o.add_option('--ex_ants',dest='ex_ants',default=None,
#            help='Antennas to exclude. Ex: [64,49,10]')
opts,args = o.parse_args(sys.argv[1:])

calpar = pickle.load(open(opts.calpar,'rb')) # the firstcal initial gains

if False:
    print 'Reading calpar and redinfo files...'
    info = omnical.info.RedundantInfoLegacy() #reading old firstcal files
    info.fromfile(opts.redinfo)
else:
    aa = aipy.cal.get_aa(opts.cal,numpy.array([.15]))
    ex_ants = [72,53] + [48,82,83,37] + [107,15,16] + [92,7] + [27,50,54] + [26,38,110,46]
    #ex_ants += [84,76,68,57,29,28,24,17,5,63,8]
    ex_ants += [29,24,28,55] + [68,17,69,13] + [76,5,77,32] + [84,40,85,14] # are all these bad?
    ex_ants += [57,63,8,74]
    info = capo.omni.aa_to_info(aa, ex_ants=ex_ants)

# if you want to plot fiducial baselines to check omnical sol, this prints out what's been calibrated
reds = info.get_reds()
reds_02 = omnical.arrayinfo.filter_reds(reds, ubls=[(64,10)])
print ','.join(['%d_%d' % (i,j) for i,j in reds_02[0]])

#XXX Eventually, it would be nice if calpar had gains with keys 'x' or 'y', instead of 'xx'/'yy'
g0 = {} 
for k in calpar.keys():
    for i in xrange(calpar[k].shape[1]):
        g0[i] = calpar[k][:,i] / numpy.abs(calpar[k][:,i])

### Loop Through Compressed Files ###
for f,filename in enumerate(args):
    print 'Reading', filename
    pol = filename.split('.')[-2] # XXX assumes 1 pol per file
    t,d,f = capo.arp.get_dict_of_uv_data([filename],antstr='cross',polstr=pol) # XXX could try reading only bls in reds
    SH = d.values()[0].values()[0].shape
    data,wgts,xtalk = {}, {}, {}
    m2,g2,v2 = {}, {}, {}
    for i in g0: g0[i] = numpy.resize(g0[i], SH)
    for bl in d:
        i,j = aipy.miriad.bl2ij(bl)
        data[(i,j)] = d[bl][pol]
        wgts[(j,i)] = wgts[(i,j)] = numpy.logical_not(f[bl][pol]).astype(numpy.int) # in case res is switched
    print '   logcal-ing' 
    m1,g1,v1 = omnical.calib.redcal(data,info,gains=g0) #logcal
    print '   lincal-ing'
    m2[pol],g2[pol[0]],v2[pol] = omnical.calib.redcal(data, info, gains=g1, vis=v1, uselogcal=False, removedegen=True)
    xtalk[pol] = capo.omni.compute_xtalk(m2[pol]['res'], wgts)
    capo.omni.to_npz(filename+'.npz', m2, g2, v2, xtalk)
