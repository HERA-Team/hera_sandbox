#! /usr/bin/env python

import omnical
import aipy
import pylab
import numpy
import capo
import pickle
import optparse
import os, sys

### Options ###
o = optparse.OptionParser()
o.set_usage('omni_run.py [options] *uvcRRE')
o.set_description(__doc__)
o.add_option('--calpar',dest='calpar',type='string',
            help='Path and name of .p calpar file.')
#o.add_option('--redinfo',dest='redinfo',type='string',
#            help='Path and name of .bin redundant info file.')
o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
            help='Path and name of calfile.')
o.add_option('--omniruntag',dest='omniruntag',default='',type='string',
            help='Tag for omni run, if wanted. Default is empty.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .omni_output npz files. Include final / in path.')
o.add_option('--xtalk',dest='xtalk',default=False,action="store_true",
            help='Include xtalk in lincal.')
#o.add_option('--ubls',dest='ubls',default=None,
#            help='Unique baselines to include. Ex: [(64,49),(64,10)]')
#o.add_option('--ex_ubls',dest='ex_ubls',default=None,
#            help='Unique baselines to exclude. Ex: [(64,49),(64,10)]')
#o.add_option('--ants',dest='ants',default=None,
#            help='Antennas to include. Ex: [64,49,10]')
#o.add_option('--ex_ants',dest='ex_ants',default=None,
#            help='Antennas to exclude. Ex: [64,49,10]')
opts,args = o.parse_args(sys.argv[1:])

### Save Options ###
CALFILE = opts.cal
CALPAR = opts.calpar # aka firstcal
RED_INFO = opts.redinfo

def compute_xtalk(res, wgts):
    xtalk = {}
    for key in res:
        r,w = numpy.where(wgts[key] > 0, res[key], 0), wgts[key].sum(axis=0)
        w = numpy.where(w == 0, 1, w)
        xtalk[key] = (r.sum(axis=0) / w).astype(res[key].dtype) # avg over time
    return xtalk

##def resize_dict(d, shape):
#def resize_dict(xtalk, wgts):
#    for k in xtalk: xtalk[k] = (xtalk[k] * wgts[k]).astype(numpy.complex64)
#    return xtalk
#    #for k in d: d[k] = numpy.resize(d[k], shape)
#    #return d

### Read Firstcal Info ###


if False:
    print 'Reading calpar and redinfo files...'
    info = omnical.info.RedundantInfoLegacy() #reading old firstcal files
    info.fromfile(RED_INFO)
    reds = info.get_reds()
else:
    aa = aipy.cal.get_aa(CALFILE,numpy.array([.15]))
    antpos = -numpy.ones((len(aa),3))
    layout = aa.ant_layout
    xs,ys = numpy.indices(layout.shape)
    for i,x,y in zip(layout.flatten(), xs.flatten(), ys.flatten()):
        antpos[i,0],antpos[i,1] = x,y
    #ants = [64,10,49,3] + [65,9,66,58] + [72,22,73,61] + [80,20,81,63]
    #ants += [88,43,89,2] + [96,53,97,21] + [104,31,105,45]
    #ants += [41,25,19,48] + [67,1,47,4] + [74,35,75,18] + [82,42,83,37]
    #ants += [90,33,91,6] + [98,15,99,16] + [106,8,107,11] 
    #ants += [92,52,93,7] + [100,62,101,44] + [108,36,109,60]
    #ants += [34,27,51,57] + [70,56,71,59] + [78,30,79,23] + [86,54,78,50]
    #ants += [94,12,95,38] + [102,0,103,26] + [110,39,111,46]
    ex_ants = [72,53] + [48,82,83,37] + [107,15,16] + [92,7] + [27,50,54] + [26,38,110,46]
    #ex_ants += [29,24,28,55] + [68,17,69,13] + [76,5,77,32] + [84,40,85,14] # are all these bad?
    ex_ants += [29,24,28,55] + [68,17,69,13] + [76,5,77] + [84,40,85,14] # are all these bad?
    #ex_ants = [i for i in range(antpos.shape[0]) if antpos[i,0] < 0]
    reds = omnical.arrayinfo.compute_reds(antpos,tol=.1) #bypasses creating AI class
    reds = omnical.arrayinfo.filter_reds(reds, ex_ants=ex_ants)
    #reds = omnical.arrayinfo.filter_reds(reds, ants=ants, ex_ants=ex_ants)
    info = omnical.info.RedundantInfo() #make RI class
    info.init_from_reds(reds,antpos)

reds_02 = omnical.arrayinfo.filter_reds(reds, ubls=[(64,10)])
print ','.join(['%d_%d' % (i,j) for i,j in reds_02[0]])
#antpos = numpy.array([a.pos for a in aa]) #nanoseconds, equatorial
#antpos = numpy.array(antpos) #antpos needs to be (#ant, 3) array with positions
#antpos = antpos*aipy.const.len_ns/100. #meters
#antpos = [aa.get_baseline(0,i,src='z') for i in range(len(aa.ants))]
#antpos = numpy.array(antpos)*aipy.const.len_ns/100. #meters
#import IPython; IPython.embed()
#print reds

d_calpar = pickle.load(open(CALPAR,'rb')) # the firstcal initial gains
#XXX Eventually, it would be nice if calpar had gains with keys 'x' or 'y', instead of 'xx'/'yy'
gains = {} #dictionary indexed by pol
for k in d_calpar.keys():
    k_1 = k[0] #XXX just taking first letter (works on 'xx' or 'yy' only)
    gains[k_1] = {} 
    for i in xrange(d_calpar[k].shape[1]):
        #gains[k_1][i] = d_calpar[k][:,i]
        gains[k_1][i] = d_calpar[k][:,i] / numpy.abs(d_calpar[k][:,i])

### Loop Through Compressed Files ###
for f,filename in enumerate(args):
    pol1,pol2 = pol = filename.split('.')[-2] # XXX assumes 1 pol per file
    tag = opts.omnipath + 'zen.'+ '.'.join(filename.split('.')[1:-1])
    print 'Reading '+str(filename)
    t,d,f = capo.arp.get_dict_of_uv_data([filename],antstr='cross',polstr=pol)
    SH = d.values()[0].values()[0].shape
    data,wgts = {}, {}
    g0 = {}
    for i in gains[pol1]: g0[i] = numpy.resize(gains[pol1][i], SH)
    for i,j in reduce(lambda x,y: x+y, reds):
        if i > j: i,j = j,i #always puts lower number first
        bl = aipy.miriad.ij2bl(i,j)
        data[(i,j)] = d[bl][pol]
        wgts[(j,i)] = wgts[(i,j)] = numpy.logical_not(f[bl][pol]).astype(numpy.int) # in case res is switched
    print '   logcal-ing' 
    m1,g1,v1 = omnical.calib.redcal(data,info,gains=g0) #logcal
    print '   lincal-ing'
    #xtalk = compute_xtalk(m1['res'], wgts)
    #xtalk = resize_dict(xtalk, m1['res'].values()[0].shape)
    #xtalk = resize_dict(xtalk, wgts)
    kwargs = {'gains':g1, 'vis':v1,'uselogcal':False,'removedegen':True}
    #if opts.xtalk: kwargs['xtalk'] = xtalk # XXX using xtalk from logcal is BAD
    m2,g2,v2 = omnical.calib.redcal(data,info,**kwargs)
    xtalk = compute_xtalk(m2['res'], wgts)

    # XXX think about commenting of next statment
    #m2['chisq'][data_with_flags==0] = 0 #flag chisq
    
    ### Save Outputs ###
    out = tag + '.omni_output'+opts.omniruntag
    print '   saving '+out+'.npz'
    d_npz = {}
    for bl in xtalk.keys(): #save xtalk
        d_npz['%s,%s,%s' % (pol,'xtalk',bl)] = xtalk[bl] # XXX conj to miriad convention
    d_npz['%s,%s' % (pol,'chisq')] = m2['chisq'] #save chisq
    for vv in v2.keys(): #save vis models
        d_npz['%s,%s,%s' % (pol,'v_log',vv)] = v1[vv]
        d_npz['%s,%s,%s' % (pol,'v_lin',vv)] = v2[vv]
    for aa in g2.keys(): #save antenna gains
        g_pre = gains[pol1][aa] #XXX taking first letter of pol
        g_pre.shape = (1,g_pre.size)
        d_npz['%s,%s,%d' % (pol1,'gains',aa)] = g2[aa] # XXX conjugate to miriad convention
    numpy.savez(out,**d_npz) 
    #import IPython; IPython.embed()
    
