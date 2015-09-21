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
o.add_option('--redinfo',dest='redinfo',type='string',
            help='Path and name of .bin redundant info file.')
o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
            help='Path and name of calfile.')
o.add_option('--omniruntag',dest='omniruntag',default='',type='string',
            help='Tag for omni run, if wanted. Default is empty.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .omni_output npz files. Include final / in path.')
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
CALPAR = opts.calpar
RED_INFO = opts.redinfo

### Read Firstcal Info ###

print 'Reading calpar and redinfo files...'

info = omnical.info.RedundantInfoLegacy() #reading old firstcal files
info.fromfile(RED_INFO)
reds = info.get_reds()
#reds = omnical.arrayinfo.filter_reds(reds,ubls=[(64,49)])
#info.init_from_reds(reds,info.get_antpos())
### need to test above ###

#aa = aipy.cal.get_aa(CALFILE,numpy.array([.15]))
#antpos = [a.pos for a in aa] #nanoseconds
#antpos = numpy.array(antpos) #antpos needs to be (#ant, 3) array with positions
#antpos = antpos*aipy.const.len_ns/100. #meters
#antpos = [aa.get_baseline(0,i,src='z') for i in range(len(aa.ants))]
#antpos = numpy.array(antpos)*aipy.const.len_ns/100. #meters
#reds = omnical.arrayinfo.compute_reds(antpos,tol=0.1) #bypasses creating AI class
#RI = omnical.info.RedundantInfo() #make RI class
#info.init_from_reds(reds,antpos)

d_calpar = pickle.load(open(CALPAR,'rb'))
#XXX Eventually, it would be nice if calpar had gains with keys 'x' or 'y', instead of 'xx'/'yy'
gains = {} #dictionary indexed by pol
for k in d_calpar.keys():
    k_1 = k[0] #XXX just taking first letter (works on 'xx' or 'yy' only)
    gains[k_1] = {} 
    for i in xrange(d_calpar[k].shape[1]):
        gains[k_1][i] = d_calpar[k][:,i]

### Loop Through Compressed Files ###

for f in range(len(args)):
    
    file = args[f]
    pol = file.split('.')[-2]
    pol1,pol2 = pol[0],pol[1]
    tag = opts.omnipath + 'zen.'+ '.'.join(file.split('.')[1:-1])
    print str(f+1)+'/'+str(len(args))+': '+'Reading '+str(file)

    t,d,f = capo.arp.get_dict_of_uv_data([file],antstr='cross',polstr=pol)
    data = {}
    reds_all = reduce(lambda x,y: x+y,reds)
    for r,red in enumerate(reds_all):
        i,j = red #red is (64,49), for example
        if i>j: i,j = j,i #always puts lower number first
        try: g_ij = gains[pol1][i]*gains[pol2][j].conj()
        except: print '!!! Pol of data file is not in calpar file.'
        g_ij = g_ij.conj() #Omnical conjugation convention is backwards
        g_ij.shape = (1,g_ij.size)
        data[(i,j)] = d[aipy.miriad.ij2bl(i,j)][pol]/g_ij #gains and data always have lower number first 
        if r == 0: #get flags from one file
            data_with_flags = data[(i,j)]
    print '   logcal-ing' 
    m,g,v = omnical.calib.redcal(data,info) #logcal
    print '   lincal-ing'
    xtalk = {}
    xtalk_flat = {}
    for key in m['res'].keys():
        xtalkavg = numpy.mean(m['res'][key],axis=0) #avg over time
        xtalk_flat[key] = xtalkavg #saved for later (not reshaped to minimize size)
        xtalk[key] = numpy.resize(xtalkavg,(m['res'][key].shape)) #must be same shape as data
    m2,g2,v2 = omnical.calib.redcal(data,info,xtalk=xtalk,gains=g,vis=v,uselogcal=False,removedegen=True) #lincal
    m2['chisq'][data_with_flags==0] = 0 #flag chisq
    
    ### Save Outputs ###
    
    out = tag + '.omni_output'+opts.omniruntag
    print '   saving '+out+'.npz'
    d_npz = {}
    for bl in xtalk_flat.keys(): #save xtalk
        d_npz['%s,%s,%s' % (pol,'xtalk',bl)] = xtalk_flat[bl]
    d_npz['%s,%s' % (pol,'chisq')] = m2['chisq'] #save chisq
    for vv in v.keys(): #save vis models
        d_npz['%s,%s,%s' % (pol,'v_log',vv)] = v[vv]
        d_npz['%s,%s,%s' % (pol,'v_lin',vv)] = v2[vv]
    for aa in g.keys(): #save antenna gains
        g_pre = gains[pol1][aa] #XXX taking first letter of pol
        g_pre.shape = (1,g_pre.size)
        d_npz['%s,%s,%d' % (pol1,'gains',aa)] = g_pre*g2[aa]
    numpy.savez(out,**d_npz) 
    #import IPython; IPython.embed()
    
