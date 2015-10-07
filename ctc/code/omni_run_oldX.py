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
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .omni_output npz files. Include final / in path.')
o.add_option('--xtalk',dest='xtalk',default=False,action="store_true",
            help='Option to use xtalk command when performing lincal. Default is False.')
o.add_option('--omniruntag',dest='omniruntag',default='',type='string',
            help='Tag for omni run, if wanted. Default is empty.')
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

#aa = aipy.cal.get_aa(CALFILE,numpy.array([.15]))
#antpos = [a.pos for a in aa] #nanoseconds
#antpos = numpy.array(antpos) #antpos needs to be (#ant, 3) array with positions
#antpos = antpos*aipy.const.len_ns/100. #meters
#reds = omnical.arrayinfo.compute_reds(antpos,tol=0.1) #bypasses creating AI class
#reds = AI.filter_reds(reds,ubls=[(64,65)])
#print reds
#RI = omnical.info.RedundantInfo() #make RI class
#RI.init_from_reds(reds,antpos)

d_calpar = pickle.load(open(CALPAR,'rb'))
gains = {} #dictionary indexed by pol
pol = d_calpar.keys()[0] #XXX only set up for 1 pol right now
for k in d_calpar.keys():
    gains[k] = {} 
    for i in xrange(d_calpar[k].shape[1]):
        gains[k][i] = d_calpar[k][:,i]

### Loop Through Compressed Files ###

for f in range(len(args)):
    
    file = args[f]
    #tag = '.'.join(file.split('.')[:-1])
    tag = opts.omnipath + 'zen.'+ '.'.join(file.split('.')[1:-1])
    print str(f+1)+'/'+str(len(args))+': '+'Reading '+str(file)

    t,d,f = capo.arp.get_dict_of_uv_data([file],antstr='cross',polstr=pol)
    data = {}
    reds_all = reduce(lambda x,y: x+y,reds)
    
    for red in reds_all:
        i,j = red #red is (64,49), for example
        if i>j: i,j = j,i #always puts lower number first
        g_ij = gains[pol][i]*gains[pol][j].conj()
        g_ij = g_ij.conj() #Omnical conjugation convention is backwards
        g_ij.shape = (1,g_ij.size)
        data[(i,j)] = d[aipy.miriad.ij2bl(i,j)][pol]/g_ij #gains and data always have lower number first 
        #XXX data should be 0 when g_ij is 0? Right now it becomes nan's
    data_with_flags = data[(reds_all[0][0],reds_all[0][1])] #flags from first file
    
    print '   logcal-ing' 
    m,g,v = omnical.calib.redcal(data,info) #logcal
    xtalk = numpy.mean(m['res'],axis=0) #avg over time
    xtalk = numpy.resize(xtalk,m['res'].shape)
    print '   lincal-ing'
    if opts.xtalk == True:
        m2,g2,v2 = omnical.calib.redcal(data,info,xtalk=xtalk,gains=g,vis=v,uselogcal=False,removedegen=True) #lincal
    else:
        m2,g2,v2 = omnical.calib.redcal(data,info,gains=g,vis=v,uselogcal=False,removedegen=True) #lincal without xtalk
    m2['chisq'][data_with_flags==0] = 0 #flags chisq where data is flagged

    ### Save Outputs ###
    out = tag + '.omni_output'+opts.omniruntag #change this manually for diff runs if needed
    print '   saving '+out+'.npz'
    d_npz = {}
    d_npz[pol] = {}
    for mm in m.keys():
        d_npz['%s,%s,%s' % (pol,'m_log',mm)] = m[mm]
        d_npz['%s,%s,%s' % (pol,'m_lin',mm)] = m2[mm]
    for vv in v.keys():
        d_npz['%s,%s,%s' % (pol,'v_log',vv)] = v[vv]
        d_npz['%s,%s,%s' % (pol,'v_lin',vv)] = v2[vv]
    for aa in g.keys():
        g_pre = gains[pol][aa]
        g_pre.shape = (1,g_pre.size)
        d_npz['%s,%s,%d' % (pol,'g_log',aa)] = g_pre*g[aa]
        d_npz['%s,%s,%d' % (pol,'g_lin',aa)] = g_pre*g2[aa]
    numpy.savez(out,**d_npz) 
    #import IPython; IPython.embed()
    
