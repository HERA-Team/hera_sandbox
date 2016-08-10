#! /usr/bin/env python

import omnical
import aipy
import numpy
import capo
import pickle
import optparse
import os, sys

### Options ###
o = optparse.OptionParser()
o.set_usage('run_omnicalv2.py [options] *uvcRRE')
o.set_description(__doc__)
o.add_option('--calpar',dest='calpar',type='string',
            help='Path to .p calpar file.')
o.add_option('--redinfo',dest='redinfo',type='string',
            help='Path to .bin redundant info file.')
o.add_option('-C',dest='cal',default='psa6622_v003',type='string',
            help='Path to calfile.')
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
for k in d_calpar.keys():
    gains[k] = {} 
    for i in xrange(d_calpar[k].shape[1]):
        gains[k][i] = d_calpar[k][:,i]

for file in range(len(args)):
    
    print 'Reading',args[file]

    t,d,f = capo.arp.get_dict_of_uv_data([file],antstr='cross',polstr='xx')
    
    data = {}
    reds_all = reduce(lambda x,y: x+y,reds)
    
    for red in reds_all:
        i,j = red #red is (64,49), for example
        if i>j: i,j = j,i #always puts lower number first
        g_ij = gains['xx'][i]*gains['xx'][j].conj()
        g_ij = g_ij.conj() ### Omnical conjugation convention is backwards ###
        g_ij.shape = (1,g_ij.size)
        data[(i,j)] = d[aipy.miriad.ij2bl(i,j)]['xx']/g_ij #gains and data always have lower number first 
                ### data should be 0 when g_ij is 0 ###
    
    m,g,v = omnical.calib.redcal(data,info) #logcal
    m2,g2,v2 = omnical.calib.redcal(data,info,xtalk=m['res'],gains=g,vis=v,uselogcal=False,removedegen=True) #lincal
    
    #write stuff into npz file with keys (diff pols)
        #multiple first_cal gains to g and g2
    #apply calibration on a per file and then per time basis in diff script
    
    import IPython; IPython.embed()
    
