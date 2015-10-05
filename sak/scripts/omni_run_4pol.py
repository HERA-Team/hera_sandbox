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
aipy.scripting.add_standard_options(o, pol=True, cal=True, ant=True)
o.add_option('--cpx',dest='calparXX',type='string',
            help='Path and name of XX .p calpar file.')
o.add_option('--cpy',dest='calparYY',type='string',
            help='Path and name of YY .p calpar file.')            
o.add_option('--redinfo',dest='redinfo',type='string',
            help='Path and name of .bin redundant info file.')
o.add_option('--xtalk',dest='xtalk',default=False,action="store_true",
            help='Option to use xtalk command when performing lincal. Default is False.')
o.add_option('--omniruntag',dest='omniruntag',default='',type='string',
            help='Tag for omni run, if wanted. Default is empty.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .omni_output npz files. Include final "/" in path.')
opts,args = o.parse_args(sys.argv[1:])

### Parse polarization cases ###

if opts.pol=='xy' or opts.pol=='yx': assert(opts.calparXX!=None and opts.calparYY!=None)
if opts.pol=='xx': assert(opts.calparYY==None)
if opts.pol=='yy': assert(opts.calparXX==None)

### Save Options ###
CALFILE = opts.cal
pol=opts.pol
CALPAR_XX = opts.calparXX
CALPAR_YY = opts.calparYY
RED_INFO = opts.redinfo

### Read Firstcal Info ###

print 'Reading calpar and redinfo files...'

info = omnical.info.RedundantInfoLegacy() #reading old firstcal files
info.fromfile(RED_INFO)
reds = info.get_reds()

d_calpar_XX = pickle.load(open(CALPAR_XX,'rb'))
d_calpar_YY = pickle.load(open(CALPAR_YY,'rb'))

gains = {}

gains['x'], gains['y'] = {}, {} 
for i in xrange(d_calpar_XX['xx'].shape[1]): gains['x'][i] = d_calpar_XX['xx'][:,i]
for i in xrange(d_calpar_YY['yy'].shape[1]): gains['y'][i] = d_calpar_YY['yy'][:,i]

"""
### Construct antenna array
uv = aipy.miriad.UV(args[0])
(j,t,j),j = uv.read()
aipy.scripting.uv_selector(uv, opts.ant, opts.pol)
aa = aipy.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
del(uv)
"""
### Loop Through Files ###

for f,file in enumerate(args):
    pol1,pol2 = opts.pol[0],opts.pol[1]
    tag = opts.omnipath + 'zen.'+ '.'.join(file.split('.')[1:-1])
    print str(f+1)+'/'+str(len(args))+': '+'Reading '+str(file)

    t,d,f = capo.arp.get_dict_of_uv_data([file],antstr='cross',polstr=opts.pol)
    data = {}
    reds_all = reduce(lambda x,y: x+y,reds)
    for r,red in enumerate(reds_all):
        i,j = red #red is (64,49), for example
        if i>j: i,j = j,i #always puts lower number first
        try: g_ij = gains[pol1][i]*gains[pol2][j].conj()
        except: print 'ERROR WITH GAINS:\n gains[pol1][i]*gains[pol2][j].conj()'
        g_ij = g_ij.conj() #Omnical conjugation convention is backwards
        g_ij.shape = (1,g_ij.size)
        #print i,j,g_ij,pol
        data[(i,j)] = d[aipy.miriad.ij2bl(i,j)][pol]/g_ij #gains and data always have lower number first 
        if r == 0:
            data_with_flags = data[(i,j)]
        #XXX data should be 0 when g_ij is 0? Right now it becomes nan's
    
    ### LOGCAL
    
    print '   logcal-ing' 
    m,g,v = omnical.calib.redcal(data,info) #logcal
    
    
    print '   lincal-ing'
    
    ### Parsing xtalk (Carina TODO)
    if opts.xtalk == True:
        xtalk = {}
        xtalk_flat = {}
        for key in m['res'].keys():
            xtalkavg = numpy.mean(m['res'][key],axis=0) #avg over time
            xtalk_flat[key] = xtalkavg #saved for later (not reshaped to minimize size)
            xtalk[key] = numpy.resize(xtalkavg,(m['res'][key].shape)) #must be same shape as data
    else: xtalk = None
    
    
    ### LINCAL
    
    m2,g2,v2 = omnical.calib.redcal(data, info, xtalk=xtalk, gains=g, vis=v, uselogcal=False, removedegen=True) #lincal
    m2['chisq'][data_with_flags==0] = 0 #flag chisq
    
    
    ### Save Outputs ###
    out = tag + '.omni_output'+opts.omniruntag
    print '   saving '+out+'.npz'
    d_npz = {}
    if opts.xtalk == True:
        for bl in xtalk_flat.keys(): #save xtalk
            d_npz['%s,%s,%s' % (pol,'xtalk',bl)] = xtalk_flat[bl]
    d_npz['%s,%s' % (pol,'chisq')] = m2['chisq'] #save chisq
    #for vv in v.keys(): #save vis models
        #d_npz['%s,%s,%s' % (pol,'v_log',vv)] = v[vv]
        #d_npz['%s,%s,%s' % (pol,'v_lin',vv)] = v2[vv]
    for A in g.keys(): #save antenna gains
        g_pre = gains[pol1][A] #XXX taking first letter of pol
        g_pre.shape = (1,g_pre.size)
        d_npz['%s,%s,%d' % (pol1,'gains',A)] = g_pre*g2[A]
    numpy.savez(out,**d_npz)
    
