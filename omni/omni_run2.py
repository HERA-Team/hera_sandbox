#! /usr/bin/env python

import omnical, aipy, numpy, capo
import pickle, optparse, os, sys

o = optparse.OptionParser()
o.set_usage('omni_run.py [options] *uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o,cal=True,pol=True)
o.add_option('--calpar',dest='calpar',type='string',
            help='Path and name of POL.p ("xx.p") calpar file.')
o.add_option('--redinfo',dest='redinfo',type='string',default='',
            help='Path and name of .bin redundant info file.')
o.add_option('--omnipath',dest='omnipath',default='',type='string',
            help='Path to save .npz files. Include final / in path.')
o.add_option('--ba',dest='ba',default=None,
            help='Antennas to exclude, separated by commas.')
opts,args = o.parse_args(sys.argv[1:])

#Dictionary of calpars and files
pols = opts.pol.split(',')
cps,files = {},{}
for i,p in enumerate(pols):
    files[p] = [] #list of files per pol
    if p in opts.calpar: #dictionary of calpars per pol
        print 'Reading', opts.calpar
        cps[p] = pickle.load(open(opts.calpar,'rb'))[p]
    else: #XXX assumes calpar naming is *pol.p
        new_cp = opts.calpar.split('.p')[0][:-2]+p+'.p'
        print 'Reading', new_cp
        cps[p] = pickle.load(open(new_cp,'rb'))[p]
for filename in args:
    for p in pols:
        fn = filename.split('.')
        fn[3] = p
        files[p].append('.'.join(fn))     


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
    #ex_ants = [5,7,8,15,16,17,24,26,27,28,29,37,38,46,48,50,51,53,55,63,68,69,72,74,76,77,82,83,84,85,92,107,110] #antennas to exclude
    #ex_ants = [72,53] + [48,82,83,37] + [107,15,16] + [92,7] + [27,50,54] + [26,38,110,46]
    #ex_ants += [29,24,28,55] + [68,17,69,13] + [76,5,77,32] + [84,40,85,14] # are all these bad?
    #ex_ants += [57,63,8,74]
    info = capo.omni.aa_to_info(aa, pols=['x'], ex_ants=ex_ants,ex_crosspols=['xy','yx'])

reds = info.get_reds()
reds_02 = omnical.arrayinfo.filter_reds(reds, ubls=[(64,10)]) #30m E/W baselines
#print ','.join(['%d_%d' % (i,j) for i,j in reds_02[0]])

print 'Getting gains from calpar'
#XXX Eventually, it would be nice if calpar had gains with keys 'x' or 'y', instead of 'xx'/'yy'
g0 = {} 
for k in calpar.keys(): #loop over pol
    g0[k] = {}
    #g0['yy'] = {}
    for i in xrange(calpar[k].shape[1]): #loop over antennas
        #g0['yy'][i] = calpar[k][:,i] / numpy.abs(calpar[k][:,i]) #gains have len(nfreq)
        g0[k][i] = numpy.conj(calpar[k][:,i] / numpy.abs(calpar[k][:,i]))
#XXX g0 pol should eventually be just 'x' or 'y' instead of 'xx'/'yy'

### Loop Through Compressed Files ###
for f,filename in enumerate(args):
    print 'Reading', filename
    pol = filename.split('.')[-2] #XXX assumes 1 pol per file
    timeinfo,d,f = capo.arp.get_dict_of_uv_data([filename],antstr='cross',polstr=pol) #XXX could try reading only bls in reds
    t_jd = timeinfo['times']
    t_lst = timeinfo['lsts']
    freqs = numpy.linspace(.1,.2,len(d[d.keys()[0]][pol][0]))
    SH = d.values()[0].values()[0].shape #shape of file data (ex: (19,203))
    data,wgts,xtalk = {}, {}, {}
    m2,g2,v2 = {}, {}, {}
    for i in g0[pol]: g0[pol][i] = numpy.resize(g0[pol][i],SH) #resize gains like data
    g0_test = {}
    g0_test[pol[0]] = g0[pol] #XXX gain pol should only be one letter!!!
    for bl in d: 
        i,j = bl
        #i,j = aipy.miriad.bl2ij(bl)
        data = d
        #data[(i,j)] = d[bl][pol]
        wgts[(j,i)] = wgts[(i,j)] = numpy.logical_not(f[bl][pol]).astype(numpy.int)
    #import IPython;IPython.embed()
    print '   Logcal-ing'
    #import IPython;IPython.embed() 
    m1,g1,v1 = capo.omni.redcal(data,info,gains=g0_test)
    #m1,g1,v1 = omnical.calib.redcal(data,info,gains=g0[pol])
    print '   Lincal-ing'
    m2[pol],g2[pol[0]],v2[pol] = capo.omni.redcal(data, info, gains=g1, vis=v1, uselogcal=False, removedegen=True)
    import IPython;IPython.embed()
    xtalk[pol] = capo.omni.compute_xtalk(m2[pol]['res'], wgts) #xtalk is time-average of residual
    print '   Saving '+opts.omnipath+filename.split('/')[-1]+'2.npz'
    m2[pol]['history'] = 'OMNI_RUN: ' + ' '.join(sys.argv) + '\n'
    capo.omni.to_npz(opts.omnipath+filename.split('/')[-1]+'2.npz', m2, g2, v2, xtalk,t_jd,t_lst,freqs)
    
