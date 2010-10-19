#!/usr/bin/env python
#
#  plot_trace.py
#  
#
#  Created by Danny Jacobs on 4/5/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m,ephem as ep
import sys, optparse,pickle,tabular as tab,warnings,logging,os
from pylab import *

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True,src=True)
o.set_usage('plot_trace.py <exp_name> <exp_number> file.uv')
#o.add_option('--snap', dest='snap', action='store_true',
#    help='Snapshot mode.  Fits parameters separately for each integration.')
o.add_option('--exclude_parms',dest='exc',
    help='comma list of parameters not to plot')
o.add_option('--prms',
    help='Parameters to plot. default=all')
o.add_option('--plot_fit',action='store_true',
    help='Plot the traces of fits.')
o.add_option('-v',dest='verb',action='store_true',
    help="Print stuff.")
o.add_option('--vv',dest='vverb',action='store_true',
    help="Print even more")    
    
    
opts, args = o.parse_args(sys.argv[1:])
if not opts.exc is None: exc = opts.exc.split(',')
else: exc = []
if opts.vverb: logging.disable(logging.DEBUG)
elif opts.verb: logging.disable(logging.INFO)
else: 
    logging.disable(logging.ERROR)
    warnings.simplefilter('ignore',Warning)
fname = sys.argv[0].split(os.sep)[-1]
log = logging.getLogger(fname)
def get_param_trace(root,ttype):
    """opens a a.fit.*param trace data set. consisting of root.pkl and root.txt
    then proceeds to mash it into a trace dict"""
    prm_list = n.genfromtxt(root+'_prms.txt')
    prm_keys = pickle.load(open(root+'_keys.pkl'))
    prm_trace = a.fit.reconstruct_prms(prm_list[0],prm_keys)
    if not root.endswith('mdl'):
        src_flags = n.genfromtxt('/'.join(root.split('/')[:-1])+'/srcs.txt',
            dtype=None,names=True,deletechars=())
    #mash variables into arrays, assume its a 2 level param dict from a catalog.
    # (Why would I have different?)
    t = n.genfromtxt(root+'_times.txt',dtype={'names':['jd','lst'],'formats':[n.float,n.float]})
    times = None
    if root.endswith('fit'):
        print prm_trace
        scores = n.loadtxt(root+'_scores.txt')
        prm_list = prm_list[scores[:,0]>0]
        t = t[scores[:,0]>0]
    for src,d in prm_trace.iteritems():
        if not root.endswith('mdl'):mask = n.logical_not(src_flags[src])
        for k in d.keys():
            if root.endswith('mdl'):
            #our bm/fit program currently misses the last time so chop those from the model too
#                print 'model hack',root
                prm_trace[src][k] = dict(zip(n.round(t[ttype],5),
                    n.array(prm_list[:-1,prm_keys[src][k][0]])))
                if times is None: times = n.round(t[ttype][:-1],5)
            else:
                prm_trace[src][k] = dict(zip(n.round(t[ttype],5),
                    n.ma.array(prm_list[:,prm_keys[src][k][0]],mask=mask)))
                if times is None: times = n.round(t[ttype],5)
#                if k=='jys' and src=='J1143+221' and root.endswith('fit'):
#                     print prm_trace[src][k],len(n.round(t[ttype],5)),len(prm_list[:,prm_keys[src][k][0]])
    return times,prm_trace
def interp_NaNs(Y):
    yleft = 0
    curnans = []
    for i,y in enumerate(Y):
        if n.isnan(y): curnans.append(i)
        else: 
            if len(curnans)>0 and not 0 in curnans:
                Y[curnans] = (yleft+y)/2.
            elif len(curnans)>0 and 0 in curnans:
                Y[curnans] = y
            curnans = []
            yleft = y
        if n.isnan(y) and i+1==len(Y):
            Y[curnans] = yleft
    return Y


if not opts.src is None:        
    srclist,cutoff,catalogs = a.scripting.parse_srcs(opts.src, opts.cat)
    #load the catalog
    if opts.cal != None:
        cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs)
    else:
        cat = a.src.get_catalog(srclist, cutoff, catalogs)        
            
#load the files to be plotted. Assume we're still dealing with output of beamform_mpi
altnames = {'jys':'S_nu','mfreq':'nu','index':'a'}
assert(len(args)>1)
root = args[0]
fit_times,fit_params = get_param_trace(root+'/fit','lst')
fit_scores = loadtxt(root+'/fit_scores.txt')

bm_times,bm_params = get_param_trace(root + '/bm','lst')

mdl_times,mdl_params = get_param_trace(args[1]+'/mdl','lst')

t = n.sort(n.array(list(set(mdl_times).intersection(set(bm_times)).intersection(set(fit_times)))))
srcs = list(set(mdl_params.keys()).intersection(set(bm_params.keys())).intersection(set(fit_params.keys())))
#print '!',srcs
#print '!!!',[bm_params['J1145+313']['dra'][ti] for ti in t]
#print [mdl_params['J1145+313']['dra'][ti] for ti in t]

if not opts.prms is None: parms = opts.prms.split(',')
else: parms = mdl_params[mdl_params.keys()[0]].keys()
nparms = len(parms)-len(exc)
m1 = n.ceil(n.sqrt(nparms))
m2 = n.int(n.ceil(nparms/m1))
print m1,m2,nparms
#print [l for l in parms if not l in exc]
#print mdl_params.keys(),bm_params.keys()
for i,parm in enumerate([l for l in parms if not l in exc]):
    if parm in exc:continue
    suptitle(root)
    subplot(m1,m2,i+1)
    handles = []
    labels = []
    for name in srcs:
        log.debug("plotting %s for src %s"%(parm,name))
        try:
            bm = n.ma.array([bm_params[name][parm][ti] for ti in t])
            mdl = n.ma.array([mdl_params[name][parm][ti] for ti in t])
        except(KeyError):
            log.info(parm+"not found in bm or mdl")
        if parm=='jys':
            try:
                log.debug("flux is scaled on time series")
                F_mdl = n.ma.array([mdl_params[name]['mfreq'][ti] for ti in t])
                F_bm = n.ma.array([bm_params[name]['mfreq'][ti] for ti in t])
                I_bm = n.ma.average(interp_NaNs(n.ma.array([bm_params[name]['index'][ti] for ti in t])))
                labels.append(name+'bmscale')
                handles.append(
                    plot(t,(bm*(F_mdl/F_bm)**I_bm - mdl)/mdl,'d')
                    )
            except(KeyError):log.info(parm+"not found in beamform")
        elif parm=='index':
            try: 
                if 'jys' in parms:
                    plot(t,(bm - mdl)/mdl,'.',label=name)
                else:
                    handles.append(plot(t,(bm - mdl)/mdl,'.',label=name))
                    labels.append(name)
            except(KeyError): log.info(parm+"not found in mdl and bm")
        if opts.plot_fit:
            try: 
#                fit = [fit_scores[:,0]>0]
                fit = n.ma.array([fit_params[name][parm][ti] for ti in t])
                handles.append(plot(t,(fit-mdl)/mdl,'+',label=name))
                labels.append(name)
            except(KeyError): log.info(parm+"not found in fit")

    if len(handles)>0:
        figlegend(handles,labels,'upper right',numpoints=1)
    title(parm)
#figure()
#print "plotting histograms"
#for i,parm in enumerate([l for l in parms if not l in exc]):
#    if parm in exc:continue
#    suptitle(root)
#    subplot(m1,m2,i+1)
#    for src in mdl_params.keys():
#        print "plotting hist %s for src %s"%(parm,src)
#        try:
#            bm = interp_NaNs(n.array([bm_params[src][parm][ti] for ti in t]))
#            mdl = interp_NaNs(n.array([mdl_params[src][parm][ti] for ti in t]))
#        except(KeyError):
#            print parm, "not found in bm or mdl"        
#        if parm=='jys' and not opts.plot_fit:
#            try:
#                print "flux is scaled on oplot histograms"
#                F_mdl = n.array([mdl_params[src]['mfreq'][ti] for ti in t])
#                F_bm = n.array([bm_params[src]['mfreq'][ti] for ti in t])
#                I_bm = n.average(interp_NaNs(n.array([bm_params[src]['index'][ti] for ti in t])))
#                hist((bm*(F_mdl/F_bm)**I_bm - mdl)/mdl,alpha=0.5,label=src+'bmscale')
#            except(KeyError): print parm,"not found in beamform"
#        elif not opts.plot_fit:
#            try: hist((bm - mdl)/mdl,alpha=0.5,label=src)
#            except(KeyError): print parm,"not found in mdl and bm"
#        if opts.plot_fit:
#            try: 
#                fit = n.array([fit_params[src][parm][ti] for ti in t])
#                hist((fit-mdl)/mdl,label=src,alpha=0.5)
#            except(KeyError): print parm,"not found in fit"
#    figlegend()
#    title(parm+" histograms")
#    xlim([-0.3,0.3])

#plot source parm histograms in their own panes
srcl = n.argsort([s[-4:] for s in srcs])
srcs = [srcs[i] for i in srcl]
m1  = n.ceil(n.sqrt(len(srcs)))
m2  = n.int(n.ceil(len(srcs)/m1))
nsub = m1*m2
print len(srcs),m1,m2
print log.debug("plotting individual source histograms")
print '#'
if not opts.src is None: print '# '+'\t'.join(('Name','Ra','Dec',
                                            'e_S_nu','v_e_S_nu',
                                            'f_e_S_nu','v_f_e_S_nu'))
else: print '# '+'\t'.join(('Name','e_S_nu','v_e_S_nu',
                                    'f_e_S_nu','v_f_e_S_nu'))
for parm in parms:
    figure()
    for i,name in enumerate(srcs):
        if not opts.src is None: 
            src = cat[name]
            ep.FixedBody.compute(src,ep.J2000)
        log.debug("plotting"+name)
        ax = subplot(m2,m1,i+1)
        try:
            bm = interp_NaNs(n.ma.array([bm_params[name][parm][ti] for ti in t]))
            mdl = interp_NaNs(n.ma.array([mdl_params[name][parm][ti] for ti in t]))
            #print bm,mdl
        except(KeyError):
            log.info("not found in bm or mdl")
        #print '+' 
        if parm=='jys' and not opts.plot_fit:
            try:
                log.debug("flux is scaled on subplot histograms")
                F_mdl = n.ma.array([mdl_params[name]['mfreq'][ti] for ti in t])
                F_bm = n.ma.array([bm_params[name]['mfreq'][ti] for ti in t])
                I_bm = n.ma.average(interp_NaNs(n.ma.array([bm_params[name]['index'][ti] for ti in t])))
                E =  (bm*(F_mdl/F_bm)**I_bm - mdl)/mdl
                hist(E.compressed(),alpha=0.5,
                    label=name+'bmscale',normed=1,histtype='step')
                print name,
                if not opts.src is None:
                    print '\t'.join(map(str,(repr(rc.ra),repr(src.dec))))+'\t',
                print '\t'.join(map(str,(n.ma.average((bm*(F_mdl/F_bm)**I_bm - mdl)),
                                n.ma.std((bm*(F_mdl/F_bm)**I_bm - mdl)),
                                n.ma.average((bm*(F_mdl/F_bm)**I_bm - mdl)/mdl),
                                n.ma.std((bm*(F_mdl/F_bm)**I_bm - mdl)/mdl))))
#                if name=='J1421+414':
#                    test =  (bm*(F_mdl/F_bm)**I_bm - mdl)/mdl
#                    print histogram(test.compressed())
                    
            except(KeyError): log.info(parm+"not found in beamform")
        elif not opts.plot_fit:
#            print "plotting",parm
            E = (bm - mdl)/mdl
            try: hist(E.compressed(),alpha=0.5,label=name,normed=1,histtype='step')
            except(KeyError): print log.info(parm+"not found in mdl and bm")
        if opts.plot_fit:
            try: 
                fit = n.ma.array([fit_params[name][parm][ti] for ti in t])
                E = (fit-mdl)/mdl
                hist(E.compressed(),label=name,alpha=0.5,normed=1,histtype='step')
            except(KeyError): log.info(parm+"not found in fit")
        text(0.2,0.85,name,transform = ax.transAxes)
        xlim([-0.3,0.3])
        xticks([-0.2,0.2])
        suptitle(parm+' '+root)
        if (i % m1): yticks([])
        if ((m2-1)*m1>(i%nsub+(m2*m1-nsub))): xticks([])
        #print '.'
    subplots_adjust(hspace=0,wspace=0)
        

        
show()
        


sys.exit()


#
#
#
#
#
#def find_ind_in_keys(id,key_prms):
#    for k,v in key_prms.iteritems():
#        if type(v) is dict:
#            parent = k
#            nparent, name = find_ind_in_keys(id,v)
# #           print nparent,name
#            if not nparent is None and not parent is None and not name is None: 
#                parent = nparent;print '.'
#            if not name is None:
##                print '!',parent,name
#                return parent, name
#        else:
#            #for k,v in key_prms.iteritems():
#  #          print v
#            if v[0]==id:return None, k;
#    return None,None
##load the fit parameters
#keyfile = args[0]+'_fit_'+args[1]+'.txt.kys' 
#fitfile = args[0]+'_fit_'+args[1]+'.txt' 
#
#keys = pickle.load(open(keyfile,'r'))
#find_ind_in_keys(1,keys)
#prms = n.loadtxt(fitfile)
#ptyp = a.fit.reconstruct_prms(prms[0][2:-3],keys)
#finals = n.where(prms[:,-3])[0]
#prm_list = prms[finals,2:-3]
#fit_times = n.array([tuple(l) for l in prms[finals,:2]],dtype={'names':['jd','lst'],'formats':[n.float,n.float]})
#scores = prms[finals,-3:]
#
##load the beam fluxes
#bmfile = args[0]+'_bm_'+args[1]+'.txt'
#bms = n.genfromtxt(bmfile,comments='#',names=True,dtype=None)
#
##load the model fluxes
#try: mdl = tab.tabarray(SVfile=args[2])
#except(IndexError): mdl = None
##find the parameters (assume gsm mode, same prms for all srcs)
#prms keys[keys.keys()[0]].values()
#parm_count = len(keys[keys.keys()[0]].values())
#
#nvalid = n.sum([1 for i in range(parm_count) if not find_ind_in_keys(i,keys)[1] is None]) #worst line of code ever!
#srcs = keys.keys()#cheating!!! to get my source list. This will not make sense tomorrow. :(
#m1 = n.int(n.round(n.sqrt(nvalid)))
#m2 = n.int(n.ceil(nvalid/m1))
#print m2,m1
#print "TIME"
#print len(mdl['lst']),len(fit_times)
##for t1,t2,t3 in zip(fit_times['lst'],bms['lst'],mdl['lst']):
##    print t1,t2,t3
#for src in srcs:
#    for i in range(prm_list[):
#        parent, name  = find_ind_in_keys(i,keys)
#        print parent,src
#        if parent!=src: continue
#        if not name is None: 
#            print parent,
#            subplot(m1,m2,i+1)
#            plot(fit_times['lst'],prm_list[:,i],'.'label=)
#            print src,prm_list[0,i]
#            print name,mdl.dtype.names,
#            #switch to vot names
#            if name in altnames.keys():
#                aname = altnames[name]
#            else: aname = name
#            print aname
#            if not mdl is None and aname in mdl.dtype.names:
#                 plot(mdl['lst'],mdl[aname],'+')
#            if aname is 'S_nu':
#                cbms = n.where(bms['src']==parent)[0]
#                print cbms,bms[aname][cbms]
#                plot(bms['lst'][cbms],bms[aname][cbms],'x')
#            title(name)
##        print fit_times[:,1],prm_list[:,i]
#
#show()
#
#
#
#
#sys.exit()
#file = args[0]
#bmnames = open(file+'bm').readlines()[0].split()
#cat = a.src.get_catalog(srcs=bmnames[2:],catalogs=['misc','three_cr'])
#bm = n.loadtxt(file+'bm')
#fitnames = open(file).readlines()[0].split()
#fit = n.loadtxt(file)
#
#mdllines = open(args[1]).readlines()
#mdl = {} 
#for l in mdllines:
#    if not l.startswith('Time'):
#        if mdl.has_key(l.split()[0]): mdl[l.split()[0]].append(l.split()[1])
#        else: mdl[l.split()[0]] = [l.split()[1]]
#
#figure()
#subplot(2,1,1)
#for i,b in enumerate(bmnames[2:]):
#    i += 2
#    plot(bm[:,1],bm[:,i],'.',label=b+'bm')
#for k,v in mdl.iteritems():
#    plot(bm[:,1],v[:-1],'+',label=k)
#title('Beamform.')
#legend()
#
#subplot(2,1,2)
#for i,b in enumerate(fitnames[2:-2]):
#    i += 2
#    print (130. /cat[b.split('_')[0]].mfreq)**cat[b.split('_')[0]].index
#    plot(fit[:,1],fit[:,i]*(0.130 /cat[b.split('_')[0]].mfreq)**cat[b.split('_')[0]].index,'.',label=b)
#for k,v in mdl.iteritems():
#    plot(bm[:,1],v[:-1],'+',label=k)
#    mdl[k] = n.array(v[:-1]).astype(n.float)
#title('fit')
#legend()    
#
#figure()
##plot a histogram of error for each source
#subplot(2,1,1)
#for i,s in enumerate(bmnames):
#    if mdl.has_key(s):
#        print bm[0,i]
#        hist((mdl[s]-bm[:,i])/mdl[s],bins=20,label=s,alpha=0.4)
#legend()
#xlim([-0.03,0.03])
#title("beamform")
#subplot(2,1,2)
#for i,s in enumerate(fitnames):
#    try:
#        name = s.split('_')[0]
#        parm = s.split('_')[1]
#        if mdl.has_key(name) and parm=='jys':
#            hist((mdl[name]-fit[:,i]*(0.130 /cat[name].mfreq)**cat[name].index)/mdl[name],bins=20,label=name,alpha=0.4)
#    except(IndexError):pass
#legend()
#title("fit")
#xlim([-0.03,0.03])
#show()
#
