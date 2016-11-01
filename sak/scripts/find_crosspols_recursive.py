#!/usr/bin/env python
import aipy, capo as C, optparse, numpy as np, matplotlib.pyplot as plt, sys
o = optparse.OptionParser()
o.set_description(__doc__)
o.set_usage('find_crosspols_recursive.py *pp.uvcRRE (only use one pol, ask for others in option)')
aipy.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('-A','--startant',dest='startant',default='0',help='antenna to start relative to (default is all antennae, i.e. start with zero)')
o.add_option('-b','--badants',dest='ba',default=None,help='bad antennae to remove, separated by commas. e.g. "-b 1,2,3"')
o.add_option('-v','--verbose',dest='verb',action='store_true',help='Toggle verbosity')
opts, args = o.parse_args(sys.argv[1:])

#parse pols
pol_strs = opts.pol.split(',')
if len(pol_strs)!=2: 
    print 'use one linpol and one crosspol'
    raise Exception

#parse files: assumes both pols are in the same directory
if pol_strs[0] in args[0]:
    p1files = args
    p2files = [w.replace(pol_strs[0],pol_strs[1]) for w in args]
else:
    p2files = args
    p1files = [w.replace(pol_strs[1],pol_strs[0]) for w in args]

#parse array data

if not opts.ba is None: badants = map(int,opts.ba.split(','))
else: badants = []
print 'reading, %s'%opts.cal
exec("import {calfile} as cal".format(calfile=opts.cal))
antpos = cal.prms['antpos']
nants = len(antpos.keys())

all_avgs = np.zeros((nants,nants))
offenders = []
if opts.verb: print 'Rel_ant : avg>2sigma'
for anchor_ant in range(int(opts.startant),nants):
    if anchor_ant in badants:
        all_avgs[anchor_ant,:] = np.nan
        continue
    #data read
    tp1,dp1,fp1 = C.arp.get_dict_of_uv_data(p1files,antstr=str(anchor_ant),polstr=pol_strs[0])
    tp2,dp2,fp2 = C.arp.get_dict_of_uv_data(p2files,antstr=str(anchor_ant),polstr=pol_strs[1])
    #data analysis
    bl_length,avg_ratios,stdevs = [],[],[]
    for ant in range(nants):
        #parse miriad keys
        if anchor_ant == ant:
            #neglect auto
            avg_ratios.append(np.nan)
            continue
        elif anchor_ant < ant: tup = (anchor_ant,ant)
        else: tup = (ant,anchor_ant)
        try:
            avg_p1 = np.nanmean(np.absolute(dp1[tup][pol_strs[0]]))
        except IndexError:
            print 'Index error on antenna %i'%ant
            continue
        avg_p2 = np.nanmean(np.absolute(dp2[tup][pol_strs[1]]))
        ratio = avg_p2/avg_p1
        avg_ratios.append(ratio)
    del(tp1);del(dp1);del(fp1);del(tp2);del(dp2);del(fp2)
    
    #crosspol'd antennae are the ones that deviate from avg by 2 sigma
    crosspols = np.where(avg_ratios > np.nanmean(avg_ratios)+2*np.nanstd(avg_ratios))[0]
    if opts.verb: print anchor_ant,':',crosspols
    offenders.append(crosspols)
    all_avgs[anchor_ant,:] = avg_ratios
