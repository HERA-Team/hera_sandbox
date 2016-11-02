#!/usr/bin/env python
import aipy, capo as C, optparse, numpy as np, matplotlib.pyplot as plt, sys
o = optparse.OptionParser()
o.set_description(__doc__)
o.set_usage('find_crosspols.py *pp.uvcRRE (only use one pol, ask for others in option)')
aipy.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('-b','--badants',dest='ba',default=None,help='bad antennae to remove, separated by commas. e.g. "-b 1,2,3"')
o.add_option('-m','--maskbad',dest='mb',action='store_true',help='Instead of plotting an "X" over bad antennae, neglect plotting them entirely.')
o.add_option('--noBL',dest='nobl',action='store_true',help='Turn off plot with BL length')
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
anchor_ant = int(opts.ant)
if not opts.ba is None: badants = map(int,opts.ba.split(','))
else: badants = []
print 'reading, %s'%opts.cal
exec("import {calfile} as cal".format(calfile=opts.cal))
antpos = cal.prms['antpos']
nants = len(antpos.keys())
x0,y0 = antpos[anchor_ant]['top_x'],antpos[anchor_ant]['top_y']

#data read
tp1,dp1,fp1 = C.arp.get_dict_of_uv_data(p1files,antstr=opts.ant,polstr=pol_strs[0])
tp2,dp2,fp2 = C.arp.get_dict_of_uv_data(p2files,antstr=opts.ant,polstr=pol_strs[1])

#data analysis
bl_length,avg_ratios,stdevs = [],[],[]
for ant in range(nants):
    #parse miriad keys
    if anchor_ant == ant: continue #neglect auto
    elif anchor_ant < ant: tup = (anchor_ant,ant)
    else: tup = (ant,anchor_ant)
    try:
        avg_p1 = np.nanmean(np.absolute(dp1[tup][pol_strs[0]]))
    except IndexError:
        print 'Index error on antenna %i'%ant
        continue
    avg_p2 = np.nanmean(np.absolute(dp2[tup][pol_strs[1]]))
    ratio = avg_p2/avg_p1
    std_p1 = np.nanstd(np.absolute(dp1[tup][pol_strs[0]]))
    std_p2 = np.nanstd(np.absolute(dp2[tup][pol_strs[1]]))
    std = ratio*np.sqrt((std_p1/avg_p1)**2. + (std_p2/avg_p2)**2.)
    x1,y1 = antpos[ant]['top_x'],antpos[ant]['top_y']
    dx,dy = np.abs(x0-x1),np.abs(y0-y1)
    L = np.sqrt(dx**2. + dy**2.)
    
    avg_ratios.append(ratio)
    stdevs.append(std)
    bl_length.append(L)
    #indiv plot
    if opts.mb and ant in badants: continue
    plt.errorbar(L,ratio,yerr=std,fmt='o',ecolor='b',color='b')
    plt.text(L,ratio,str(ant))
    if ant in badants: plt.plot(L,ratio,'kx',ms=10)

if opts.verb: print np.where(avg_ratios > np.nanmean(avg_ratios)+2*np.nanstd(avg_ratios))
#format
plt.xlabel('Basline length [m]')
plt.xlim(0,300)
plt.ylabel(r'$\langle | V_{a,j} | \rangle_{t,\nu}$',size=20)
plt.suptitle('Relative to antenna %i'%anchor_ant,size=15)
if not opts.nobl: plt.show()
plt.close()

import IPython;IPython.embed()
