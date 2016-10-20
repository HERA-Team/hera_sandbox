#!/usr/bin/env python
import aipy, capo as C, optparse, numpy as np, matplotlib.pyplot as plt, sys
o = optparse.OptionParser()
o.set_description(__doc__)
aipy.scripting.add_standard_options(o, cal=True, ant=True, pol=True)
o.add_option('-b','--badants',dest='ba',default=None,help='bad antennae to remove, separated by commas. e.g. "-b 1,2,3"')
opts, args = o.parse_args(sys.argv[1:])

#Array data
ant = int(opts.ant)
if not opts.ba is None: badants = map(int,opts.ba.split(','))
else: badants = []
print 'reading, %s'%opts.cal
exec("import {calfile} as cal".format(calfile=opts.cal))
antpos = cal.prms['antpos']
nants = len(antpos.keys())
x0,y0 = antpos[ant]['top_x'],antpos[ant]['top_y']

#get data
t,d,f = C.arp.get_dict_of_uv_data(args,antstr=opts.ant,polstr=opts.pol)

bl_length,averages,stdevs = [],[],[]
#plot, skipping bad antennae
for a in [x for x in xrange(nants) if x not in badants]:
    if ant<a: tup = (ant,a)
    elif ant>a: tup = (a,ant)
    else: continue
    D = np.absolute(d[tup][opts.pol])
    averages.append(np.nanmean(D))
    stdevs.append(np.nanstd(D))
    x1,y1 = antpos[a]['top_x'],antpos[a]['top_y']
    dx,dy = np.sqrt(np.power(x0-x1,2.)),np.sqrt(np.power(y0-y1,2.))
    L = np.sqrt(np.power(dx,2.) + np.power(dy,2.))
    bl_length.append(L)
plt.errorbar(bl_length,averages,yerr=stdevs,fmt='o',ecolor='b')
plt.xlabel('Basline length [m]')
plt.ylabel(r'$\langle | V_{a,j} | \rangle_{t,\nu}$',size=20)
plt.yscale('log', nonposy='clip')
plt.show()
    
    
    
