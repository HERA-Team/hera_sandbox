#! /usr/bin/env python
import sys, numpy as np, aipy, optparse, capo
from matplotlib import pyplot as plt

o = optparse.OptionParser()
o.set_usage('angleVar.py [options] *.uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('--ba', dest='badants', default='', help='comma separated list of bad antennae')
opts,args = o.parse_args(sys.argv[1:])

#get antenna positions
print 'reading, %s'%opts.cal
exec("import {calfile} as cal".format(calfile=opts.cal))
antpos = cal.prms['antpos']

#get initial info from first file
uv = aipy.miriad.UV(args[0])
nants,nchan = uv['nants'],uv['nchan']
uv.select('antennae',0,1,include=True) #XXX assumes ants 0 and 1 are in uv file
_T=[]
for p,d in uv.all(): 
    _,_t,_ = p
    _T.append(_t)
ntimes = len(_T)
del(uv)
assert(nants == len(antpos.keys())) #check the cal file and uv files are going to cooperate

var_stor = np.zeros((nants),dtype='complex128')
flg_stor = np.zeros_like(var_stor)

for uv in args:
    print '    Reading %s...'%uv
    times,data,flags = capo.arp.get_dict_of_uv_data([uv],'cross',opts.pol) 
    for i in range(nants):
        for j in range(nants):
            if i==j: continue #neglect autos
            try: 
                var_stor[i]+=np.var(np.angle(data[(i,j)][opts.pol].T)).real #sigma_i1^2 + sigma_i2^2 + ...
                flg_stor[i] += 1. #N
            except KeyError:
                if i<j: print 'KeyError on (%i,%i)'%(i,j) #this should not happen
                continue

sigma_stor = np.sqrt(var_stor/flg_stor)

#parse options
if not len(opts.badants)==0: badants=map(int,opts.badants.split(','))
else: badants=[]

#XXX There's an array manipulation way to do this
nobad_sigma_stor = []
for i in range(nants):
    if i not in badants:
        nobad_sigma_stor.append(sigma_stor[i])
good_avg = np.nanmean(nobad_sigma_stor)
good_std = np.nanstd(nobad_sigma_stor)

#Plot
plt.plot(range(nants),sigma_stor,'bo')

out = []
for i in range(nants):
    if sigma_stor[i]>=good_avg+good_std and not i in badants:
        out.append(i)
    if i in badants:
        plt.plot(i,sigma_stor[i],'kx',ms=10)
plt.grid()
plt.axhline(good_avg,color='k')
plt.fill_between(range(nants),good_avg,good_avg+good_std,alpha=0.3)
plt.xlim(-0.5,112)
plt.show()
print out

