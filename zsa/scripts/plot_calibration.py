#! /usr/bin/env python 
import numpy as n
import pylab as p
import sys, glob, os
import capo, aipy as a

args = sys.argv[1:]
aa = a.cal.get_aa('psa6240_v003', .1,.1,.1)
seps=['0,1','-1,1','1,1']

sep2ij = {}
ij2sep = {}
s = ''
for sep in seps:
    sep2ij[sep] = capo.dfm.grid2ij(aa.ant_layout)[0][sep].split(',')
    s += capo.dfm.grid2ij(aa.ant_layout)[0][sep] + ','
    
    for ij in sep2ij[sep]:
        ij2sep[ij] = sep

conj = {}
toconj = capo.dfm.grid2ij(aa.ant_layout)[1]
for k in toconj.keys():
    conj['%d_%d'%a.miriad.bl2ij(k)] = toconj[k]


chan = [120]
#times = {}
#data = {}
#flags = {}

_d = {}
_w = {}
curtime = 0
lij = []
for f in args:
    c = 0
    print f
    uv = a.miriad.UV(f)
    a.scripting.uv_selector(uv, ants=s, pol_str='xx')
    _d[f] = {}
    _w[f] = {}
    for sep in seps:
        _d[f][sep] = []
        _w[f][sep] = []
    for (uvw,t,(i,j)),d,w in uv.all(raw=True):
        print t
        if i==j: continue#dont include auto correlations
        if (curtime != 0) and (curtime != t): break #only include the first time
        curtime = t
        ijinsep = '%d_%d'%(i,j); ijinseprev = '%d_%d'%(j,i)
        lij.append(ijinsep)
        if ijinsep in ij2sep.keys():
            c+=1 
            print c, i,j
            if conj[ijinsep]:
                _d[f][ij2sep[ijinsep]].append(n.conjugate(d.take(chan)))
            else:
                _d[f][ij2sep[ijinsep]].append(d.take(chan))

print 'Finished reading data'

for f in args:
    for key in _d[f].keys():
        _d[f][key] = n.array(_d[f][key])
print _d

colors = ['c','m','y','b','g','r','k']
shapes = ['o','p','<','>','^','v','x','d','D']          
print len(_d.keys())

for e,f in enumerate(args):
    for i, key in enumerate(_d[f].keys()):
        x,y = divmod(i,7); print x,y
        dd = _d[f][key]
        print len(dd)
        print f, colors[y+e], shapes[x+e]
        p.scatter(dd.real, dd.imag,c=colors[y+e], marker=shapes[x+e], label=key+f[-1], s=30)

try:
    #to get npz file.
    chunk = args[0].split('.')[1] + '.' + args[0].split('.')[2][0] 
    print os.path.dirname(args[0])
    npzfile = glob.glob(os.path.dirname(args[0]) + '/data*%s*xx*.npz'%chunk)
    npz = n.load(npzfile[0])
    sepdict = {'0,1':'sep24', '-1,1':'sep6', '1,1':'sep35'}

    times = npz['times']
    time_index = n.where(n.round(times,4) == n.round(curtime,4))
    #import IPython 
    #IPython.embed()


    for sep in sepdict.keys():
        print time_index , chan
        mdldata = npz[sepdict[sep]][time_index,chan]
        print mdldata
        p.scatter(mdldata.real, mdldata.imag, c='black', marker='o', label=sep, s=40)

except:
    print 'Could not find omnical files'

p.legend()
p.show()
