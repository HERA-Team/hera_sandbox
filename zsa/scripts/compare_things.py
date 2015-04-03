#!/usr/bin/env python
import numpy as n
import capo as c
import pylab as p
import sys, glob

#afiles = glob.glob('/data2/home/aparsons/omnical_results/data_psa6246_2456246.2_xx_add7.omnical.npz')
afiles = glob.glob('/data2/home/aparsons/omnical_results/omnical_lstbin.npz')
#bfiles = glob.glob('/data2/home/zakiali/psa_live/forlstbinning_omnical_2/data_psa6246_2456246.2_xx_add7.omnical.npz')
bfiles = glob.glob('/data2/home/zakiali/psa_live/forlstbinning_omnical_2/omnical_lstbin.npz')

ad = [n.load(f) for f in afiles]
bd = [n.load(f) for f in bfiles]


adict = {} 
bdict = {}
dicts = [adict, bdict]

for i,g in enumerate([ad,bd]):
    for file in g:
        for key in file.keys():
            if not (key in dicts[i]):
                dicts[i][key] = file[key]
            else:
                dicts[i][key] = n.concatenate([dicts[i][key], file[key]])

for key in adict.keys():
    print key
    try:
        adict[key] = n.where(n.isnan(adict[key]), 0, adict[key])
        bdict[key] = n.where(n.isnan(bdict[key]), 0, bdict[key])
    except:
        continue
         
import IPython.Shell; IPython.Shell.IPShellEmbed('')()
