#! /usr/bin/env python
import capo.hex as hx, capo.arp as arp, capo.red as red, capo.omni as omni
import numpy as n, pylab as p, aipy as a
import sys

args = sys.argv[1:]

connection_file='HERA_connections.csv'
#hera info assuming a hex of 19 and 128 antennas
#info = hx.hera_to_info(3, 128, connections=connection_file, ex_ants=[9,22,89,53], ubls=[(80,104),(53,80)])
info = hx.hera_to_info(3, 128, connections=connection_file, ubls=[(80,104),(53,80)])
reds = info.get_reds()

#Read in data here.
#generate antenna string for data reader.
ant_string =','.join(map(str,info.subsetant))
times, data, flags = arp.get_dict_of_uv_data(args, '('+ant_string+')_('+ant_string+')', 'xx', verbose=True)
dataxx = {}
flagsxx = {}
for k in data.keys():
    dataxx[k] = data[k]['xx']
for k in flags.keys():
    flagsxx[k] = flags[k]['xx']
fqs = n.linspace(.1,.2,1024)

#gets phase solutions per frequency.
fc = omni.FirstCal(dataxx,fqs,info)
sols = fc.run()
dataxx_c = {}
for (a1,a2) in dataxx.keys():
    dataxx_c[(a1,a2)] = dataxx[(a1,a2)]*omni.get_phase(fqs,sols[a1])*n.conj(omni.get_phase(fqs,sols[a2]))

import IPython; IPython.embed()
   





