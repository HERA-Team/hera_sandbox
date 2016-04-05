#! /usr/bin/env python
import capo.hex as hx, capo.arp as arp, capo.red as red
import numpy as n, pylab as p, aipy as a
import sys

args = sys.argv[1:]

connection_file='/home/zakiali/src/mycapo/zsa/calfiles/HERA_connections.csv'
#hera info assuming a hex of 19 and 128 antennas
info = hx.hera_to_info(3, 128, connections=connection_file, ex_ants=[9,22])
reds = info.get_reds()
#seps are the representative baseline for a given 
#separation starting with the smallest antenna value.
use_hbls = [(0,1), (0,4)]

reds_dict = {}
paper_reds_dict = {}
paper_reds_conj_dict = {} #to conjugate a baseline or no in the paper layout
paper_to_string = []
red_to_paper = {} #keeps track of the key mapping of reds_dict to paper_reds_dict
for r in reds:
    #have key be the first baseline of a given type. Should be in numerical order.
    if not r[0] in use_hbls: continue
    reds_dict[r[0]] = r
    paper_bls = [ (paper_ants[bl[0]],paper_ants[bl[1]]) for bl in r ]
    red_to_paper[r[0]] = paper_bls[0]
    paper_reds_dict[paper_bls[0]] = paper_bls
    paper_to_string += paper_bls
    for bl in paper_bls:
        
        if bl[0] > bl[1]: 
            paper_reds_conj_dict[a.miriad.ij2bl(*bl)] = True
        else:
            paper_reds_conj_dict[a.miriad.ij2bl(*bl)] = False    

#Read in data here.
#generate antenna string for data reader.
#paper_to_string = [paper_reds_dict[r] for r in paper_reds_dict.keys()]
print paper_to_string
ants= ','.join(map(str,n.unique(n.array(paper_to_string).flatten())))
ant_string = '('+ants+')_('+ants+')'
info, data, flags = arp.get_dict_of_uv_data(args, ant_string, 'xx', verbose=True)
fqs = n.linspace(.1,.2,1024)

#will hold delay for each pairing of antnnas
d = {}
NANTS = len(n.unique(n.array(paper_to_string).flatten()))
print NANTS
for bl in reds_dict.keys():
    print reds_dict[bl]
    print paper_reds_dict[red_to_paper[bl]]

d1 = data[(80,104)]['xx']
d2 = n.conj(data[(96,104)]['xx'])
phase = red.redundant_bl_cal_simple(d1,d2,fqs,verbose=True)



#Do calibration here.



