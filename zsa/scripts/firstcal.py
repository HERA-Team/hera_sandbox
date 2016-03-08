#! /usr/bin/env python
import capo.hex as hx, capo.arp as arp
import numpy as n, pylab as p
import aipy as a



connection_file='../calfiles/HERA_connections.csv'
paper_ants = hx.get_paper_ants(connection_file)
info = hx.hex_to_info(3,ex_ants=[4,8,9,10,12,13,14] )#hex info with hex number 3 = 19 antennas
reds = info.get_reds()
reds_dict = {}
paper_reds_dict = {}
paper_reds_conj_dict = {}
paper_to_string = []
for r in reds:
    #have key be the first baseline of a given type. Should be in numerical order.
    reds_dict[r[0]] = r
    paper_bls = [ (paper_ants[bl[0]],paper_ants[bl[1]]) for bl in r ]
    paper_reds_dict[paper_bls[0]] = paper_bls
    paper_to_string += paper_bls
    for bl in paper_bls:
        
        if bl[0] > bl[1]: 
            paper_reds_conj_dict[a.miriad.ij2bl(*bl)] = True
        else:
            paper_reds_conj_dict[a.miriad.ij2bl(*bl)] = False    

#Read in data here.
ant_string = ','.join(map(str,n.unique(n.array(paper_to_string).flatten())))
info, data, flags = arp.get_dict_of_uv_data(args, ant_string, 'xx')

import IPython; IPython.embed()


#Do calibration here.



