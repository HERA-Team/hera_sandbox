import numpy as np
import capo.omni as omni

def hex_to_info(hexnum, pols=['x'], **kwargs):
    '''Given a hex number generate redundancies, assuming antenna
       numbering starts at 0 in the bottom left of the hex.
       Returns an omnical info class'''
    nant = 3*hexnum**2 - 3*hexnum + 1
    antpos = -np.ones((nant*len(pols),3))
    ant=0
    for row in range(hexnum-1, -(hexnum), -1):
        for col in range(2*hexnum-abs(row)-1):
            for z, pol in enumerate(pols):
                x = ((-(2*hexnum-abs(row))+2)/2.0 + col)
                y = row*-1*np.sqrt(3)/2
                i = omni.Antpol(ant,pol,nant)
                antpos[i,0],antpos[i,1],antpos[i,2] = x,y,z
            ant+=1
    reds = omni.compute_reds(nant, pols, antpos[:nant],tol=.1)
    ex_ants = []
    kwargs['ex_ants'] = kwargs.get('ex_ants',[]) + ex_ants
    reds = omni.filter_reds(reds, **kwargs)
    info = omni.RedundantInfo(nant)
    info.init_from_reds(reds,antpos)
    return info

def hera_to_info(hexnum, nant, pols=['x'], connections=None, **kwargs):
    '''Go from hera array to info. Use the fact that the 
        hex < nant in ZA to get the reds and ordering correct.'''
    antpos = -np.ones((nant*len(pols),3))
    if connections:
        paper_ants = get_paper_ants(connections)
    hex_index = 0
    for row in range(hexnum-1, -(hexnum), -1):
        for col in range(2*hexnum-abs(row)-1):
            for z, pol in enumerate(pols):
                z = 2**z
                x = ((-(2*hexnum-abs(row))+2)/2.0 + col)
                y = row*-1*np.sqrt(3)/2
                i = omni.Antpol(paper_ants[hex_index],pol,nant)
                antpos[i,0],antpos[i,1],antpos[i,2] = x,y,z
            hex_index+=1
    reds = omni.compute_reds(nant, pols, antpos[:nant],tol=.1)
    ex_ants = [omni.Antpol(i,nant).ant() for i in range(antpos.shape[0]) if antpos[i,-1] < 0] # look at z (pol) component to see if flagged antenna. Positions can be positive and negative.
    #import IPython; IPython.embed()
    kwargs['ex_ants'] = kwargs.get('ex_ants',[]) + ex_ants
    reds = omni.filter_reds(reds, **kwargs)
    info = omni.RedundantInfo(nant)
    info.init_from_reds(reds,antpos)
    return info

    
 
def get_paper_ants(connection_file):
    '''Get the PAPER antenna numbers for a hex grid.
       Assumes that the antenna numbers are from 0-N
       Uses a connections file
            16 17 18
           12 13 14 15
          7  8  9  10 11
            3  4  5  6
             1  2  3....for example.    
            '''
    connections = np.recfromcsv(connection_file)
    hex_stations = connections['station']
    hex_stations_argsort = np.argsort(hex_stations)
    paper_ants = connections['paper'][hex_stations_argsort]
    return paper_ants

