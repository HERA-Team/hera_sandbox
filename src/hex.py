import numpy as np
import capo.omni as omni

def hex_to_info(hexnum, pols=['x'], **kwargs):
    '''Given a hex number generate redundancies, assuming antenna
       numbering starts at 0 in the bottom left of the hex.
       Returns an omnical info class'''
    nant = 3*hexnum**2 - 3*hexnum + 1
    hindices = np.arange(1, nant+1)
    positions = []
    antpos = -np.ones((nant*len(pols),3))
    ant=0
    for row in range(hexnum-1, -(hexnum), -1):
        for col in range(2*hexnum-abs(row)-1):
            for z, pol in enumerate(pols):
                x = ((-(2*hexnum-abs(row))+2)/2.0 + col)
                y = row*-1*np.sqrt(3)/2
                i = Antpol(ant,pol,nant)
                antpos[i,0],antpos[i,1],antpos[i,2] = x,y,z
            ant+=1
    reds = omni.compute_reds(nant, pols, antpos[:nant],tol=.1)
    ex_ants = []
    kwargs['ex_ants'] = kwargs.get('ex_ants',[]) + ex_ants
    reds = omni.filter_reds(reds, **kwargs)
    info = omni.RedundantInfo(nant)
    info.init_from_reds(reds,antpos)
    return info
 

