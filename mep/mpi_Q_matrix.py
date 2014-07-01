import aipy as a, numpy as n, pylab as p
import capo as C
import useful_functions as uf
import global_sky_model as gsm
from scipy import special
import matplotlib as mpl


def get_single_Q_element((tx,ty,tz),amp,baseline,l,m):
	# split baseline
    bx,by,bz = baselines
	# compute spherical harmonic
    Ynorm = special.sph_harm(0,0,0,0)
    theta = n.arctan(ty/tx) # using math convention of theta=[0,2pi], phi=[0,pi]
    phi = n.arccos(n.sqrt(1-tx*tx-ty*ty))
    Y = n.array(special.sph_harm(m,l,theta,phi))/Ynorm #using math convention of theta=[0,2pi], phi=[0,pi]    
    #fringe pattern
    phs = n.exp(-2j*n.pi*fq * (bx*tx+by*ty+bz*tz)) 
    # not sure if I actually need this stuff
    valid = n.logical_not(tx.mask)
    Y.shape = phs.shape = amp.shape = im.uv.shape
    amp = n.where(valid, amp, 0)
    phs = n.where(valid, phs, 0)
    Y = n.where(valid, Y, 0) 
    Q_element = n.sum(amp*Y*phs)/n.sum(amp)
    return Q_element


