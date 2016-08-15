#input a list of signal loss power spectra (as output by boots_to_pspsec)
#compute a loss-calibrated limit
#output as a new pspec file
from glob import glob
import argparse
from capo.eor_results import read_bootstraps_dcj,average_bootstraps
from capo.pspec import dk_du
from capo import cosmo_units
import numpy as np

parser = argparse.ArgumentParser(description='Calculate a limit from a range of injected pspecs')
parser.add_argument('pspec_files', metavar='pspsec.npz', type=str, nargs='+',default='.',
                    help='Directory containing injections')
args = parser.parse_args()

def G(x,mx,dx):
    return 1/(dx*np.sqrt(2*np.pi)) * np.exp(-1/2.*((x-mx)/dx)**2)
def G_mc(x_lim,dist_x,dist_sig,Ntrials=1000):
    #calculate the probability of drawing a point larger than
    # x_lim from the gaussian distribution dist_x +/- dist_sig
    #
    mc_x = np.random.normal(dist_x,dist_sig,size=Ntrials)
    return np.sum(mc_x>x_lim)/float(Ntrials)

#read in the bootstrapped pspecs
pspec_channels = ['pCv_fold','pCv_fold_err', #weighted data pspec
                    'pC_fold','pC_fold_err', #weighted data+inj pspec
                    'pI_fold','pI_fold_err'] #unweighted inj pspec
pspecs = {}
for filename in args.pspec_files:
    F = np.load(filename)
    for pspec_channel in pspec_channels:
        F[pspec_channel]
        try:
            pspecs[pspec_channel].append(F[pspec_channel])
        except(KeyError):
            pspecs[pspec_channel] = [F[pspec_channel]]
for pspec_channel in pspec_channels:
    pspecs[pspec_channel] = np.array(pspecs[pspec_channel])
Ninj,Nk = pspecs[pspec_channel].shape
k = F['k']

print "found injects:",Ninj
print "found k bins:",Nk


#limit option #1. the net probability that pC is above pcV
probs = np.zeros((Ninj,Nk))
for inj in xrange(Ninj):
    for k_ind in xrange(Nk):
        lossy_limit = pspecs['pCv_fold'][k_ind]+pspecs['pCv_fold_err'][k_ind]
        probs[inj,k_ind] = G_mc(lossy_limit, #limit
                        pspecs['pC_fold'][k_ind],pspecs['pC_fold_err'][k_ind])
for k_ind in xrange(Nk):
    loglog(k,probs[:,k_ind],'-')
show()
