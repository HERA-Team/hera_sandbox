#input a list of signal loss power spectra (as output by boots_to_pspsec)
#compute a loss-calibrated limit
#output as a new pspec file
from glob import glob
import argparse,os
from capo.eor_results import read_bootstraps_dcj,average_bootstraps
from capo.pspec import dk_du
from capo import cosmo_units
import numpy as np
from matplotlib.pyplot import *
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
def injpath_to_injlevel(filepath):
    injdir = os.path.split(os.path.dirname(filepath))[-1]
    return float(injdir.split('_')[-1])
#read in the bootstrapped pspecs
pspec_channels = ['pCv_fold','pCv_fold_err', #weighted data pspec
                    'pC_fold','pC_fold_err', #weighted data+inj pspec
                    'pI_fold','pI_fold_err'] #unweighted inj pspec
pspecs = {}
#sort the input files. makes things easier later
injlevels = [injpath_to_injlevel(filename) for filename in args.pspec_files]
fileorder = np.argsort(injlevels)
filenames = [args.pspec_files[i] for i in fileorder]
for filename in filenames:
    print filename
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
        lossy_limit = pspecs['pCv_fold'][0,k_ind]+pspecs['pCv_fold_err'][0,k_ind]
        if k_ind==5:
            print "lossy_limit: ", lossy_limit
            print "pC lower limit",pspecs['pC_fold'][inj,k_ind]-pspecs['pC_fold_err'][inj,k_ind]
        probs[inj,k_ind] = G_mc(lossy_limit, #limit
                        pspecs['pC_fold'][inj,k_ind],
                        pspecs['pC_fold_err'][inj,k_ind])
for k_ind in xrange(Nk):
    loglog(pspecs['pI_fold'][:,k_ind],probs[:,k_ind],'-',label=k[k_ind])
legend(loc='best')
show()
