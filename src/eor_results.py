import numpy as n,os
from capo import pspec
import glob
import ipdb
#measurements

def PAPER_32_all():
    '''
    Results from  PAPER 32  Jacobs et.al 2014

    outputs results[z] = n.array([k, Delta^2, 2-sigma upper, 2-sigma lower])
    '''

    PAPER_RESULTS_FILES = glob.glob(os.path.dirname(__file__)+'/data/psa32_apj/pspec_*.npz')
    PAPER_RESULTS_FILES.sort()
    freqs= []
    for filename in PAPER_RESULTS_FILES:
        try:
            freqs.append(n.load(filename)['freq']*1e3)
        except(KeyError):
            try:
                dchan = int(filename.split('/')[-1].split('.')[0].split('_')[2])-int(filename.split('/')[-1].split('.')[0].split('_')[1])
                chan = int(filename.split('/')[-1].split('.')[0].split('_')[1]) + dchan/2.
                freqs.append(chan/2. + 100) #a pretty good apprximation of chan 2 freq for 500kHz channels
            except: continue
    freqs = n.array(freqs)
    zs = pspec.f2z(freqs*1e-3)
    results = {}
    for files,z in zip(PAPER_RESULTS_FILES,zs):
        f=n.load(files)
        results[z] = n.array([f['k'],f['k3pk'],f['k3pk']+f['k3err'],f['k3pk']-f['k3err']]).T

    return results

def PAPER_64_all():
    '''
    Results from  PAPER 64  Ali et.al 2015

    outputs results[z] = n.array([k, Delta^2, 2-sigma upper, 2-sigma lower])
    '''

    PAPER_RESULTS_FILES = glob.glob(os.path.dirname(__file__)+'/data/psa64_apj/pspec_*.npz')
    PAPER_RESULTS_FILES.sort()
    freqs = []
    for filename in PAPER_RESULTS_FILES:
        try:
            freqs.append(n.load(filename)['freq']*1e3)
        except(KeyError):
            try:
                dchan = int(filename.split('/')[-1].split('.')[0].split('_')[2])-int(filename.split('.')[0].split('_')[1])
                chan = int(filename.split('/')[-1].split('.')[0].split('_')[1]) + dchan/2.
                freqs.append(chan/2. + 100) #a pretty good apprximation of chan 2 freq for 500kHz channels
            except: continue
    freqs = n.array(freqs)
    zs = pspec.f2z(freqs*1e-3)
    #zs = n.array([8.31])
    results = {}
    for files,z in zip(PAPER_RESULTS_FILES,zs):
        f=n.load(files)
        results[z] = n.array([f['k'],f['k3pk'],f['k3pk']+f['k3err'],f['k3pk']-f['k3err']]).T

    return results

def MWA_32T_all():
    """
    From Josh Dillon April 10, 2014
    Delta^2  |  2sig bottom bar  |  2sig top bar  |  k  |  z | "Detection"
    
    or
    
    2sig top bar   |  [ignore this entry]  |  2sig top bar  |  k  |  z | "2sigUL"
    
    in the units you requested.
    
    (my request The ideal format would be a text file listing three things:
    Delta^2 [mk^2], k [hMpc^-1], z [z])
    
    
    TREATMENT:
    I will store the delta^2 +/- values.  If he doesn't give me a delta^2 or lower limit, I will put 
    Delta^2 = 0 and bottom = -top
    
    return format will be dict[z] = n.array([[k,Delta^2,top,bottom]]) all in mK^2
    
    """
    

    MWA_RESULTS_FILE=os.path.dirname(__file__)+'/data/MWA_32T_Results_Summary.txt'
    lines = open(MWA_RESULTS_FILE).readlines()
    result = {}
    for line in lines:
        if line.endswith('2sigUL'):#if its just an upper limit
            l = map(float,line.split()[:-1])
            z = l[-1]
            try:
                result[z].append([l[3],0,l[0],-l[0]])
            except(KeyError):
                result[z] = [[l[3],0,l[0],-l[0]]]
        else:#if its a "detection" which just means its not consistent with 0, not really EoR but probably foregrounds
            #or systematics
            l = map(float,line.split()[:-1])
            z = l[-1]
            try:
                result[z].append([l[3],l[0],l[0],-l[0]])
            except(KeyError):
                result[z] = [[l[3],l[0],l[2],l[1]]]
    for z in result.keys():
        result[z] = n.array(result[z])
    return result

def MWA_128_all():
    '''
    MWSA_128 data from dillion 2015
    return format will be dict[z] = n.array([[k,Delta^2,top,bottom]]) all in mK^2
    '''
    MWA_RESULTS_FILE=glob.glob(os.path.dirname(__file__)+'/data/mwa128/*.dat')
    zs= []
    results = {}
    for files in MWA_RESULTS_FILE:
        name = files.split('/')[-1]
        nums = name.split('=')[1].split('.')[:2]
        z= float(nums[0]+'.'+nums[1])

        data = n.genfromtxt(files, delimiter=' ')
        results[z] = data[:,[0,1,5,4]]

    return results
def z_slice(redshift,pspec_data):
    """
    input a power spectrum data dict output of MWA_32T_all() or GMRT_2014_all()
    returns a slice along k for the input redshift 
    example
    z,pspec[k,k3pK] = z_slice(MWA_32T_all())
    """
    zs = n.array(pspec_data.keys())
    closest_z = zs[n.abs(zs-redshift).argmin()]
    return closest_z, pspec_data[closest_z]
def k_slice(k,pspec_data):
    """
    input a power spectrum data dict output of MWA_32T_all() or GMRT_2014_all()
    returns a slice along z for the input redshift 
    example
    zs,pspec[k,k3pK] = k_slice(MWA_32T_all())
    """

    zs = n.array(pspec_data.keys())
    k_is = [n.abs(pspec_data[redshift][:,0]-k).argmin() for redshift in zs]
    ks = [pspec_data[redshift][k_i,0] for k_i in k_is]
    power = n.vstack([pspec_data[redshift][k_i,:] for k_i in k_is])
    return zs,power
    

def MWA_32T_at_z(redshift):
    """
    Input a redshift
    Return the Delta^2(k) power spectrum at the redshift nearest the input redshift
    returns: z,array([k,Delta^2,2sig_upper_limit,2sig_lower_limit])
    """
    data = MWA_32T_all()
    zs = n.array(data.keys())
    closest_z = zs[n.abs(zs-redshift).argmin()]
    return closest_z, data[closest_z]
def MWA_32T_at_k(k):
    """
    Input a k
    Returns the Delta^2(k) power spectrum vs redshift at that k. (all pspec in mk^2)
    return format: z,array([k,Delta^2,2sig_upper_limit,2sig_lower_limit])
    """
    data = MWA_32T_all()
    zs = n.array(data.keys())
    k_is = [n.abs(data[redshift][:,0]-k).argmin() for redshift in zs]
    ks = [data[redshift][k_i,0] for k_i in k_is]
    power = n.vstack([data[redshift][k_i,:] for k_i in k_is])
    return zs,power
def GMRT_2014_all():
    return {8.6:n.array(
            [[0.1, 0,2e5,  -2e5],
            [0.13, 0,4e5,  -4e5],
            [0.16, 0,1e5,  -1e5],
            [0.19, 0,1.9e5,-1.9e5],
            [0.275,0,2e5,  -2e5],
            [0.31, 0,4e5,  -4e5],
            [0.4,  0,6e5,  -6e5],
            [0.5,  0,8e4,  -8e4]])}
