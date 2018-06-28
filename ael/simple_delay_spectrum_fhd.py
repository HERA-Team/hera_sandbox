
### Do delay spectrum analysis without using hera_pspec

import capo.pspec
from scipy.io import readsav
import numpy as np
import sys, yaml, glob, os, time
from eorsky import comoving_voxel_volume 
from pyuvdata import UVData

### Define cosmological functions (adapted from capo):

def jy2T(f, bm):
    '''Return [mK] / [Jy] for a beam size vs. frequency (in Hz)
        f = frequencies (Hz!)
        bm = primary beam integral at corresponding f (Omega_p)
    '''
    c_cmps = 2.99792458e10   # cm/s
    k_boltz = 1.380658e-16   # erg/K
    lam = c_cmps / f   #cm
    return 1e-23 * lam**2 / (2 * k_boltz * bm) * 1e3

def dL_df(z, omega_m=0.266):
    '''[ Mpc]/Hz, from Furlanetto et al. (2006)'''
    #### IF EVERYTHING'S OFF BY A LOT, CHECK THIS FIRST.
    return (1.7 / 0.1) * ((1+z) / 10.)**.5 * (omega_m/0.15)**-0.5 * 1e3 * 1e-9 * 0.7

def dL_dth(z):
    '''[Mpc]/radian, from Furlanetto et al. (2006)'''
    arcmin = (np.pi/180.)*(1/60.)
    return 1.9 * (1./arcmin) * ((1+z) / 10.)**.2 * 0.7

def dk_deta(z):
    '''2pi * [Mpc^-1] / [Hz^-1]'''
    return 2*np.pi / dL_df(z) 

def X2Y(z):
    '''[Mpc^3] / [str * Hz] scalar conversion between observing and cosmological coordinates'''
    return dL_dth(z)**2 * dL_df(z)

### Parse command line args

if len(sys.argv) > 1:
    config_file = str(sys.argv[1])
else:
    print("Command takes one argument: config_file")
    sys.exit(1)

with open(config_file, 'r') as cfile:
    try:
        cfg = yaml.load(cfile)
    except yaml.YAMLError as exc:
        raise(exc)

data_cfg = cfg['data']
pspec_cfg = cfg['pspec']


#-------------------------------------------------------------------------------
# Get obs (and Op, Opp) from first metadata
#-------------------------------------------------------------------------------
if 'obs' in data_cfg.keys():
    obsfile = data_cfg['obs']
else:
    #Take the first recognized obs file.
    obsfile = glob.glob( os.path.join(data_cfg['root'], 'metadata/*_obs.sav') )[0]
 
obs = readsav(obsfile)['obs']
Op = obs['PRIMARY_BEAM_AREA'][0]   # Shape = (npol, nfreq), type=float32. This is from FHD.
Opp = obs['PRIMARY_BEAM_SQ_AREA'][0]

n_freq_obs = obs['n_freq'][0]
freq_center_obs = obs['freq_center'][0]
freq_res_obs = obs['freq_res'][0]
off = (n_freq_obs/2.)*freq_res_obs
freq_obs = np.linspace( freq_center_obs-off, freq_center_obs+off, n_freq_obs)
Bandwidth = freq_obs.size * freq_res_obs

#-------------------------------------------------------------------------------
# Prepare list of data files
#-------------------------------------------------------------------------------
files=[]
for i in range(len(data_cfg['subdirs'])):
    files += glob.glob( os.path.join(data_cfg['root'], 
                                     data_cfg['subdirs'][i], 
                                     data_cfg['template']) )
for f in files:
    print(f)

print("Found %d files." % len(files))

#### Locating the redundant baselines
#if data.antenna_positions is None:
#    try:
#        if not os.path.exists(data_cfg['antpos_file']):
#            path = os.path.dirname(sys.argv[1])
#            aposf = os.path.join(path, data_cfg['antpos_file'])
#            if not os.path.exists(aposf):
#                print aposf
#                raise OSError("Cannot find antenna positions file")
#            data_cfg['antpos_file'] = aposf
#        antpos = np.load(data_cfg['antpos_file'])['antpos']
#        data.antenna_positions = antpos
#    except IOError:
#        print("Antenna positions npz file must be given.")

# Load all miriad datafiles into UVData objects


### Select baselines from parameters:

bl = eval(pspec_cfg['baseline'])
bls = [bl]


#-------------------------------------------------------------------------------
# Load data files into memory
#-------------------------------------------------------------------------------
print("Loading data files...")
t0 = time.time()

files.sort()

dsets = []
data = None
count = 0
if len(files) == 0:
    raise ValueError("Couldn't find any files")
for f in files:
    print count, f
    _d = UVData()
    _d.read_miriad(f, bls=bls)
#    _d.polarization_array = uv.utils.polnum2str(_d.polarization_array)
    if _d.vis_units == "JY":
        _d.vis_units = 'Jy'
    _d.object_name = 'Zenith'
    if data is None: data = _d
    else: data += _d
    count+=1
#    dsets.append(_d)
print("Loaded data in %1.1f sec." % (time.time() - t0))


# Get k_parallel
Zs = 1420e6/data.freq_array[0] - 1
Zmean = np.mean(Zs)
etas = np.fft.fftfreq(data.Nfreqs, d=data.channel_width)
k_parallel = dk_deta(Zmean) * etas


vis_xx = data.data_array[:,0,:,0]
vis_yy = data.data_array[:,0,:,1]

# NOTE  Assuming Op/Opp of pI are the same as Op/Opp of XX

Op_I = Op[0]# + Op[1]
Opp_I = Opp[0]# + Opp[1] + 2*Op[0]*Op[1]  #Both have shape (Nfreqs)

vis_I = vis_xx + vis_yy

#Convert to mK*str
vis_I = vis_I * jy2T(data.freq_array, Op_I)

# FFT

_visI = np.fft.ifft(vis_I, axis=1)

dspec_instr = abs(_visI)**2


# Cosmological Normalization:
Oeff = Op_I**2/Opp_I
scalar = X2Y(Zmean) * Oeff * Bandwidth

dspec_I = dspec_instr * scalar

nside=512.
om = 4*np.pi/(12.*nside**2)
dnu = np.diff(freq_obs)[0]/1e6   #MHz
ref_line = comoving_voxel_volume(Zmean, dnu, om) * 4e6


f = readsav('/users/alanman/FHD/catalog_data/eor_power_1d.idlsave')
maxpow = max(f['power'])
ref_line = maxpow


##Next --- construct vis, Op, and Opp for pseudo-stokes-I.
##         Convert Jy to mK
##         Do FFTs, squaring, etc.
##         Save to outfile: dspec_ii (ntimes, Nk) , k_parallel, cosmo scalars, Op, Opp, vis_list
 
