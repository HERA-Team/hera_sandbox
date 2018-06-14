#!/usr/bin/env python

#SBATCH -J hera_pspec
#SBATCH --mem=5G
#SBATCH -t 3:00:00
#SBATCH -n 4
#SBATCH --mail-type=FAIL     
#SBATCH --mail-user=adam_lanman@brown.edu


"""
Run HERA OQE power spectrum estimation code on sets of redundant baselines.
Modified from pspec_red.py to use FHD results.
"""
import numpy as np
import hera_pspec as hp
from hera_pspec.utils import log as hlog, load_config
from hera_cal import redcal
import pyuvdata as uv
import os, sys, glob, time
from scipy.io import readsav


def log(string,**kwargs):
    hlog(string,**kwargs)
    sys.stdout.flush()

# Default settings for pspec calculation
pspec_defaults = {
    'overwrite':                False,
    'little_h':                 False,
    'exclude_auto_bls':         False,
    'exclude_permutations':     False,
}

###TODO:
###     Figure out where visibility-based addition is happening in the pipeline.
###              (it should be possible to get separate pspecs from different files and average them as pspecs.
###     Once I'm more confident in my results, start sending things to Adrian.
###     Upgrade the plotter a little
    
#-------------------------------------------------------------------------------
# Settings
#-------------------------------------------------------------------------------

# Get configuration filename from cmdline argument
if len(sys.argv) > 1:
    cfg_file = str(sys.argv[1])
else:
    print("Command takes one argument: config_file")
    sys.exit(1)

# Load configuration file
cfg = load_config(cfg_file)
data_cfg = cfg['data']
pspec_cfg = cfg['pspec']


#-------------------------------------------------------------------------------
# Prepare list of data files
#-------------------------------------------------------------------------------

files = []
for i in range(len(data_cfg['subdirs'])):
    files += glob.glob( os.path.join(data_cfg['root'], 
                                     data_cfg['subdirs'][i], 
                                     data_cfg['template']) )
for f in files:
    print(f)

log("Found %d files." % len(files))

#-------------------------------------------------------------------------------
# Load data files into memory
#-------------------------------------------------------------------------------
log("Loading data files...")
t0 = time.time()

# Load all miriad datafiles into UVData objects
dsets = []
data = None
for f in files:
    _d = uv.UVData()
    _d.read_miriad(f)
#    _d.polarization_array = uv.utils.polnum2str(_d.polarization_array)
    if _d.vis_units == "JY":
        _d.vis_units = 'Jy'
    _d.object_name = 'Zenith'
    if data is None: data = _d
    else: data += _d
#    dsets.append(_d)
log("Loaded data in %1.1f sec." % (time.time() - t0), lvl=1)

if 'antnum_select' in data_cfg.keys():
    data.select(antenna_nums=data_cfg['antnum_select'])
    inds = data_cfg['antnum_select']
    data.antenna_numbers = data.antenna_numbers[inds]
    if len(data.antenna_names) == data.Nants_telescope:
        data.antenna_names = np.array(data.antenna_names)[inds]
    else:
        data.antenna_names = np.array([ 'ANT'+num for num in data_cfg['antnum_select']])
    data.Nants_telescope= len(data.antenna_numbers)


#-------------------------------------------------------------------------------
# If antenna positions are missing, load from an npz file.
#-------------------------------------------------------------------------------

if data.antenna_positions is None:
    try:
        if not os.path.exists(data_cfg['antpos_file']):
            path = os.path.dirname(sys.argv[1])
            aposf = os.path.join(path, data_cfg['antpos_file'])
            if not os.path.exists(aposf):
                print aposf
                raise OSError("Cannot find antenna positions file")
            data_cfg['antpos_file'] = aposf
        antpos = np.load(data_cfg['antpos_file'])['antpos']
        data.antenna_positions = antpos
    except IOError:
        print("Antenna positions npz file must be given.")

# Convert to pseudostokes I
print 'XX', np.mean(data.data_array[:,:,:,0])
data = hp.pstokes.construct_pstokes(data, data)
print 'pI', np.mean(data.data_array[:,:,:,0])

# Generate beam

if 'obs' in data_cfg.keys():
    obsfile = data_cfg['obs']
else:
    #Take the first recognized obs file.
    obsfile = glob.glob( os.path.join(data_cfg['root'], 'metadata/*_obs.sav') )[0]
 
obs = readsav(obsfile)['obs']
Op = obs['PRIMARY_BEAM_AREA'][0]   # Shape = (npol, nfreq), type=float32. This is from FHD.
Opp = obs['PRIMARY_BEAM_SQ_AREA'][0]

polnames = ['XX','YY','XY','YX']
#polnums = [-5,-6,-7,-8]
Op_dict, Opp_dict = {}, {}
for pi in range(4):
    if Opp[pi] is None or Op[pi] is None:
        continue
    Op_dict[polnames[pi]] = Op[pi]
    Opp_dict[polnames[pi]] = Opp[pi]

### Construct pseudostokes I beam
#Op_dict['pI'] = 0.5*(Op_dict['XX'] + Op_dict['YY'])
#Opp_dict['pI'] = 0.25*(Opp_dict['XX'] + Opp_dict['YY'] + 2*Op_dict['XX']*Op_dict['YY'])  #!!!Assumes XX and YY are identical for now.
Op_dict['pI'] = 0.5*(Op_dict['XX'] + Op_dict['YY'])
Opp_dict['pI'] = 0.25*(Opp_dict['XX'] + Opp_dict['YY'] + 2*Opp_dict['YY'])  #!!!Assumes XX and YY are identical for now.
print np.mean(Op_dict['pI']), np.mean(Op_dict['XX'])
print np.mean(Opp_dict['pI']), np.mean(Opp_dict['YY'])

n_freq_obs = obs['n_freq'][0]
freq_center_obs = obs['freq_center'][0]
freq_res_obs = obs['freq_res'][0]
off = (n_freq_obs/2.)*freq_res_obs
obs_freq = np.linspace( freq_center_obs-off, freq_center_obs+off, n_freq_obs)
### Round to nearest kHz
obs_freq = np.floor(obs_freq/1e3)*1e3
beam = hp.PSpecBeamFromArray(Op_dict, Opp_dict, obs_freq)

log("Generated beam from obs: %s" % obsfile)

# Convert data files from Jy to mK if requested
if 'convert_jy_to_mk' in data_cfg.keys():
    if data_cfg['convert_jy_to_mk']:
        freqs = data.freq_array.flatten()
        data.data_array[:,:,:,0] *= beam.Jy_to_mK(freqs, pol='pI')[None, None, :]
#        for poli, polstr in enumerate(Op_dict.keys()):
#            data.data_array[:,:,:,poli] *= beam.Jy_to_mK(freqs,pol=polstr)[None, None, :]
#        for i in range(len(dsets)):
#            freqs = dsets[i].freq_array.flatten()
#            for poli, polstr in enumerate(Op_dict.keys()):
#                dsets[i].data_array[:,:,:,poli] *= beam.Jy_to_mK(freqs,pol=polstr)[None, None, :]
#            dsets[i].vis_units = 'mK'

#-------------------------------------------------------------------------------
# Calculate power spectrum and package output into PSpecContainer
#-------------------------------------------------------------------------------

for key in pspec_defaults.keys():
    if key in pspec_cfg.keys(): pspec_defaults[key] = pspec_cfg[key]

antpos, ants = data.get_ENU_antpos()
antpos = dict(zip(ants, antpos))
red_bls = redcal.get_pos_reds(antpos, bl_error_tol=1.0)#, low_hi=True)
## FIXME: Use only the first redundant baseline group for now (shortest baselines)
#
### Fix conjugation error:
# (might not be necessary since we're not averaging visibilities)
bls = []
for bl in red_bls[0]:
    if bl[0] > bl[1]:
        bl_ind = data.antnums_to_baseline(bl[1],bl[0])
        inds = data.baseline_array == bl_ind
        data.data_array[inds,:,:,:] = data.data_array[inds,:,:,:].conj()
        bl = (bl[1], bl[0])
    bls.append(bl)
print("Baselines: %s" % bls)

bls1, bls2, blpairs = hp.utils.construct_blpairs(bls, exclude_auto_bls=pspec_cfg['exclude_auto_bls'], exclude_permutations=pspec_cfg['exclude_permutations'])

## Add together data files into one dataset (for data lengths less than 24 hours, equivalent to lst binning)
#data = None
#for d in dsets:
#    if data is None: data = d
#    else: data += d
#
ds = hp.PSpecData(dsets=[data, data], wgts=[None, None], beam=beam)

ps_store = hp.PSpecContainer(pspec_cfg['output'], mode='rw')


ps_ii = ds.pspec(bls1, bls2, (0,0), ("pI","pI"),
                      input_data_weight=pspec_cfg['weight'],
                      norm=pspec_cfg['norm'], 
                      taper=pspec_cfg['taper'], 
                      spw_ranges=None,
                      little_h=pspec_cfg['little_h'])

avg_pspec_ii = ps_ii.average_spectra(blpair_groups=[blpairs], inplace=False)

np.savez('pspec_results', ps_ii=ps_ii, avg_pspec_ii=avg_pspec_ii)

ps_store.set_pspec(group=pspec_cfg['groupname'], psname='ps_ii', pspec=ps_ii, overwrite=pspec_defaults['overwrite'])

## Loop over pairs of datasets
#dset_idxs = range(len(ds.dsets))
#for i in dset_idxs:
#    for j in dset_idxs:
#        if i >= j: continue
#        
#        # Name for this set of power spectra
#        pspec_name = "pspec_dset(%d,%d)" % (i,j)
#        
#        # Calculate power spectra for all baseline pairs (returns UVPSpec)
#        ps = ds.pspec(bls1, bls2, (i,j), ('XX','XX'),
#                      input_data_weight=pspec_cfg['weight'],
#                      norm=pspec_cfg['norm'], 
#                      taper=pspec_cfg['taper'], 
#                      spw_ranges=None,
#                      little_h=pspec_defaults['little_h'])
#        
#        # Store power spectra in container
#        ps_store.set_pspec(group=pspec_cfg['groupname'], psname=pspec_name, 
#                           pspec=ps, overwrite=pspec_defaults['overwrite'])
#
ps_store.tree()

