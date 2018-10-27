## Generic delay spectrum calculation.
## Using a config file, load files and save a delay spectrum

import numpy as np
import sys
import os
import glob
import yaml
import optparse
import time
from pyuvdata import UVData
from scipy.io import readsav
from sandbox_cosmo_funcs import *

### Parse command line args

o = optparse.OptionParser()
o.add_option('-r', "--random_vis", action='store_true', help="Replace data with random visibilities of comparable variance.")
o.add_option('-z', "--zero_means", action='store_true', dest='zero', help="Subtract visibility mean.", default=False)
o.add_option('-p', "--partition_averaging", dest='part', action='store_true', help="Partition the time series, average partitions")
o.add_option('-b', '--Nboots', dest='Nboots', help="Use bootstrap averaging instead, specifying number of boots here.", default=None)
o.add_option('-s', '--max_step', dest='max_step', help="Use decimation steps, up to maximum step length specified here.", default=None)
o.add_option('-c', '--covariance', dest='covar', action='store_true', help='Calculate and return the covariance of the data array')

opts,args = o.parse_args(sys.argv[1:])

partition = opts.part
over_k = not opts.part
multislice=False
if opts.Nboots is not None:
    Nboots = int(opts.Nboots)
    boot = True
    over_k = False
    partition = False
if opts.max_step is not None:
    max_step = int(opts.max_step)
    multislice = True

if len(args) == 1:
    config_file = str(args[0])
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
out_cfg = cfg['out']
try:
    catalog_cfg = cfg['catalog']
except KeyError:
    catalog_cfg = {}

out_cfg.update(opts.__dict__)

save_dict={}    # results to save to npz

#-------------------------------------------------------------------------------
# Prepare list of data files
#-------------------------------------------------------------------------------


files=[]
if not 'subdirs' in data_cfg:
    data_cfg['subdirs'] = ['']
for i in range(len(data_cfg['subdirs'])):
    files += glob.glob( os.path.join(data_cfg['root'], 
                                     data_cfg['subdirs'][i], 
                                     data_cfg['template']) )
print("Found %d files." % len(files))


# Select baselines from parameters:
bl = eval(pspec_cfg['baseline'])
bls = [bl]

#-------------------------------------------------------------------------------
# Load data files into memory
#-------------------------------------------------------------------------------
print("Loading data files...")
t0 = time.time()

#def sortfunc(s):
#    s = s.split('/')[-1]
#    try:
#        return int(s.split('_')[2])
#    except KeyError:
#        # FHD files have an index in this column
#        return s
#
#files.sort(key=sortfunc)

dsets = []
data = None
count = 0
if len(files) == 0:
    raise ValueError("Couldn't find any files")
for f in files:
    print count, f
    sys.stdout.flush()
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

if 'chan_range' in pspec_cfg:
    cmin,cmax = map(lambda x: int(x.strip()), pspec_cfg['chan_range'].split("_"))
    chans = np.arange(cmin, cmax)
else:
    chans = np.arange(data.Nfreqs)

data.select(freq_chans=chans)
Bandwidth = data.freq_array[0, -1] - data.freq_array[0,0]
print 'Bandwidth (MHz): ' , Bandwidth/1e6

if opts.random_vis:
    shape = (data.Nblts, data.Nspws)
    for fi in range(data.Nfreqs):
        sigma = np.sqrt(np.var(data.data_array[:,:,fi,:]))
        mu = np.mean(data.data_array[:,:,fi,:])
        if opts.zero:
            if fi == 0: print("!!! Setting mean to 0")
            mu = 0.0
        if fi == 0: print 'random_vis: mu={}, sigma={}'.format(mu, sigma)
        imfac=1
#        if 'conj_y' in pspec_cfg:
#            imfac = -1
        vis_I = (np.random.normal(np.real(mu), sigma, size=shape) + imfac*(1j)*np.random.normal(np.imag(mu), sigma, size=shape))/np.sqrt(2.)
        data.data_array[:,:,fi,:] = np.repeat(vis_I[...,np.newaxis], data.Npols, axis=-1) #Same random data in all pols
#        data.data_array[:,:,fi,:] = (np.random.normal(np.real(mu), sigma, size=shape) + imfac*(1j)*np.random.normal(np.imag(mu), sigma, size=shape))/np.sqrt(2.)

#-------------------------------------------------------------------------------
# Get obs (and Opp) from first metadata
#-------------------------------------------------------------------------------
if 'bsq_int' in data.extra_keywords:
    Opp_I = data.extra_keywords['bsq_int'] * np.ones(data.Nfreqs)
elif 'opp' in data_cfg:
    Opp_I = data_cfg['opp'] * np.ones_like(data.Nfreqs)
else:
    if 'obs' in data_cfg.keys():
        obsfile = data_cfg['obs']
    else:
        #Take the first recognized obs file.
        obsfile = glob.glob( os.path.join(data_cfg['root'], 'metadata/*_obs.sav') )[0]
        obs = readsav(obsfile)['obs']
        Opp_I = obs['PRIMARY_BEAM_SQ_AREA'][0][0]

Opp_I = Opp_I[chans]
 
#Op = obs['PRIMARY_BEAM_AREA'][0]   # Shape = (npol, nfreq), type=float32. This is from FHD.

#n_freq_obs = obs['n_freq'][0]
#freq_center_obs = obs['freq_center'][0]
#freq_res_obs = obs['freq_res'][0]
#off = (n_freq_obs/2.)*freq_res_obs
#freq_obs = np.linspace( freq_center_obs-off, freq_center_obs+off, n_freq_obs)
#Bandwidth = freq_obs.size * freq_res_obs
#print("Obs bandwidth (MHz): ", Bandwidth/1e6)


#-------------------------------------------------------------------------------
# Delay transformation in stokes I
#-------------------------------------------------------------------------------


if data.Npols == 1:
    # eorsky, stokes I
    vis_I = data.data_array[:,0,:,0]
#    inds = np.where(vis_I == 0.0)
#    vis_I = np.delete(vis_I, inds[0], axis=0)
elif data.Npols == 2:
    # FHD, XX and YY
    vis_xx = data.data_array[:,0,:,0]
    vis_yy = data.data_array[:,0,:,1]
#    # Remove flagged visibilities
        ## Flagged visibilities should have been removed on the file read. This will throw off averaging if not careful.
#    inds = np.where(vis_xx == 0.0)
#    vis_xx = np.delete(vis_xx, inds[0], axis=0)
#    vis_yy = np.delete(vis_yy, inds[0], axis=0)
    if ('conj_y' in pspec_cfg) and not opts.random_vis:
        print("Conj_y line 191")
        vis_yy = np.conj(vis_yy)
    vis_I = vis_xx + vis_yy
#data.Ntimes = vis_I.shape[0]
#data.time_array = np.delete(data.time_array, inds=0)

if opts.zero and not opts.random_vis:
    print("!!! Setting per-freq means to 0 for data")
    vis_I = (vis_I.T - np.mean(vis_I, axis=1).T).T
print('Vis_I Mean: ', np.mean(vis_I))
#Convert to K*str
vis_I = vis_I * jy2Tstr(data.freq_array[0])

# FFT
dnu = np.diff(data.freq_array[0])[0]
_visI = np.fft.ifft(vis_I, axis=1)
print("dnu", dnu)
dspec_instr = _visI*_visI.conj()

# Get k_parallel
Zs = 1420e6/data.freq_array[0] - 1
Zmean = np.mean(Zs)
etas = np.fft.fftfreq(data.Nfreqs, d=data.channel_width)
k_parallel = dk_deta(Zmean) * etas

# Cosmological Normalization:
scalar = X2Y(Zmean) * (Bandwidth / Opp_I)
dspec_I =  dspec_instr * scalar


#-------------------------------------------------------------------------------
# Reference Line
#-------------------------------------------------------------------------------

if 'nside' in data.extra_keywords:
    if 'skysig' in data.extra_keywords:
        catalog_cfg['Nside'] = data.extra_keywords['nside']
        data.extra_keywords['nside'] = data.extra_keywords['skysig']
if 'Nside' in catalog_cfg:
    if 'sigma' in catalog_cfg:
        nside = catalog_cfg['Nside']
        sigma = catalog_cfg['sigma']
        om = 4*np.pi/(12.*nside**2)
        dnu = np.diff(data.freq_array[0])[0]   #Hz
        ref_line = X2Y(Zmean) * dnu * om * sigma**2
#        ref_line = comoving_voxel_volume(Zmean, dnu, om) * (sigma)**2
elif 'furlanetto_flat' in catalog_cfg:
    f = readsav('/users/alanman/FHD/catalog_data/eor_power_1d.idlsave')
    maxpow = max(f['power'])
    ref_line = maxpow/0.7**3   #Furlanetto's curve includes little-h
    ref_line *= 1e6     # mK^2 to K^2
elif 'npz' in catalog_cfg:
    # Setting reference level to the most averaged spectrum in an npz file.
    f = np.load(catalog_cfg['npz'])
    try:
        ref_avg_mode = f['avg_mode']
    except KeyError:
        ref_avg_mode = 'over_k'

    ref_avgs = f['avgs']
    if ref_avg_mode == 'over_k':
        print ref_avgs.shape
        ref_line = np.average(ref_avgs[-1])
    if ref_avg_mode == 'partitioned':
        print ref_avgs.shape
        ref_line = np.average(ref_avgs[-1, -1])
else:
    ref_line = np.average(dspec_I, axis=0)

try:
    save_dict['ref_line'] = ref_line
except NameError:
    pass


#-------------------------------------------------------------------------------
# Averaging pspecs
#-------------------------------------------------------------------------------

# Zero-visibilities have been removed and Ntimes/time array have been updated accordingly.
inttime = float(data.integration_time[0])
daytominutes = 24*60.
avg_lens_time = (data.time_array - data.time_array[0]) * daytominutes
avg_lens_time = avg_lens_time[2:]  #
avg_lens_nsamp = avg_lens_time/(inttime/60.)
avg_lens_nsamp = np.ceil(avg_lens_nsamp).astype(int)

# Don't need to calculate the variance for every single averaging length.
# Select a cadence for averaging lengths

calc_rate = int(4 / (inttime/60.)) # Minutes to integrations
avg_lens_nsamp = avg_lens_nsamp[::calc_rate]
avg_lens_time = avg_lens_time[::calc_rate]

print('Calc rate: ', calc_rate)

avgs = []
var_stds = []
variances = []
means = []
lens = avg_lens_nsamp
nsamps = []

def varfunc(arr):
    return np.var(arr)
    #return np.sqrt(np.mean( np.abs( arr - ref_line)**2))

if multislice:
        # Take steps of given lengths. Do the progressive averaging.
        # The goal is to find the sampling length at which there is a 1/N relationship
        variances = {}
        avgs = {} 
        step_lengths = np.arange(1, max_step)
        for s in step_lengths:
            dspec_slice = dspec_I[::s,:]
            avgs[s] = [np.mean(dspec_slice[0:l], axis=0) for l in lens]
            variances[s] = [varfunc(avgs[s][l]) for l in lens]
        ## TODO now save... figure out how to process.
else:
    for i,l in enumerate(lens):
        if partition:
            parts = np.array([dspec_I[j*l:(j+1)*l,:] for j in range(int(np.floor(data.Ntimes/float(l))))])
            parts = [ p for p in parts if len(p) > 1]
            avgs.append(np.array([np.mean(p, axis=0) for p in parts if len(p) > 0]))
            means.append(np.mean(avgs[-1][0]))
            vs = np.apply_along_axis(varfunc, 1, avgs[-1])
            if np.any(np.isnan(vs)):
                print(l,avgs[-1].shape)
#                import IPython; IPython.embed()
                raise ValueError("Variances are nan") 
            variances.append(np.mean(vs))
            var_stds.append(np.sqrt(np.var(vs)))
            nsamps.append(len(parts))
        elif over_k:
            avgs.append(np.mean(dspec_I[0:l,:], axis=0))
            variances.append(varfunc(avgs[-1]))
            var_stds.append(np.sqrt(np.var(variances[-1])))
            means.append(np.mean(avgs[-1]))
            nsamps.append(dspec_I.shape[1])
        elif boot:
            # Use bootstrapping method.
            r_inds = np.random.randint(0, data.Ntimes, size = Nboots)
            boots = [ dspec_I[np.arange(ri,ri+l)%(data.Ntimes),:] for ri in r_inds]  #Modulus should ensure circular indexing.
            inds = [ np.arange(ri, ri+l)%data.Ntimes for ri in r_inds]
            avgs.append(np.mean(boots, axis=1))
            vs = np.apply_along_axis(varfunc, 1, avgs[-1])
            variances.append(np.mean(vs))
            var_stds.append(np.sqrt(np.var(vs)))
            nsamps.append(Nboots)
            

#avgs_second = [np.mean(dspec_I[0:l,:], axis=0) for l in xrange(dspec_I.shape[0])]
print("Mean pspec_instr: {:.4e}".format(np.mean(dspec_instr)))
print("Scalar=", np.mean(scalar))
print("Scalar*mean_pspc_instr=", np.mean(scalar)*np.mean(dspec_instr))
print "Avg/Ref: ", np.mean(avgs[-1])/ref_line
print 'Integration time = {}'.format(data.integration_time[0])
print 'Compare: {} samples vs. {} = length of dspec_I'.format(lens[-1], dspec_I.shape[0])

save_dict.update(dict(k_parallel=k_parallel, avgs=avgs, avg_lens_nsamp=lens,
                      avg_lens_time=avg_lens_time, dspec_instr=dspec_instr,
                      scalar=scalar, bl=bls, nsamps=nsamps, variances=variances,
                      var_stds=var_stds, means=means))

save_results(save_dict, data_cfg, out_cfg)
