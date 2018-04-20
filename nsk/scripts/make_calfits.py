from pyuvdata import UVCal
import numpy as np

def make_calfits(fname, data_array, freq_array, time_array, jones_array, ants, flag_array=None,
                 channel_width=0.0, gain_convention='multiply', history='', telescope_name='HERA',
                 x_orientation='east', integration_time=10.0, freq_range=None, clobber=False):
    """
    make a calfits file from data_array etc. objects
   
    fname : str, filename

    data_array : ndarray, shape=(Nants, Nfreqs, Ntimes, Npols)

    freq_array : ndarray, shape=Nfreqs

    time_array : ndarray, shape=Ntimes

    jones_array : ndarray, shape=Npols

    ants : ndarray, shape=Nants
    """
    # specify ant params
    Nants_data = len(ants)
    Nants_telescope = len(ants)
    ant_array = np.array(ants, np.int)
    antenna_names = np.array(ants, np.str)
    antenna_numbers = np.array(ants, np.int)

    # frequency params
    Nfreqs = len(freq_array)
    if freq_range is None:
        freq_range = np.array([freq_array[0], freq_array[-1]])
    freq_array = freq_array.reshape(1, -1) 

    # pol params
    Njones = len(jones_array)

    # time params
    Ntimes = len(time_array)
    time_range = np.array([time_array[0], time_array[-1]])

    # spw params
    Nspws = 1
    spw_array = np.array([0])
    data_array = data_array[:, np.newaxis, :, :, :]
    if flag_array is not None:
        flag_array = flag_array[:, np.newaxis, :, :, :].astype(np.bool)

    # data params
    if data_array.shape[2] > 1:
        gain_array = data_array
        delay_array = None
        cal_type = 'gain'
    else:
        gain_array = None
        delay_array = data_array
        cal_type = 'delay'

    if flag_array is None:
        flag_array = np.array(np.zeros_like(data_array), np.bool)
    quality_array = np.zeros_like(data_array, np.float64) 

    # make blank uvc
    uvc = UVCal()
    params = ['Nants_data', 'Nants_telescope', 'ant_array', 'antenna_names', 'antenna_numbers', 'Nfreqs',
              'freq_array', 'Njones', 'Ntimes', 'time_range', 'Nspws', 'spw_array', 'data_array', 'gain_array',
              'delay_array', 'jones_array', 'time_array', 'freq_array', 'cal_type', 'flag_array', 'quality_array', 'channel_width', 'gain_convention',
              'history', 'telescope_name', 'x_orientation', 'integration_time', 'freq_range']
    for p in params:
        uvc.__setattr__(p, locals()[p])

    uvc.set_redundant()
    uvc.write_calfits(fname, clobber=clobber)

