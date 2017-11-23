from pyuvdata import uvutils
import numpy as np

def get_antpos(uvd, center=True, pick_data_ants=False):
    """
    get antenna positions in topocentric coordinates in units of meters

    uvd : pyuvdata.UVData object
    center : bool, if True, subtract median of array position
    """
    antpos = uvutils.ENU_from_ECEF((uvd.antenna_positions + uvd.telescope_location).T, *uvd.telescope_location_lat_lon_alt).T
    ants = uvd.antenna_numbers

    if pick_data_ants:
        data_ants = np.unique(np.concatenate([uvd.ant_1_array, uvd.ant_2_array]))
        telescope_ants = uvd.antenna_numbers
        select = map(lambda x: x in data_ants, telescope_ants)
        antpos = antpos[select, :]
        ants = telescope_ants[select]

    if center is True:
        antpos -= np.median(antpos, axis=0)

    return antpos, ants

