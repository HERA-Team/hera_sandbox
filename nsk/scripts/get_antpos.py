import pyuvdata.uvutils as uvutils

def get_antpos(uvd):
    return uvutils.ENU_from_ECEF((uvd.antenna_positions + uvd.telescope_location).T, *uvd.telescope_location_lat_lon_alt).T


