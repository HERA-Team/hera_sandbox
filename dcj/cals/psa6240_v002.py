"""
Cal File for 2012-2013 observing season.
"""
import aipy as a, numpy as n
from capo import EoR2012 as EOR 

prms = {
    'dly_coeffs': n.array(
        [[  0.  ,  0.  ,-10.11,  8.25,  5.64, -0.5 ,-18.21,-16.12],
         [  0.  ,  5.98,-13.64, -7.1 , -4.3 ,-16.03, -4.47,-16.3 ],
         [ 10.67, -5.63,  3.02,-10.09,  0.6 ,-15.3 ,-16.4 ,-19.23],
         [ 10.03,-11.03,-13.06,  2.34,  1.75,-10.77,-11.32,-25.84],
         [  5.81,-10.16, -6.85,  0.  ,-20.13,-13.41,-17.54,-24.47],
         [-12.94, -4.42, -9.29,-12.23,-18.93, -5.75, -4.13,-24.88],
         [  8.22,  0.45, -1.28, 16.97,-21.84,-22.09,-12.38,-26.46],
         [-12.81, -8.38, -5.57, -7.54,-20.48,-17.75,-25.85,-16.06]]), # PSA6240.17, xx
    'dly_xx_to_yy': n.array(
        [[  0.  , -9.55,-21.49, -1.61, -9.06,-12.11,-17.88,-19.54],
         [-12.5 , -1.94,-26.47, -5.45,-18.6 ,-16.18,-19.95,-17.17],
         [  0.82,-17.5 , -7.33,-24.22, -6.28,-25.92, -6.95, -6.33],
         [ -3.56,  0.76,-13.87, -5.02, -9.93,-22.25,-22.76,-35.87],
         [ -5.46,-11.33,-18.08,  0.  ,-24.2 ,-29.23,-18.17,-36.04],
         [-19.83,-15.44,-10.53,-26.86,-25.85,-15.9 ,-15.99,-21.4 ],
         [-20.9 , -9.89,-11.08,  6.37,-34.94,-27.21,-19.38,-38.45],
         [-22.48,-14.31,-16.6 ,-19.01,-32.48,-23.04,-23.94,-40.45]]), # psa6240.17, yy - xx
#    'tau_ew': 2.71,
#     'tau_ew': -4.84,#41_49 (0,1)
#     'tau_ew':-4.837, #1st row (0,1)
     'tau_ew':2.222, #good on 41_49,47_49 0.66, bls 0,4 like at 0.64
#    'tau_ew':-5.26, #preferrerd for (0,4)
#     'tau_ns': -3.56,#default
#     'tau_ns':-2.554,#49_58
     'tau_ns':0.944, #49_58 0.82
#      'tau_ns':-6.094, #preferred by 3_49
#    'gain': .0036096,
    'gain':0.0658,
    'amp_coeffs': n.array(
        [[ 1.   , 1.228, 1.127, 1.119, 1.114, 0.807, 1.019, 1.197],
         [ 1.051, 1.283, 1.159, 1.143, 1.102, 1.232, 0.773, 1.202],
         [ 1.18 , 1.083, 1.17 , 1.154, 1.18 , 1.119, 1.124, 1.104],
         [ 1.261, 1.105, 1.15 , 1.191, 1.157, 1.12 , 1.161, 1.156],
         [ 0.841, 1.003, 0.951, 1.   , 1.276, 1.092, 1.274, 0.976],
         [ 1.302, 0.982, 1.291, 0.979, 1.187, 0.993, 1.404, 0.824],
         [ 0.908, 0.959, 1.182, 0.653, 1.241, 0.958, 1.392, 0.878],
         [ 1.402, 0.866, 1.319, 0.967, 1.427, 0.985, 1.179, 0.858]]), # psa6240.17, xx
    'amps': {},
    'twist': n.array([0]*64),#n.array([.1746] * 32),
}

def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = EOR.prms['loc']
    antennas = []
    nants = len(EOR.prms['antpos'])
    for i in range(nants):
        beam = EOR.prms['beam'](freqs, nside=32, lmax=20, mmax=20, deg=7)
        try: beam.set_params(EOR.prms['bm_prms'])
        except(AttributeError): pass
        phsoff = {'x':[0.,0.], 'y':[0.,0.]}
        amp  = prms['amps'].get(i, 4e-3); amp = {'x':amp,'y':amp}
        bp_r = EOR.prms['bp_r'][i]; bp_r = {'x':bp_r, 'y':bp_r}
        bp_i = EOR.prms['bp_i'][i]; bp_i = {'x':bp_i, 'y':bp_i}
        twist = prms['twist'][i]
        antennas.append(a.pol.Antenna(0., 0., 0., beam, phsoff=phsoff,
                amp=amp, bp_r=bp_r, bp_i=bp_i, pointing=(0.,n.pi/2,twist)))
    aa = EOR.AntennaArray(EOR.prms['loc'], antennas, tau_ew=prms['tau_ew'], tau_ns=prms['tau_ns'],
        gain=prms['gain'], amp_coeffs=prms['amp_coeffs'],
        dly_coeffs=prms['dly_coeffs'], dly_xx_to_yy=prms['dly_xx_to_yy'], ant_layout=EOR.prms['ant_layout'])
    for i in range(nants):
        pos = EOR.prms['antpos'][i]
        i = str(i)
        aa.set_params({i:pos})
    return aa

src_prms = {
'cen':{ 'jys':10**3.282102, 'index':  0.235166 , },
'cyg':{ 'jys':10**3.566410, 'index':  -0.266315 , },
'hyd':{ 'jys':10**2.448816, 'index':  -0.866462 , },
#'pic':{ 'jys':10**2.714456, 'index':  -0.436361 , },
'pic':{'jys':450, 'index':-1.},
'for':{'jys':447,'index':-1.15},
'vir':{ 'jys':10**2.200725, 'index':  0.202425 , },
'Sun': {'a1': 0.00644, 'index': 1.471, 'a2': 0.00586, 'jys': 55701.96, 'th': -0.000512},
#'for': {'a1': 0.00851, 'a2': 0.00413, 'jys': 907.09, 'th': 0.230},
}

def get_catalog(srcs=None, cutoff=None, catalogs=['helm','misc']):
    '''Return a catalog containing the listed sources.'''
    custom_srcs = ['J1347-603','J1615-610', 'J1336-340', 'J1248-412', 'J1531-423', 'J1359-415']
    if srcs is None:
        cat = a.src.get_catalog(srcs=srcs, cutoff=cutoff, catalogs=catalogs)
    else:
        cat = a.src.get_catalog(srcs=[s for s in srcs if not s in custom_srcs],
            cutoff=cutoff, catalogs=catalogs)
        for src in [s for s in srcs if s in custom_srcs]:
            cat[src] = a.fit.RadioFixedBody(0., 0., janskies=0., mfreq=.15, name=src)
    cat.set_params(src_prms)
    return cat

