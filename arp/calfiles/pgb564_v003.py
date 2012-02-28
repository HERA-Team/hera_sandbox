'''This is a calibration file for data collected at PAPER in Green Bank
on JD 2454564.'''

import aipy as a, numpy as n, pyfits

class BeamFlaps(a.fit.BeamFlat):
    def __init__(self, freqs, fitsfile):
        a.fit.BeamFlat.__init__(self, freqs)
        self.from_fits(filename=fitsfile)
    def from_fits(self, filename, hdunum=1):
        # Figure out how many cols there are (is there a better way?)
        hdu = pyfits.open(filename)[hdunum]
        ncol = 1
        while True:
            try: i = hdu.data.field(ncol)
            except(IndexError): break
            ncol += 1
        nside = hdu.header['NSIDE']
        del(hdu)
        # Make a spectral index HPM for each additional col in fits file
        self.hmap = []
        for i in range(ncol):
            h = a.healpix.HealpixMap(nside=nside)
            h.from_fits(filename, hdunum=hdunum, colnum=i)
            self.hmap.append(h)
    def response(self, top):
        """Return beam response across active band for specified topocentric 
        coordinates (x=E,y=N,z=UP). x,y,z may be multiple coordinates.  
        Returns 'x' pol (rotate pi/2 for 'y')."""
        top = [a.healpix.mk_arr(c, dtype=n.double) for c in top]
        px,wgts = self.hmap[0].crd2px(*top, **{'interpolate':1})
        poly = n.array([n.sum(h.map[px] * wgts, axis=-1) for h in self.hmap])
        rv = n.polyval(poly, n.reshape(self.afreqs, (self.afreqs.size, 1)))
        return rv

prms = {
    'loc': ('38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
    'antpos':
        [[   0.06,    0.05,    0.03],
         [ 213.80, -136.50, -260.92],
         [ 196.48, -809.82, -253.66],
         [-254.86, -673.84,  308.06],
         [-284.72, -462.60,  350.35],
         [-277.44, -361.49,  342.45],
         [-174.90, -102.47,  218.02],
         [ -75.13,  -21.22,   96.93]],
    'delays': [ 0.37,-10.08,  5.91,  6.51,  4.50,  1.92,  2.18,  0.09],
    'offsets': [
        [0,0,0,0,0.0], 
        [0,0,0,0,0.5], 
        [0,0,0,0,0.501], 
        [0,0,0,0,0.5], 
        [0,0,0,0,0.5],
        [0,0,0,0,0.999], 
        [0,0,0,0,0.499], 
        [0,0,0,0,0.501],
    ],
    'amps': [ .00245, .00239, .00272, .00284, .00271, .00211, .00251, .00203],
    'bp_r': n.array([
        [ -1.0217e+06,  6.2546e+05, -1.4296e+05,  1.4452e+04, -5.4406e+02],
        [ -6.1083e+05,  3.7618e+05, -8.6584e+04,  8.8220e+03, -3.3460e+02],
        [ -8.1105e+05,  4.9770e+05, -1.1407e+05,  1.1567e+04, -4.3676e+02],
        [ -8.2055e+05,  5.0252e+05, -1.1503e+05,  1.1657e+04, -4.4007e+02],
        [ -5.9154e+05,  3.6374e+05, -8.3638e+04,  8.5167e+03, -3.2288e+02],
        [ -8.7310e+05,  5.3409e+05, -1.2209e+05,  1.2355e+04, -4.6579e+02],
        [ -5.3423e+05,  3.3100e+05, -7.6588e+04,  7.8379e+03, -2.9820e+02],
        [ -7.8334e+05,  4.7941e+05, -1.0958e+05,  1.1081e+04, -4.1703e+02],
    ]),
    'bp_i': n.array([
        [0.000],
        [0.000],
        [0.000],
        [0.000],
        [0.000],
        [0.000],
        [0.000],
        [0.000],
    ]),
    'beam': BeamFlaps,
}

def get_aa(freqs):
    '''Return the AntennaArray to be used fro simulation.'''
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    assert(len(prms['delays']) == nants and len(prms['amps']) == nants \
        and len(prms['bp_r']) == nants and len(prms['bp_i']) == nants)
    for i in range(len(prms['antpos'])):
        beam = prms['beam'](freqs, fitsfile='bmflaps.fits')
        pos = prms['antpos'][i]
        dly = prms['delays'][i]
        amp = prms['amps'][i]
        bp_r = prms['bp_r'][i]
        bp_i = prms['bp_i'][i]
        off = prms['offsets'][i]
        antennas.append(
            a.fit.Antenna(pos[0],pos[1],pos[2], beam, delay=dly, offset=off,
                amp=amp, bp_r=bp_r, bp_i=bp_i)
        )
    aa = a.fit.AntennaArray(prms['loc'], antennas)
    return aa

src_prms = {
        #'cyg': {
        #    'str': 10500,
        #    'index':-0.69
        #},
        #'cas': {
        #    #'str': 9150,
        #    #'str': 12200,
        #},
        'Sun': {
            'str': 53400,
            'index':2.00,
            'angsize':0.00896,
        },
        #'vir': {
        #    #'str': 1446,
        #    'str': 1700,
        #    #'index': -0.86,
        #    'index': -0.75
        #},
        #'crab': {
        #    #'str':  1838,
        #    'str': 1400,
        #    #'index': -0.30,
        #    'index': -0.29
        #},
}
