'''This is a calibration file for data collected at PAPER in Boolardy
on JD 2454303.'''

import aipy as a, numpy as n, glob

class AntennaArray(a.fit.AntennaArray):
    def __init__(self, *args, **kwargs):
        self.sh_bp_poly = [-232718.7113099774, 143295.36518061417, -32767.774150412501, 3290.3465007911891, -121.51333810472337]
        self.blpolys = {
            a.miriad.ij2bl(0,1): [-605821024356.96509, 929327909307.20325, -639459462462.92554, 259893778099.61926, -69088708767.947281, 12551566807.948889, -1578126290.8351059, 135589226.64294812, -7618297.2864271775, 252763.07267842517, -3760.4444961469708] ,
            a.miriad.ij2bl(1,2): [-363950777552.90387, 551299949563.00793, -374450804986.1629, 150166447205.28934, -39373210391.942276, 7052081956.5478392, -873738998.22397363, 73937559.130901143, -4089407.1546051479, 133482.11906468633, -1952.4558572234528] ,
            a.miriad.ij2bl(0,3): [-314400863696.44781, 475591742451.27417, -322626614858.6228, 129239594502.42848, -33853433151.757805, 6058476863.01019, -750134823.4077388, 63445823.632708505, -3507908.504287676, 114480.33574243639, -1674.4753298044598] ,
            a.miriad.ij2bl(0,2): [-320523940797.76831, 490784507932.59937, -337114830020.28723, 136786033418.61569, -36305448645.863884, 6585943658.0970211, -826896530.1481384, 70950328.350660875, -3981373.5992053226, 131933.53181877444, -1960.4797399176532] ,
            a.miriad.ij2bl(2,3): [-549824022273.61963, 844460153112.9895, -581769954209.30811, 236732885152.59372, -63007167261.236298, 11460311654.46904, -1442617910.4162033, 124091388.95223087, -6980358.4309154879, 231864.4113735327, -3453.4937543959531] ,
            a.miriad.ij2bl(1,3): [-437464873820.13055, 662973515098.25513, -450534294799.39746, 180780502542.15082, -47429654063.030296, 8500930209.4925852, -1054060540.6233016, 89273844.940607563, -4942442.8821271891, 161502.8750442483, -2365.2396877128194] ,
        }
        a.fit.AntennaArray.__init__(self, *args, **kwargs)
        self._update_passband()
    def _update_passband(self):
        afreqs = self.ants[0].beam.afreqs
        self._passband = 1
        bm_resp = self[0].bm_response((0,0,1))
        self._beamgain = (bm_resp * n.conjugate(bm_resp)).squeeze()
        self._passband *= self._beamgain
        self._passband *= n.polyval(self.sh_bp_poly, afreqs)**2
        self._passband = n.clip(self._passband, 0, n.Inf)
        self.bl_gain = {}
        for bl in self.blpolys:
            self.bl_gain[bl] = n.polyval(self.blpolys[bl], afreqs)
            self.bl_gain[bl] = n.clip(self.bl_gain[bl], 0, n.Inf)
    def update(self):
        a.fit.AntennaArray.update(self)
        self._update_passband()
    def passband(self, i, j):
        bl = self.ij2bl(i,j)
        g = 1
        #g = self.bl_gain[bl]**2
        return a.fit.AntennaArray.passband(self, i, j) * self._passband * g
    def bm_response(self, i, j, pol='xx'):
        return a.fit.AntennaArray.bm_response(self,i,j,pol=pol)/self._beamgain


class BeamNoFlaps(a.ant.Beam):
    """A specific beam model for the PAPER experiment.  This model is for
    a single dipole element with no flaps."""
    def __init__(self, freqs, **kwargs):
        """The axes of the Cs matrices are polynomials in frequency going
        to the right, and polys in cos(2*az) going down."""
        #self.CsAm = n.array([        # N. Gugliucci 08/07
        #    [ 2.3541326  ,-0.0039133512 , 1.6055088e-05,-2.7468911e-08],
        #    [-0.46909345 , 0.0084471178 ,-5.1260711e-05, 1.0299793e-07],
        #    [ 0.32753617 ,-0.0081176326 , 5.8822952e-05,-1.3303273e-07],
        #    [-0.046844105, 0.00076223627,-3.5474502e-06, 4.4132035e-09],
        #    [-0.073523813, 0.0018151892 ,-1.3435102e-05, 3.1225928e-08],
        #    [ 0.047340855,-0.00085097424, 4.9799933e-06,-9.5164123e-09],
        #])
        #self.CsXo = n.array([        # N. Gugliucci 08/07
        #   [-121.29224, 1.9851554 ,-0.011876889  , 2.5222526e-05],
        #   [ 76.969303,-1.3947796 , 0.0085644354 ,-1.7448153e-05],
        #   [-36.638691, 0.93699466,-0.0068616164 , 1.5544311e-05],
        #   [ 10.189859,-0.18212180, 0.00098309486,-1.6152395e-06],
        #   [ 5.9997050,-0.15737420, 0.0012090764 ,-2.8862905e-06],
        #   [-5.6561847, 0.10468756,-0.00063126068, 1.2444705e-06],
        #])
        self.CsSd = n.array([        # N. Gugliucci 08/07
            [ 143.84525,-1.1088605  , 0.0048397670 ,-7.1054741e-06],
            [-104.00886, 1.9980993  ,-0.013304344  , 2.8955473e-05],
            [ 28.304230,-0.75088201 , 0.0056338561 ,-1.2898564e-05],
            [-8.7149717, 0.16108215 ,-0.00090283393, 1.5386691e-06],
            [-3.4672940, 0.091929321,-0.00071502397, 1.7311496e-06],
            [ 3.4123240,-0.063083812, 0.00038093617,-7.5356570e-07],
        ])
        a.ant.Beam.__init__(self, freqs)
        self._update_coeffs()
    def _update_coeffs(self):
        mhz_freqs = 1e3 * self.afreqs # GHz -> MHz
        mhz_freqs = n.array([
           n.ones_like(mhz_freqs),
           mhz_freqs,
           mhz_freqs**2,
           mhz_freqs**3,
        ])
        #self.BAm = n.dot(self.CsAm, mhz_freqs)
        #self.BXo = n.dot(self.CsXo, mhz_freqs)
        self.BSd = n.dot(self.CsSd, mhz_freqs)
    def select_chans(self, active_chans):
        a.ant.Beam.select_chans(self, active_chans)
        self._update_coeffs()
    def response(self, xyz):
        """Return the beam response across the active band for the specified
        topocentric coordinates (with z = up, x = east). 2nd axis should be 
        multiple coordinates.  Returns 'x' pol (rotate pi/2 for 'y')."""
        az,alt = a.coord.top2azalt(xyz)
        zang = n.pi/2 - alt
        if zang.size == 1:
            zang = n.array([zang]); zang.shape = (1,)
            az = n.array([az]); az.shape = (1,)
        A = n.array([0,2,4,6,8,10],dtype=n.float)
        A.shape = (1,) + A.shape; az.shape += (1,); zang.shape += (1,)
        A = n.cos(n.dot(az, A))
        A[:,0] = 0.5
        #a1 = n.dot(A, self.BAm)
        #a2 = n.dot(A, self.BXo)
        a3 = n.dot(A, self.BSd)
        #z = (180*zang/n.pi - a2) / a3
        z = (180*zang/n.pi) / a3
        #rv = n.sqrt(a1 * n.exp(-z**2/2))
        rv = n.sqrt(n.exp(-z**2/2))
        return rv.transpose()
    def get_params(self, prm_list):
        return {}
    def set_params(self, prm_list):
        return

class Antenna(a.fit.Antenna):
    def bm_response(self, top, pol='x'):
        """Return response of beam for specified polarization."""
        top = n.array(top)
        top = {'y':top, 'x':n.dot(self.rot_pol_y, top)}[pol]
        x,y,z = top
        return self.beam.response((x,y,z))

prms = {
    'loc': ('-26:44:12.74', '116:39:59.33'), # Boolardy, Australia
    'antpos':
        #[[   0.,    0.,    0.],
        # [-105.39, 139.44,-201.57],
        # [ -68.28, 381.97,-134.25],
        # [  62.76, 462.75, 119.78]],
        #---------------------------
        #[[  -2.59,  -0.07,  -0.54],
        # [-103.07, 138.19,-198.50],
        # [ -71.14, 381.77,-134.54],
        # [  56.90, 463.28, 117.57]],
        #---------------------------
        [[ -2.357, -0.088, -0.322],
         [  -103.138, 138.251, -198.001],
         [  -70.618, 381.601, -134.217],
         [  57.521, 462.655, 117.856]],

    'delays': [1.156, 0.503,  0.936, 1.058],
    'offsets': [0.003, 0.003,-0.062,-0.080],
    #'delays': [2.440, 1.311,  2.103, 2.245],
    #'offsets': [0.003, 0.073, 0.005, 0.003],
    'amps': [0.0176, 0.0166, 0.0163, 0.0183],
    #'amps': [1] * 4,
    #'bp_r':n.array([[-232718.7113099774, 143295.36518061417, -32767.774150412501, 3290.3465007911891, -121.51333810472337]]*4),
    'bp_r': n.array([
        [-9.4993629754934084, 2.4188473808454298, 0.85794454024305478],
        [-8.9490935732613472, 1.8780901359059294, 0.92187063301210914],
        [-4.0490730882252954, 1.3148666879751385, 0.88892127101152751],
        [-0.56177062767531205, 1.8241785585891628, 0.74296843853565719],
    ]),
    #'bp_r': n.array([[0,0,1]]*4),
    'bp_i': n.array([[0.000]]*4),
    'beam': BeamNoFlaps,
}

class _PeggedBody:
    def __init__(self):
        self._pegs = {}
    def peg(self, times, amps, inds, ras, decs):
        for t,a,i,r,d in zip(times,amps,inds,ras,decs):
            self._pegs[t] = (a,i,r,d)
    def fix_str(self, observer):
        t = str(observer.get_jultime())
        if self._pegs.has_key(t):
            self.janskies *= self._pegs[t][0]
            self.index += self._pegs[t][1]
    def fix_pos(self, observer):
        t = str(observer.get_jultime())
        if self._pegs.has_key(t):
            self._def_ionref = self.ionref
            self.ionref = [self._pegs[t][2], self._pegs[t][3]]
        else:
            try: self.ionref = self._def_ionref
            except(AttributeError): pass

class RadioFixedBody(a.fit.RadioFixedBody, _PeggedBody):
    def __init__(self, *args, **kwargs):
        a.fit.RadioFixedBody.__init__(self, *args, **kwargs)
        _PeggedBody.__init__(self)
    def compute(self, observer):
        _PeggedBody.fix_pos(self, observer)
        a.fit.RadioFixedBody.compute(self, observer)
        _PeggedBody.fix_str(self, observer)

class RadioSpecial(a.fit.RadioSpecial, _PeggedBody):
    def __init__(self, *args, **kwargs):
        a.fit.RadioSpecial.__init__(self, *args, **kwargs)
        _PeggedBody.__init__(self)
    def compute(self, observer):
        _PeggedBody.fix_pos(self, observer)
        a.fit.RadioSpecial.compute(self, observer)
        _PeggedBody.fix_str(self, observer)

pegged_srcs = {
    'cyg': glob.glob('24*.cyg'),
    #'cas': glob.glob('24*.cas'),
    #'vir': glob.glob('24*.vir'),
    #'crab': glob.glob('24*.crab'),
    'Sun': glob.glob('24*.Sun'),
}

def get_aa(freqs):
    '''Return the AntennaArray to be used fro simulation.'''
    beam = prms['beam'](freqs)
    try: beam.set_params(prms)
    except(AttributeError): pass
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    assert(len(prms['delays']) == nants and len(prms['amps']) == nants \
        and len(prms['bp_r']) == nants and len(prms['bp_i']) == nants)
    for pos, dly, off, amp, bp_r, bp_i in zip(prms['antpos'], prms['delays'],
            prms['offsets'], prms['amps'], prms['bp_r'], prms['bp_i']):
        antennas.append(
            a.fit.Antenna(pos[0],pos[1],pos[2], beam, phsoff=[dly,off],
                amp=amp, bp_r=bp_r, bp_i=bp_i)
        )
    aa = AntennaArray(prms['loc'], antennas)
    return aa

def get_catalog(srcs=None, cutoff=None):
    '''Return a catalog containing the listed sources.'''
    cat = a.src.get_catalog(srcs=[s for s in srcs if not s.startswith('J')],
        cutoff=cutoff, fixedbody=RadioFixedBody, special=RadioSpecial)
    for src in [s for s in srcs if s.startswith('J')]:
        cat[src] = a.fit.RadioFixedBody(0.,0.,janskies=0., mfreq=.15, name=src)
    cat.set_params(src_prms)
    for src in cat:
        if src in pegged_srcs and len(pegged_srcs[src]) > 0:
            dat = []
            for filename in pegged_srcs[src]:
                dat += [L.split() for L in open(filename).readlines()]
            times = [L[0] for L in dat]
            amps = [float(L[1]) for L in dat]
            inds = [float(L[2]) for L in dat]
            ras = [float(L[3]) for L in dat]
            decs = [float(L[4]) for L in dat]
            cat[src].peg(times, amps, inds, ras, decs)
    return cat

src_prms = {
        'Sun': {
            'a1': 0.00546, 'a2': 0.00587, 'th': -1.723,
            'str': 36878, 'index':2.616,
        },
        'vir': {
            'str': 832, 'index':-0.212,
        },
        '3c353': {
            'str': 180, 'index': 1.08,
        },
        '3c392': {
            'str': 480, 'index': 1.46,
            'a1':.00425, 'a2':.00482, 'th':-0.002
        },
        'hyd': {
            'str': 280, 'index': 0.40,
        },
        'sgr': {
            'str': 122, 'index':-4.21,
        },
        'crab': {
            'str': 1049, 'index':-0.85,
        },
        'pic': {
            'str': 200.,
            'index': -0.76,
        },
        'for': {
            'str': 144.,
            'index': -2.45,
        },
        'cen': {
            'a1': 0.00, 'a2': 0.00407, 'th': 0.637,
            'str': 1011., 'index':0.752,
        },
        #'J1660': {
        #    'ra': 4.49985380753,
        #    'dec': -1.15258589749,
        #    'str':122,
        #    'index':0,
        #},
        #'J2169' : {
        #    'ra':5.44748143695,
        #    'dec':-1.25357396202,
        #    'str':60,
        #    'index':0,
        #},
        #'J2216' : {
        #    'ra':'22:14',
        #    'dec':'-16:49', 
        #    'str':200, 
        #    'index':0,
        #},
        'J1360' : {
            'ra':'13:45:30.3',
            'dec':'-60:11:38',
            'str':190, 
            'index':0.57,
        },
        'J0842' : {
            'ra':'08:22:43.1',
            'dec':'-42:59:41.5',
            'str':118, 
            'index':-2.90,
        },
}
