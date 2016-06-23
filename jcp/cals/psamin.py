"""
Calfile for data taken from 2456240. - 2456378.
"""
import aipy as a, numpy as n

class AntennaArray(a.pol.AntennaArray):
    def __init__(self, *args, **kwargs):
        a.pol.AntennaArray.__init__(self, *args, **kwargs)
        self.ant_layout = kwargs.pop('ant_layout')
        self.amp_coeffs = kwargs.pop('amp_coeffs')
        self.gain = kwargs.pop('gain')
        self.update_gains()
        self.dly_coeffs = kwargs.pop('dly_coeffs')
        self.dly_xx_to_yy = kwargs.pop('dly_xx_to_yy')
        self.tau_ew = kwargs.pop('tau_ew')
        self.tau_ns = kwargs.pop('tau_ns')
        self.update_delays()
    def update_gains(self):
        gains = self.gain * self.amp_coeffs
        for i,gain in zip(self.ant_layout.flatten(), gains.flatten()):
            self[i].set_params({'amp_x':gain})
            self[i].set_params({'amp_y':gain})
    def update_delays(self):
        ns,ew = n.indices(self.ant_layout.shape)
        dlys = ns*self.tau_ns + ew*self.tau_ew + self.dly_coeffs
        dlys_xx_to_yy = ns*self.tau_ns + ew*self.tau_ew + self.dly_xx_to_yy
        for i,tau_x,tau_y in zip(self.ant_layout.flatten(), dlys.flatten(), dlys_xx_to_yy.flatten()):
            self[i].set_params({'dly_x':tau_x})
            self[i].set_params({'dly_y':tau_y})
    def update(self):
        self.update_gains()
        self.update_delays()
        a.pol.AntennaArray.update(self)
    def get_params(self, ant_prms={'*':'*'}):
        try: prms = a.pol.AntennaArray.get_params(self, ant_prms)
        except(IndexError): return {}
        for k in ant_prms:
            if k == 'aa':
                if not prms.has_key('aa'): prms['aa'] = {}
                for val in ant_prms[k]:
                    if   val == 'tau_ns': prms['aa']['tau_ns'] = self.tau_ns
                    elif val == 'tau_ew': prms['aa']['tau_ew'] = self.tau_ew
                    elif val == 'gain': prms['aa']['gain'] = self.gain
            else:
                try: top_pos = n.dot(self._eq2zen, self[int(k)].pos)
                # XXX should multiply this by len_ns to match set_params.
                except(ValueError): continue
                if ant_prms[k] == '*':
                    prms[k].update({'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]})
                else:
                    for val in ant_prms[k]:
                        if   val == 'top_x': prms[k]['top_x'] = top_pos[0]
                        elif val == 'top_y': prms[k]['top_y'] = top_pos[1]
                        elif val == 'top_z': prms[k]['top_z'] = top_pos[2]
        return prms
    def set_params(self, prms):
        changed = a.pol.AntennaArray.set_params(self, prms)
        for i, ant in enumerate(self):
            ant_changed = False
            top_pos = n.dot(self._eq2zen, ant.pos)
            try:
                top_pos[0] = prms[str(i)]['top_x']
                ant_changed = True
            except(KeyError): pass
            try:
                top_pos[1] = prms[str(i)]['top_y']
                ant_changed = True
            except(KeyError): pass
            try:
                top_pos[2] = prms[str(i)]['top_z']
                ant_changed = True
            except(KeyError): pass
            if ant_changed: ant.pos = n.dot(n.linalg.inv(self._eq2zen), top_pos) / a.const.len_ns
            changed |= ant_changed
        try: self.tau_ns, changed = prms['aa']['tau_ns'], 1
        except(KeyError): pass
        try: self.tau_ew, changed = prms['aa']['tau_ew'], 1
        except(KeyError): pass
        try: self.gain, changed = prms['aa']['gain'], 1
        except(KeyError): pass
        if changed: self.update()
        return changed

prms = {
    'loc': ('-30:43:17.5', '21:25:41.9'), # KAT, SA (GPS)
    'antpos': {
        0:  {'top_x':  600.0,  'top_y':    0.00,  'top_z':    0.0},
        1:  {'top_x': 1000.0,  'top_y':    0.00,  'top_z':    0.0},
        2:  {'top_x': 1400.0,  'top_y':    0.00,  'top_z':    0.0},
        3:  {'top_x': 1800.0,  'top_y':    0.00,  'top_z':    0.0},
        
        4:  {'top_x':  400.0,  'top_y':  346.41,  'top_z':    0.0},
        5:  {'top_x':  800.0,  'top_y':  346.41,  'top_z':    0.0},
        6:  {'top_x': 1200.0,  'top_y':  346.41,  'top_z':    0.0},
        7:  {'top_x': 1600.0,  'top_y':  346.41,  'top_z':    0.0},
        8:  {'top_x': 2000.0,  'top_y':  346.41,  'top_z':    0.0},
        
        9:  {'top_x':  200.0,  'top_y':  692.82,  'top_z':    0.0},
        10: {'top_x':  600.0,  'top_y':  692.82,  'top_z':    0.0},
        11: {'top_x': 1000.0,  'top_y':  692.82,  'top_z':    0.0},
        12: {'top_x': 1400.0,  'top_y':  692.82,  'top_z':    0.0},
        13: {'top_x': 1800.0,  'top_y':  692.82,  'top_z':    0.0},
        14: {'top_x': 2200.0,  'top_y':  692.82,  'top_z':    0.0},
        
        15: {'top_x':    0.0,  'top_y': 1039.23,  'top_z':    0.0},
        16: {'top_x':  400.0,  'top_y': 1039.23,  'top_z':    0.0},
        17: {'top_x':  800.0,  'top_y': 1039.23,  'top_z':    0.0},
        18: {'top_x': 1600.0,  'top_y': 1039.23,  'top_z':    0.0},
        19: {'top_x': 2000.0,  'top_y': 1039.23,  'top_z':    0.0},
        20: {'top_x': 2400.0,  'top_y': 1039.23,  'top_z':    0.0},
        
        21: {'top_x':  200.0,  'top_y': 1385.64,  'top_z':    0.0},
        22: {'top_x':  600.0,  'top_y': 1385.64,  'top_z':    0.0},
        23: {'top_x': 1000.0,  'top_y': 1385.64,  'top_z':    0.0},
        24: {'top_x': 1400.0,  'top_y': 1385.64,  'top_z':    0.0},
        25: {'top_x': 1800.0,  'top_y': 1385.64,  'top_z':    0.0},
        26: {'top_x': 2200.0,  'top_y': 1385.64,  'top_z':    0.0},
        
        27: {'top_x':  400.0,  'top_y': 1732.05,  'top_z':    0.0},
        28: {'top_x':  800.0,  'top_y': 1732.05,  'top_z':    0.0},
        29: {'top_x': 1200.0,  'top_y': 1732.05,  'top_z':    0.0},
        30: {'top_x': 1600.0,  'top_y': 1732.05,  'top_z':    0.0},
        31: {'top_x': 2000.0,  'top_y': 1732.05,  'top_z':    0.0},
        
        32: {'top_x':  600.0,  'top_y': 2078.46,  'top_z':    0.0},
        33: {'top_x': 1000.0,  'top_y': 2078.46,  'top_z':    0.0},
        34: {'top_x': 1400.0,  'top_y': 2078.46,  'top_z':    0.0},
        35: {'top_x': 1800.0,  'top_y': 2078.46,  'top_z':    0.0},
    
    #to be excluded
        36: {'top_x': 12000.0,  'top_y': 10039.23,  'top_z':    0.0},
#        37: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
#        38: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
#        39: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
#        40: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
#        41: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
#        42: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
#        43: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
#        44: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
#        45: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
#        46: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
#        47: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
#        48: {'top_x': 1200.0,  'top_y': 1039.23,  'top_z':    0.0},
    },
    'ant_layout': n.array(
        [[ 32, 33, 34, 35, 36, 36, 36],
         [ 27, 28, 29, 30, 31, 36, 36],
         [ 21, 22, 23, 24, 25, 26, 36],
         [ 15, 16, 17, 36, 18, 19, 20],
         [ 36,  9, 10, 11, 12, 13, 14],
         [ 36, 36,  4,  5,  6,  7,  8],
         [ 36, 36, 36,  0,  1,  2,  3]]
         ),
    'dly_coeffs': n.array(
       [[ 0., 0., 0., 0., 0., 0., 0. ],
        [ 0., 0., 0., 0., 0., 0., 0. ],
        [ 0., 0., 0., 0., 0., 0., 0. ],
        [ 0., 0., 0., 0., 0., 0., 0. ],
        [ 0., 0., 0., 0., 0., 0., 0. ],
        [ 0., 0., 0., 0., 0., 0., 0. ],
        [ 0., 0., 0., 0., 0., 0., 0. ]]
       ),
    #dly_xx_to_yy is actually the yy delays and not an offset from xx. misnomer.
    'dly_xx_to_yy': n.array(
        [[ 0., 0., 0., 0., 0., 0., 0. ],
         [ 0., 0., 0., 0., 0., 0., 0. ],
         [ 0., 0., 0., 0., 0., 0., 0. ],
         [ 0., 0., 0., 0., 0., 0., 0. ],
         [ 0., 0., 0., 0., 0., 0., 0. ],
         [ 0., 0., 0., 0., 0., 0., 0. ],
         [ 0., 0., 0., 0., 0., 0., 0. ]]
         ), # PSA6240-PSA6378 one 11/5/2013

    'tau_ew': 0.00,

    'tau_ns': 0.00,
    #'delays': {},
    'gain': 1.,
    'amp_coeffs': n.array(
        [[ 1.   , 1.   , 1.   , 1.   , 1.   , 1.   , 1.   ],
         [ 1.   , 1.   , 1.   , 1.   , 1.   , 1.   , 1.   ],
         [ 1.   , 1.   , 1.   , 1.   , 1.   , 1.   , 1.   ],
         [ 1.   , 1.   , 1.   , 1.   , 1.   , 1.   , 1.   ],
         [ 1.   , 1.   , 1.   , 1.   , 1.   , 1.   , 1.   ],
         [ 1.   , 1.   , 1.   , 1.   , 1.   , 1.   , 1.   ],
         [ 1.   , 1.   , 1.   , 1.   , 1.   , 1.   , 1.   ]]
         ), # psa6240.17, xx
    'amps': {},
    'bp_r': n.array([[1.]] * 49), # from pic (danny's paper).
    'bp_i': n.array([[0., 0., 0.]] * 49),
    'beam': a.fit.BeamAlm,
   'twist': n.array([0]*49),#n.array([.1746] * 32),
    'bm_prms': {}
}
#
def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    for i in range(nants):
        beam = prms['beam'](freqs, nside=49, lmax=20, mmax=20, deg=7)
        try: beam.set_params(prms['bm_prms'])
        except(AttributeError): pass
        phsoff = {'x':[0.,0.], 'y':[0.,0.]}
        amp = prms['amps'].get(i, 4e-3); amp = {'x':amp,'y':amp}
        bp_r = prms['bp_r'][i]; bp_r = {'x':bp_r, 'y':bp_r}
        bp_i = prms['bp_i'][i]; bp_i = {'x':bp_i, 'y':bp_i}
        twist = prms['twist'][i]
        antennas.append(a.pol.Antenna(0., 0., 0., beam, phsoff=phsoff,
                amp=amp, bp_r=bp_r, bp_i=bp_i, pointing=(0.,n.pi/2,twist)))
    aa = AntennaArray(prms['loc'], antennas, tau_ew=prms['tau_ew'], tau_ns=prms['tau_ns'],
        gain=prms['gain'], amp_coeffs=prms['amp_coeffs'],
        dly_coeffs=prms['dly_coeffs'], dly_xx_to_yy=prms['dly_xx_to_yy'], ant_layout=prms['ant_layout'])
    for i in range(nants):
        pos = prms['antpos'][i]
        i = str(i)
        aa.set_params({i:pos})
    return aa


