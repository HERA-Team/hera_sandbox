import aipy as a, numpy as n

class AntennaArray(a.fit.AntennaArray):
    def get_params(self, ant_prms={'*':'*'}):
        prms = a.fit.AntennaArray.get_params(self, ant_prms)
        for k in ant_prms:
            top_pos = n.dot(self._eq2zen, self[int(k)].pos)
            if ant_prms[k] == '*':
                prms[k].update({'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]})
            else:
                for val in ant_prms[k]:
                    if   val == 'top_x': prms[k]['top_x'] = top_pos[0]
                    elif val == 'top_y': prms[k]['top_y'] = top_pos[1]
                    elif val == 'top_z': prms[k]['top_z'] = top_pos[2]
        return prms
    def set_params(self, prms):
        changed = a.fit.AntennaArray.set_params(self, prms)
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
            if ant_changed: ant.pos = n.dot(n.linalg.inv(self._eq2zen), top_pos)
            changed |= ant_changed
        return changed

L = 1000 / a.const.len_ns 
dL = 1000 / a.const.len_ns #10m packed grid spacing
antpos = []
cen_y, cen_z = 0, 0
for cen_x in n.arange(-16, 16, 1)*L:
    dx = 0
    for dy in n.arange(32)*dL:
        antpos.append((cen_x+dx, cen_y+dy, cen_z))

prms = {
    'loc': ('38:25:59.24',  '-79:51:02.1'), # Green Bank, WV
    'antpos': antpos,
    'beam': a.fit.Beam2DGaussian,
}

def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    for i in range(nants):
        #sigma = 0.085 for 5lambda dish
        beam = prms['beam'](freqs, xwidth=.085, ywidth=.085)
        antennas.append(a.fit.Antenna(0, 0, 0, beam))
    aa = AntennaArray(prms['loc'], antennas)
    p = {}
    for i in range(nants):
        top_pos = prms['antpos'][i]
        p[str(i)] = {'top_x':top_pos[0], 'top_y':top_pos[1], 'top_z':top_pos[2]}
    aa.set_params(p)
    return aa

def get_catalog(*args, **kwargs): return a.src.get_catalog(*args, **kwargs)
