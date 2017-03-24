import aipy as a, numpy as n, os

class AntennaArray(a.pol.AntennaArray):
    def __init__(self, *args, **kwargs):
        a.pol.AntennaArray.__init__(self, *args, **kwargs)
        self.array_params = {}
    def get_ant_params(self, ant_prms={'*':'*'}):
        prms = a.fit.AntennaArray.get_params(self, ant_prms)
        for k in ant_prms:
            top_pos = n.dot(self._eq2zen, self[int(k)].pos)
            if ant_prms[k] == '*':
                prms[k].update({'top_x':top_pos[0],'top_y':top_pos[1],'top_z':top_pos[2]})
            else:
                for val in ant_prms[k]:
                    if   val == 'top_x': prms[k]['top_x'] = top_pos[0]
                    elif val == 'top_y': prms[k]['top_y'] = top_pos[1]
                    elif val == 'top_z': prms[k]['top_z'] = top_pos[2]
        return prms
    def set_ant_params(self, prms):
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
    def get_arr_params(self):
        return self.array_params
    def set_arr_params(self, prms):
        for param in prms:
            self.array_params[param] = prms[param]
            if param == 'dish_size_in_lambda':
                FWHM = 2.35*(0.45/prms[param]) #radians
                self.array_params['obs_duration'] = 60.*FWHM / (15.*a.const.deg)# minutes it takes the sky to drift through beam FWHM
            if param == 'antpos':
                bl_lens = n.sum(n.array(prms[param])**2,axis=1)**.5
                self.array_params['uv_max'] = n.ceil(n.max(bl_lens)) #longest baseline
        return self.array_params

#===========================ARRAY SPECIFIC PARAMETERS==========================

#Set antenna positions here; for regular arrays like Hera we can use an algorithm; otherwise antpos should just be a list of [x,y,z] coords in light-nanoseconds
""" #HERA
nside = 7. #hex number
L = 1400 / a.const.len_ns 
dL = 1212 / a.const.len_ns #close packed hex
antpos = []
cen_y, cen_z = 0, 0
for row in n.arange(nside):
    for cen_x in n.arange((2*nside-1)-row):
        dx = row/2
        antpos.append(((cen_x + dx)*L, row*dL, cen_z))
        if row != 0:
            antpos.append(((cen_x + dx)*L, -row*dL, cen_z))
""" 
#PAPER128
antlayout = n.array(
[[ 64,10,49,3,41,25,19,48,29,24,28,55,34,27,51,57 ],
[ 65,9,66,58,67,1,47,4,68,17,69,13,70,56,71,59 ],
[ 72,22,73,61,74,35,75,18,76,5,77,32,78,30,79,23 ],
[ 80,20,81,63,82,42,83,37,84,40,85,14,86,54,87,50 ],
[ 88,43,89,2,90,33,91,6,92,52,93,7,94,12,95,38 ],
[ 96,53,97,21,98,15,99,16,100,62,101,44,102,0,103,26 ],
[ 104,31,105,45,106,8,107,11,108,36,109,60,110,39,111,46 ]]).flatten()
#badants = [2,10,15,22,31,33,42,43,47,58,64,72,91,97,105,107,34,84,100,56,7] #S1E1 omni_v3_xtalk
badants = [34,84,100,56,7] #S1E2 omni_v3_xtalk
badinds = []
for ba in badants:
    badinds.append(n.where(antlayout == ba)[0][0]) #where bad antennas are located in the array
badinds_x = n.array(badinds)/16
badinds_y = n.array(badinds) - 16*badinds_x
badcoords = []
for b, ba in enumerate(badinds_x):
    badcoords.append((badinds_x[b],badinds_y[b])) #x,y indices of bad antennas in the array
antpos = []
factor = 3.335640962 #m to light-ns
ypos = 0.0
for x in range(7):
    xpos = 0.0
    for y in range(16):
        if (x,y) in badcoords: pass #skip over baselines w/bad antennas
        else: antpos.append((xpos,ypos,0.0))
        xpos += 15.0*factor  
    ypos += 4.0*factor
#antpos = [(0.0, 0.0, 0.0),(30.0*factor, 0.0, 0.0)]
prms = {
    'name': os.path.basename(__file__)[:-3], #remove .py from filename
    'loc': ('-30:43:17.5', '21:25:41.9'), # KAT, SA (GPS)
    'antpos': antpos,
    'beam': a.fit.Beam2DGaussian,
    'dish_size_in_lambda': 0.82, #1., #in units of wavelengths at 150 MHz = 2 meters; this will also define the observation duration (ex: HERA dishes are 14-m = 7 wavelengths)
    'Trx': 200e3 #5e5 #receiver temp in mK
}

#=======================END ARRAY SPECIFIC PARAMETERS==========================

def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    for i in range(nants):
        beam = prms['beam'](freqs, xwidth=(0.45/prms['dish_size_in_lambda']), ywidth=(0.45/prms['dish_size_in_lambda'])) #as it stands, the size of the beam as defined here is not actually used anywhere in this package, but is a necessary parameter for the aipy Beam2DGaussian object
        antennas.append(a.fit.Antenna(0, 0, 0, beam))
    aa = AntennaArray(prms['loc'], antennas)
    p = {}
    for i in range(nants):
        top_pos = prms['antpos'][i]
        p[str(i)] = {'top_x':top_pos[0],'top_y':top_pos[1],'top_z':top_pos[2]}
    aa.set_ant_params(p)
    aa.set_arr_params(prms) 
    return aa

def get_catalog(*args, **kwargs): return a.src.get_catalog(*args, **kwargs)
