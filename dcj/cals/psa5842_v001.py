"""
This cal file applies to 2455842 to the end of that season

The calibration parameters,aa,positions,and antenna code came from psa746_v008 from capo/arp/calfiles

The beam is now imported from capo.EoR2012
"""


import aipy as a, numpy as n,ephem
from capo import EoR2012 as EOR
import generic_catalog
import logging
loglevel = logging.CRITICAL
logging.basicConfig(level=loglevel)
log = logging.getLogger('psa5842_v001')
log.setLevel(loglevel)
class Antenna(a.fit.Antenna):
    def __init__(self, *args, **kwargs):
        self.dphsoff = n.array([0., 0.])
        a.fit.Antenna.__init__(self,*args,**kwargs)
        self._pos = self.pos.copy()
        self.dpos = n.array([0., 0., 0.])
        self.update()
        self.lat = ephem.degrees(kwargs['lat'])
        self._eq2zen = a.coord.eq2top_m(0., self.lat)
        self.top_pos = n.dot(self._eq2zen,self.pos)
    def _update_phsoff(self):
        self.phsoff = n.polyval(self._phsoff + self.dphsoff, self.beam.afreqs)
    def update(self):
        a.fit.Antenna.update(self)
        #self._bm_gain = a.fit.Antenna.bm_response(self, (0,0,1), pol='x')[0,0]
        #self._bm_gain = n.average(a.fit.Antenna.bm_response(self, (0,0,1), pol='x')[:,0])
        self._bm_gain = 1
        self.pos = self._pos + self.dpos
    def bm_response(self, *args, **kwargs):
        #print '(', self._bm_gain, ')',
        return a.fit.Antenna.bm_response(self,*args,**kwargs) / self._bm_gain
    def get_params(self, prm_list=['*']):
        """Return all fitable parameters in a dictionary."""
        x,y,z = self._pos
        aprms = {'x':x, 'y':y, 'z':z, 'dly':self._phsoff[-2],
            'off':self._phsoff[-1], 'phsoff':self._phsoff}
        aprms['ddly'] = self.dphsoff[-2]
        aprms['doff'] = self.dphsoff[-1]
        aprms['dx'] = self.dpos[0]
        aprms['dy'] = self.dpos[1]
        aprms['dz'] = self.dpos[2]
        aprms['bp_r'] = list(self.bp_r)
        aprms['bp_i'] = list(self.bp_i)
        aprms['amp'] = self.amp
        aprms.update(self.beam.get_params(prm_list))
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, prms):
        """Set all parameters from a dictionary."""
        changed = False
        self.beam.set_params(prms)
        try: self._pos[0], changed = prms['x'], True
        except(KeyError): pass
        try: self._pos[1], changed = prms['y'], True
        except(KeyError): pass
        try: self._pos[2], changed = prms['z'], True
        except(KeyError): pass
        try: self.dpos[0], changed = prms['dx'], True
        except(KeyError): pass
        try: self.dpos[1], changed = prms['dy'], True
        except(KeyError): pass
        try: self.dpos[2], changed = prms['dz'], True
        except(KeyError): pass
        try: self._phsoff[-2], changed = prms['dly'], True
        except(KeyError): pass
        try: self._phsoff[-1], changed = prms['off'], True
        except(KeyError): pass
        try: self.dphsoff[-2], changed = prms['ddly'], True
        except(KeyError): pass
        try: self.dphsoff[-1], changed = prms['doff'], True
        except(KeyError): pass
        try: self._phsoff, changed = prms['phsoff'], True
        except(KeyError): pass
        try: self.bp_r, changed = prms['bp_r'], True
        except(KeyError): pass
        try: self.bp_i, changed = prms['bp_i'], True
        except(KeyError): pass
        try: self.amp, changed = prms['amp'], True
        except(KeyError): pass
        if changed: self.update()
        return changed

class AntennaArray(a.fit.AntennaArray):
    def set_jultime(self, t=None):
        a.fit.AntennaArray.set_jultime(self, t=t)
        if t == None: return
        # Implement time-dependent cal
#        if t > 2455747.0:  # corr restart preceding psa747 run (relative to psa746)
#            for i in range(48,64):
#                self[i].set_params({'ddly':5.056})
#        else: # default to psa746 phase calibration
#            for i in range(48,64): self[i].set_params({'ddly':0})
    def get_params(self, ant_prms={'*':'*'}):
        try: prms = a.fit.AntennaArray.get_params(self, ant_prms)
        except(IndexError): return {}
        for k in ant_prms:
            try: top_pos = n.dot(self._eq2zen, self[int(k)].pos)
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
        
prms = {
    #'loc': ('-30:45:40', '21:24:24.5'), # KAT, SA (Google)
    'loc': ('-30:43:17.5', '21:25:41.9'), # KAT, SA (GPS)
    'antpos':{
            #32:[-112.893563,109.228967,-182.880941], # original
            #47:[82.998742,-157.005822,143.375763], # original
            #53:[228.948582,78.038958,393.662509], # original
             0:n.array([ 147.93527, 337.12580, 264.51808]) + n.array([  0.0    ,  0.0    , 0.0    ]),
             1:n.array([ -122.020309402, -268.023014624, -199.815623761,]),
             2:n.array([ 178.182697611, -281.228145672, 314.256932747,]),
             3:n.array([ -27.26782,-368.33731, -35.49412]) + n.array([ -0.29182,  0.29581,-0.06470]),
             4:n.array([ -136.271417956, -63.8538192223, -223.81696088,]),
             5:n.array([ -185.394536198, 61.2146784634, -307.596987596,]),
             6:n.array([  81.62552,-398.52015, 151.70308]) + n.array([ -0.90456,  0.31155,-0.03229]),
             7:n.array([  61.18198, 423.06854, 115.99715]) + n.array([  0.04571,  0.00224, 0.06347]),
             8:n.array([ 151.744912834, -230.476259414, 266.100117346]),
             9:n.array([-121.15655,-329.37685,-197.06224]) + n.array([  0.07513,  0.00316,-0.14179]),
            10:n.array([ -31.90961,-420.62509, -43.56460]) + n.array([ -0.54834,  0.39122,-0.24124]),
            11:n.array([-181.51436,-188.96090,-301.06291]) + n.array([ -0.01202,  0.07766,-0.02666]),
            12:n.array([ 160.72973, 208.27653, 286.82157]) + n.array([ -0.00755,  0.06142,-0.02066]),
            13:n.array([ -77.8707662359, 265.443210445, -122.866728234,]),
            14:n.array([  93.49461, 405.98665, 171.30307]) + n.array([  0.01422, -0.06648, 0.09750]),
            15:n.array([ 137.15781,-345.85204, 247.10124]) + n.array([ -0.38124,  0.29203,-0.17590]),
            16:n.array([ 72.4889166113, -363.847833856, 135.520493153]),
            17:n.array([ -169.390511149, 113.82049335, -280.122814575]),
            18:n.array([ -174.729494724, -51.6560694198, -289.267372556]),
            19:n.array([ 33.674525626, -74.865368323, 69.2972343811,]),
            20:n.array([ 221.456027623, -108.390005006, 391.362891776]),
            21:n.array([ 210.58399,-180.86686, 373.06684]) + n.array([ -1.00517,  0.11158,-0.28632]),
            22:n.array([ -55.3656528503, -405.662993034, -84.3890558675]),
            23:n.array([ -72.5558848039, 377.393814897, -113.678876716]),
            24:n.array([ -90.34611,   4.21680,-144.20799]) + n.array([  0.01805,  0.10158,-0.09630]),
            25:n.array([ -24.79049,-153.74222, -31.28959]) + n.array([ -0.35529,  0.45195,-0.14708]),
            26:n.array([ 210.10061544, 187.975673579, 370.683657743]),
            27:n.array([ -21.84807, 312.60274, -26.72816]) + n.array([  0.04485,  0.29455,-0.25724]),
            28:n.array([ -18.8534131847, 166.106071174, -21.0403737068]),
            29:n.array([  88.43837, -21.20718, 162.84732]) + n.array([  0.06337,  0.22358,-0.09150]),
            30:n.array([-136.73690, 313.93707,-223.87046]) + n.array([  0.00693,  0.14085,-0.14483]),
            31:n.array([ 230.19780727, 47.1563467525, 406.219047894,]),
            32:n.array([-121.89484, 147.01857,-178.53714]) + n.array([-12.34044,-21.29123,-1.62317]),
            33:n.array([ 121.35534,-319.42959, 209.57574]) + n.array([-13.46091,-21.08721,-2.04411]),
            34:n.array([  -1.18600, 298.79078,  -1.57273]) + n.array([-12.15281,-21.09783,-1.02046]),
            35:n.array([-150.75421,-224.46078,-258.59405]) + n.array([-13.33191,-21.31901,-2.09637]),
            36:n.array([-148.16634, 285.39014,-254.15270]) + n.array([-12.14134,-20.96592,-1.06129]),
            37:n.array([ 65.4745545234, -398.989793689, 135.937771903]),
            38:n.array([ 183.23862, 145.04638, 314.99738]) + n.array([-12.53417,-21.24913,-1.65605]),
            39:n.array([ 201.11005, 270.60894, 345.38803]) + n.array([-12.49116,-21.34547,-1.58327]),
            40:n.array([-187.75317, 101.63458,-322.33070]) + n.array([-13.027  ,-20.968  ,-1.872  ]),
            41:n.array([  32.85944,-311.36127,  57.49240]) + n.array([-13.45655,-21.32004,-2.25554]),
            42:n.array([ 111.79179,-360.75226, 193.12456]) + n.array([-13.38514,-21.23868,-1.85556]),
            43:n.array([ 185.29648,  12.47387, 318.94840]) + n.array([-12.90096,-21.24585,-1.81436]),
            44:n.array([  66.84088, 269.98916, 115.13990]) + n.array([-12.32095,-21.10367,-1.21838]),            
            45:n.array([ 208.32754,-181.02402, 358.71376]) + n.array([-13.77787,-21.21314,-2.35549]),
            46:n.array([ 222.40198, 114.55998, 382.32980]) + n.array([-12.61714,-21.02361,-1.42653]),
            47:n.array([  87.62899,-157.64380, 134.49244]) + n.array([-12.57070,-20.86601,-1.67331]),
            48:n.array([-123.36405,   7.56840,-211.39198]) + n.array([-12.63378,-21.12748,-1.07278]),
            49:n.array([  42.32481,-394.59655,  73.80015]) + n.array([-13.86259,-21.23736,-2.03213]),
            50:n.array([ 155.42810, 103.98180, 267.54514]) + n.array([-12.63677,-21.21007,-1.64005]),
            51:n.array([   4.00271, 454.85825,   7.08648]) + n.array([-11.63465,-20.89803,-1.78326]),
            52:n.array([28.8976025965, 358.993678483, 69.7233597137,]),
            53:n.array([ 247.12661,  75.95094, 404.94444]) + n.array([-13.14037,-20.72110,-1.99411]),
            54:n.array([ 195.692326034, 148.9067559, 358.382640028]),
            55:n.array([  22.16270, 221.12001,  38.44946]) + n.array([-13.027  ,-20.968  ,-1.872  ]),
            56:n.array([ -85.96290, 360.45682,-147.01823]) + n.array([-11.24916,-21.18504,-1.59575]),
            57:n.array([ -22.18217, 447.51766, -37.58554]) + n.array([-11.87513,-21.21941,-1.20133]),
            58:n.array([ -40.13290,-349.20766, -68.17466]) + n.array([-13.60638,-21.21031,-2.56979]),
            59:n.array([ -38.86438, 362.86645, -66.27003]) + n.array([-11.69585,-21.02930,-1.18640]),
            60:n.array([ 121.811369523, 377.231448971, 229.602800651]),
            61:n.array([ -94.7836220969, -297.64232068, -141.49713256]),
            62:n.array([-161.60804, 226.51205,-277.24339]) + n.array([-12.05381,-20.93013,-1.23832]),
            63:n.array([170.275635,-299.76472, 293.55448]) + n.array([-13.44311,-21.21483,-2.31634]),
    },
    'delays': {
         0:  0.00000,  1:  1.90251,  2: -5.04434,  3:  1.82363,
         4: -9.59716,  5:  6.38722,  6: -7.75454,  7:  3.19018,
         8: -8.43995,  9: -0.03911, 10:-11.88403, 11: -6.97286,
        12:  2.48693, 13:-14.78799, 14:-10.46206, 15: -3.52409,
        16: 16.87281, 17: -0.92220, 18: -0.51058, 19:  3.25037,
        20:  2.60934, 21:  1.78428, 22:  5.27815, 23:-14.37782,
        24: -1.90652, 25:-15.40243, 26:-11.23833, 27:  3.14987,
        28:  3.87429, 29:  7.40496, 30: -2.08892, 31:-13.72179,
        32: 22.01993, 33: -2.60862, 34: -7.10621, 35: -8.65217,
        36: -7.53194, 37: -5.13895, 38: 10.45674, 39: -7.99611,
        40:  0.00000, 41: -0.52698, 42: -0.80384, 43:-10.87221,
        44: -7.79094, 45: -2.03461, 46: -7.14113, 47:-14.80853,
        48: -1.98915, 49: -2.47736, 50: -5.714,#score=0.7
        51: -3.63866,
        52: -8.13744, 53:  6.05448, 54: -1.43166, 55:  0.00000,
        56: -4.51703, 57: -2.00401, 58: -4.97970, 59: -4.43863,
        60: -2.48776, 61: -8.93252, 62: -9.89454, 63: -7.03737,
    },
    'offsets': { },
    'amps': {
  0:0.00339455692388,  1:0.00370372788408,  2:0.00376527818632,  3:0.00387665,
  4:0.00396057349076,  5:0.00410252231494,  6:0.00348871724217,  7:0.00392404211306,
  8:0.00375883047688,  9:0.00334462819991, 10:0.00337093254841, 11:0.00341326452831,
 12:0.00442460216876, 13:0.00389868,       14:0.00354483945983, 15:0.00407328010273,
 16:0.00403397,       17:0.00439585577112, 18:0.00431250874167, 19:0.00404608945244,
 20:0.00413586961583, 21:0.00403113024808, 22:0.00359282557633, 23:0.00394237335571,
 24:0.00366244535577, 25:0.00380304422308, 26:0.00326923586611, 27:0.00287590169933,
 28:0.00403816252674, 29:0.00401518091533, 30:0.00409408588555, 31:0.00359266622740,

 32:0.00382517,       33:0.00403864141466, 34:0.00358246,       35:0.00342863956488,
 36:0.00448295132661, 37:0.00372165748921, 38:0.00442515606401, 39:0.00435455776743,
 40:0.003,            41:0.00369075183672, 42:0.00357438,       43:0.00420578225324,
 44:0.00382302025597, 45:0.00367334705810, 46:0.00380154752356, 47:0.00333605,
 48:0.00375050553983, 49:0.00362668689090, 50:0.00349, #score=0.7
 51:0.00343844386287,
 52:0.00331377056150, 53:0.00391876407669, 54:0.00427390096123, 55:0.003,
 56:0.00392878,       57:0.00379400167760, 58:0.00356638,       59:0.00387608,
 60:0.00382425542649, 61:0.00336631508857, 62:0.00388450863696, 63:0.00370594585437,
    },

        'bp_r': n.array([[-167333390752.98276, 198581623581.65594, -102487141227.4993, 30027423590.548084,
        -5459067124.669095, 630132740.98792362, -45056600.848056234, 1822654.0034047314, -31892.9279846797]] * 64) *
        1.0178**0.5,
   # 'bp_r': n.array([[-546778459030.53168, 664643788581.23596, -352000715429.32422, 106069000024.00294, -19886868672.0816, 2375187771.2150121, -176441928.4305163, 7452103.7565970663, -136981.43950786022]] * 64), # from J2214-170 in Helmboldt
    #'bp_r': n.array([[-546778459030.53168, 664643788581.23596, -352000715429.32422, 106069000024.00294, -19886868672.0816, 2375187771.2150121, -176441928.4305163, 7452103.7565970663, -136981.43950786022]] * 64) * 1.0178**0.5, # from J2214-170 in Helmboldt, with adjustment for tempgain.py gain adjustment
    'bp_i': n.array([[0., 0., 0.]] * 64),
    'beam': a.fit.BeamAlm,
   'twist': n.array([0]*64)
   }



def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    for i in range(nants):
        beam = EOR.prms['beam'](freqs, nside=32, lmax=20, mmax=20, deg=7)
        try: beam.set_params(EOR.prms['bm_prms'])
        except(AttributeError): pass
        pos = prms['antpos'][i]
        dly = prms['delays'].get(i, 0.)
        off = prms['offsets'].get(i, 0.)
        amp = prms['amps'].get(i, 4e-3)
        bp_r = prms['bp_r'][i]
        bp_i = prms['bp_i'][i]
        twist = prms['twist'][i]
        antennas.append(Antenna(pos[0],pos[1],pos[2], beam, phsoff=[dly,off],
                amp=amp, bp_r=bp_r, bp_i=bp_i, pointing=(0.,n.pi/2,twist),lat=prms['loc'][0]))
        if i >= 32:
            antennas[-1].set_params({
               #'dx': -13.0274536654,
               #'dy': -20.9677193273,
               #'dz': -1.87182968283,
               #'ddly': -1.36854425919,
            })
    #aa = a.fit.AntennaArray(prms['loc'], antennas)
    aa = AntennaArray(prms['loc'], antennas)
    return aa


src_prms = {
'cen':{ 'jys':10**3.282102, 'index':  0.235166 , },
'cyg':{ 'jys':10**3.566410, 'index':  -0.266315 , },
'hyd':{ 'jys':10**2.448816, 'index':  -0.866462 , },
#'pic':{ 'jys':10**2.714456, 'index':  -0.436361 , },
'pic':{'jys':450, 'index':-1.},
'for':{'jys':267,'index':-1.15},#score=0.70
#{'jys':447,'index':-1.15},#old value
'vir':{ 'jys':10**2.200725, 'index':  0.202425 , },
'Sun': {'a1': 0.00644, 'index': 1.471, 'a2': 0.00586, 'jys': 55701.96, 'th': -0.000512},
#'for': {'a1': 0.00851, 'a2': 0.00413, 'jys': 907.09, 'th': 0.230},
}

def get_catalog(srcs=None, cutoff=None, catalogs=['helm','misc']):
    '''Return a catalog containing the listed sources.'''
    custom_srcs = ['J1347-603','J1615-610', 'J1336-340', 'J1248-412', 'J1531-423', 'J1359-415']
    srclist =[]
    for c in catalogs:
        log.info("looking for %s in a local file"%(c,))
        this_srcs = generic_catalog.get_srcs(srcs=srcs,
              cutoff=cutoff,catalogs=[c],loglevel=loglevel)
        if len(this_srcs)==0:
            log.warning("no sources found with genericfile, trying built in catalog")
            tcat = a.src.get_catalog(srcs=srcs, 
                   cutoff=cutoff, catalogs=[c])
            srclist += [tcat[src] for src in tcat]
        else: srclist += this_srcs    
    cat = a.fit.SrcCatalog(srclist)
    #Add specials.  All fixed radio sources must be in catalog, for completeness
    if not srcs is None:
        for src in srcs:
            if src in src_prms.keys():
                if src in custom_srcs:
                    cat[src] = a.fit.RadioSpecial(src,**src_prms[src])
    return cat    
#    if srcs is None:
#        cat = a.src.get_catalog(srcs=srcs, cutoff=cutoff, catalogs=catalogs)
#    else:
#        cat = a.src.get_catalog(srcs=[s for s in srcs if not s in custom_srcs],
#            cutoff=cutoff, catalogs=catalogs)
#        for src in [s for s in srcs if s in custom_srcs]:
#            cat[src] = a.fit.RadioFixedBody(0., 0., janskies=0., mfreq=.15, name=src)
    cat.set_params(src_prms)
    return cat
        
