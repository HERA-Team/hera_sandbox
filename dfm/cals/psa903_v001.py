"""
This is a cal file for data from psa32 maximum redundancy. Time-dependent gains are incorporated.
"""
import aipy as a
import numpy as n
import bm_prms as bm
from glob import glob
import ephem

prefix = '/data3/paper/eor_reduce'
npzfiles = [file.split('/')[-1] for file in glob(prefix+'/*_cal1.npz')]
npzfiles.sort()

class AntennaArray(a.pol.AntennaArray):
    def __init__(self, *args, **kwargs):
        a.pol.AntennaArray.__init__(self, *args, **kwargs)
        self.ant_layout = kwargs.pop('ant_layout')
        self.gain = kwargs.pop('gain')
        self.npzfiles = npzfiles
        a.pol.AntennaArray.update(self)
    def update_gains(self,npz):
        for i,gain in zip(npz['inpnums'].flatten(), npz['C_gain'].flatten()):
            print 'updating antenna %d with gain %0.2f'%(i,10**gain)
            self[i].set_params({'amp':prms['amps'].get(i,1.)*10**gain}) 
    def update_delays(self,npz):
        for i,tau in zip(npz['inpnums'].flatten(), npz['C'].flatten()):
            print 'updating antenna %d with delay %0.2f'%(i,tau)
            self[i].set_params({'dly':prms['delays'].get(i,0.)+tau}) 
    def update(self,npz):
        self.update_gains(npz)
        self.update_delays(npz)
        a.pol.AntennaArray.update(self)
    
    def set_jultime(self,t):
        #When this is called, we want to search for a cal0.npz file, eg, and apply the gainsolution therewithin
        a.pol.AntennaArray.set_jultime(self,t)
        for file in self.npzfiles: #which file corresponds to our jultime?
            filetime = float(file.split('.')[1]+'.'+file.split('.')[2])
            if t > filetime: pass
            else: break
        if t - filetime <= 11./(60.*24.): #make sure we're in the same file.
            file = prefix+'/'+file
            self.update(n.load(file))
        else: pass

    def get_params(self, ant_prms={'*':'*'}):
        try: prms = a.pol.AntennaArray.get_params(self, ant_prms)
        except(IndexError): return {}
        for k in ant_prms:
            if k == 'aa':
                if not prms.has_key('aa'): prms['aa'] = {}
                for val in ant_prms[k]:
                    if val == 'gain': prms['aa']['gain'] = self.gain
            else:
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
        try: self.gain, changed = prms['aa']['gain'], 1
        except(KeyError): pass
        if changed: self.update()
        return changed

prms = {
    'loc': ('-30:43:17.5', '21:25:41.9'), # KAT, SA (GPS)
    'antpos': {
        3  :[-1.8703164651063193, -700.38179398735906, -2.486720328795109] ,
        2  : [4.8446646530159176, -700.39598247763035, 9.0712755708770896] ,
        1  : [11.591949094991772, -700.51959903136083, 20.550883906201953] ,
        0  : [18.303883632714129, -700.54155018781171, 32.147239675343215] ,

        7  :[-0.79691390069416523, -299.77026013347034, -1.1324501024272999] ,
        6  : [5.8589125156210056, -299.78452412641195, 10.430549257939724] ,
        5  : [12.570956582720362, -299.96597679156088, 22.013562464192603] ,
        4  : [19.315903718703861, -300.06359266903411, 33.543205414232105] ,

        11 :[-1.3878432448048694, -500.38082812384613, -1.7545471394250993] ,
        10 : [5.3138669086039991, -500.49730834340397, 9.946881319840271] ,
        9  : [12.000410085638356, -500.52004427391256, 21.494870295016263] ,
        8  : [18.700416704947475, -500.51760577469884, 32.967807350450215] ,

        15 :[-0.30041897106453369, -99.592047095062725, -0.4986783217372045] ,
        14 : [6.4059753694621326, -99.95224098313237, 10.967587452252552] ,
        13 : [13.103671269187693, -100.16568752920908, 22.495562579989407] ,
        12 : [19.71347252253916, -100.15113820962031, 33.991849122807949] ,

        19 : [-1.4502397412457093, -600.18812435901839, -1.9179935464574072] ,
        18 : [5.1717744224385322, -600.38441009883991, 9.5482722259790709] ,
        17 : [11.92743084773735, -600.50952069915058, 21.087922100514781] ,
        16 : [18.660932830697366, -600.37151947354073, 32.734312482817408] ,

        23 :[-0.65004684427314674, -199.90294892149257, -0.70882370043213849] ,
        22 : [6.085058720866952, -199.91607722453426, 10.832493993299661] ,
        21 : [12.764558079120691, -200.09371555153467, 22.343790916649397] ,
        20 : [19.490548357060103, -200.03876633995446, 33.933475405899749] ,

        27 :[-1.1241140538220571, -400.00575563876032, -1.5594121443679183] ,
        26 : [5.627745048786398, -400.21358435343956, 10.011926318192666] ,
        25 : [12.334663628936832, -400.32135868916214, 21.53990144903608] ,
        24 : [19.066037795876635, -400.28464748951757, 33.111239910043388] ,

        31 :[0.0, 0.0, 0.0] ,
        30 :[6.6895727415641932, 0.035200138018290757, 11.576341924808688] ,
        29 : [13.380230127377081, 0.14164134349682395, 23.136005645230252] ,
        28 : [20.113783156528598, 0.12912235844973025, 34.702340643989459] ,
    },
    'ant_layout': n.array(
        [[0, 16, 8, 24, 4, 20, 12, 28],
         [1, 17, 9, 25, 5, 21, 13, 29],
         [2, 18, 10, 26, 6, 22, 14, 30],
         [3, 19, 11, 27, 7, 23, 15, 31]]),
    'delays': {},
    'amps': {},
    'offsets': { 8: 0.5 },
    'gain': .0036096,
    'bp_r': n.array([[-546778459030.53168, 664643788581.23596, -352000715429.32422, 106069000024.00294, -19886868672.0816, 2375187771.2150121, -176441928.4305163, 7452103.7565970663, -136981.43950786022]] * 64) * 1.0178**0.5, # from J2214-170 in Helmboldt, with adjustment for tempgain.py gain adjustment
    'bp_i': n.array([[0., 0., 0.]] * 32),
    'beam': a.fit.BeamAlm,
   'twist': n.array([0]*32),#n.array([.1746] * 32),
}

def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    for i in range(nants):
        for pi in ('x'):
            beam = bm.prms['beam'](freqs, nside=32, lmax=20, mmax=20, deg=7)
            try: beam.set_params(bm.prms['bm_prms'])
            except(AttributeError): pass
            pos = prms['antpos'][i]
            dly = prms['delays'].get(i, 0.)
            off = prms['offsets'].get(i, 0.)
            amp = prms['amps'].get(i, 1.)
            bp_r = prms['bp_r'][i]
            bp_i = prms['bp_i'][i]
            twist = prms['twist'][i]
            antennas.append(a.pol.Antenna(pos[0], pos[1], pos[2], beam, num=i, pol=pi, phsoff=[dly,off],
                    amp=amp, bp_r=bp_r, bp_i=bp_i, pointing=(0.,n.pi/2,twist),lat=prms['loc'][0]))
    aa = AntennaArray(prms['loc'], antennas, gain=prms['gain'], ant_layout=prms['ant_layout'])
    return aa

src_prms = {
'cen':{ 'jys':10**3.282102, 'index':  0.235166 , },
'cyg':{ 'jys':10**3.566410, 'index':  -0.266315 , },
'hyd':{ 'jys':10**2.448816, 'index':  -0.866462 , },
'pic':{ 'jys':10**2.714456, 'index':  -0.436361 , },
'vir':{ 'jys':10**2.200725, 'index':  0.202425 , },
'Sun': {'a1': 0.00644, 'index': 1.471, 'a2': 0.00586, 'jys': 55701.96, 'th': -0.000512},
'for': {'a1': 0.00851, 'a2': 0.00413, 'jys': 907.09, 'th': 0.230},
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

