"""
"""
import aipy as a, numpy as n
import bm_prms as bm

class Antenna(a.fit.Antenna):
    def __init__(self,*args,**kwargs):
        a.fit.Antenna.__init__(self,*args,**kwargs)
        self.dB = kwargs.pop('dB')
    def get_params(self, pol, prm_list=['*']):
        x,y,z = self.pos
        if self.dp:
           pi = {'x':0,'y':1}[pol]
           aprms = {'x':x, 'y':y, 'z':z, 'dly':self._phsoff[pi][-2], 
               'off':self._phsoff[pi][-1], 'phsoff':self._phsoff[pi]}
           aprms['bp_r'] = list(self.bp_r[pi])
           aprms['bp_i'] = list(self.bp_i[pi])
           aprms['amp'] = self.amp[pi]
           aprms['dB'] = self.dB[pi]
        else:
           aprms = {'x':x, 'y':y, 'z':z, 'dly':self._phsoff[-2], 
               'off':self._phsoff[-1], 'phsoff':self._phsoff}
           aprms['bp_r'] = list(self.bp_r)
           aprms['bp_i'] = list(self.bp_i)
           aprms['amp'] = self.amp
           aprms['dB'] = self.dB
        aprms.update(self.beam.get_params(prm_list)) 
        prms = {}
        for p in prm_list:
            if p.startswith('*'): return aprms
            try: prms[p] = aprms[p]
            except(KeyError): pass
        return prms
    def set_params(self, pol, prms):
            """Set all parameters from a dictionary."""
            changed = False
            self.beam.set_params(prms)
            if self.dp:
                pi = {'x':0,'y':1}[pol]
                try: self.pos[0], changed = prms['x'], True
                except(KeyError): pass
                try: self.pos[1], changed = prms['y'], True
                except(KeyError): pass
                try: self.pos[2], changed = prms['z'], True
                except(KeyError): pass
                try: self._phsoff[pi][-2], changed = prms['dly'], True
                except(KeyError): pass
                try: self._phsoff[pi][-1], changed = prms['off'], True
                except(KeyError): pass
                try: self._phsoff[pi], changed = prms['phsoff'], True
                except(KeyError): pass
                try: self.bp_r[pi], changed = prms['bp_r'], True
                except(KeyError): pass
                try: self.bp_i[pi], changed = prms['bp_i'], True
                except(KeyError): pass
                try: self.amp[pi], changed = prms['amp'], True
                except(KeyError): pass
                try: self.dB[pi], changed = prms['dB'], True
                except(KeyError): pass
            else:
                try: self.pos[0], changed = prms['x'], True
                except(KeyError): pass
                try: self.pos[1], changed = prms['y'], True
                except(KeyError): pass
                try: self.pos[2], changed = prms['z'], True
                except(KeyError): pass
                try: self._phsoff[-2], changed = prms['dly'], True
                except(KeyError): pass
                try: self._phsoff[-1], changed = prms['off'], True
                except(KeyError): pass
                try: self._phsoff, changed = prms['phsoff'], True
                except(KeyError): pass
                try: self.bp_r, changed = prms['bp_r'], True
                except(KeyError): pass
                try: self.bp_i, changed = prms['bp_i'], True
                except(KeyError): pass
                try: self.amp, changed = prms['amp'], True
                except(KeyError): pass
                try: self.dB, changed = prms['dB'], True
                except(KeyError): pass 
            if changed: self.update()
            return changed

class AntennaArray(a.fit.AntennaArray):
    def __init__(self,*args,**kwargs):
        a.fit.AntennaArray.__init__(self,*args,**kwargs)
        self.gain = kwargs.pop('gain')
        self.update_gains()
    def update_gains(self):
        for i in range(len(self)):
            self[i].amp = self.gain * 10 ** (self[i].dB / 10.) 
    def update(self):
        self.update_gains()
        a.fit.AntennaArray.update(self)

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
    'delays' : {
        # scores using pic,for
        #Overall, xx=59%; yy=62%
        #Column 0
         0: [  0.000,  -0.146],
         1: [ -3.433, -13.144],
         2: [  8.549,   1.508],
         3: [ 13.848,   1.504],
        #Column 4
         4: [ 15.892,   1.840],
         5: [  6.678,  -6.811],
         6: [  7.758,   1.582],
         7: [ 14.620,   3.751],
        #Column 2
         8: [ -5.381, -16.459],
         9: [ -7.033, -19.521],
        10: [ 10.684,   0.886],
        11: [ -5.649,  -5.339],
        #Column 6
        12: [ -4.403,  -3.207],
        13: [ 11.350,  -3.667],
        14: [ -0.906,   9.449],
        15: [  6.252,  -4.582],
        #Column 1
        16: [  2.462,  -6.240],
        17: [  9.967,   0.411],
        18: [ -1.858, -12.046],
        19: [ -6.087,   7.305],
        #Column 5
        20: [ 12.263,   1.189],
        21: [ -3.101,  -2.670],
        22: [ -6.140, -16.030],
        23: [  3.763,  -7.149],
        #Column 3
        24: [ 11.104,   0.253],
        25: [  0.665,   3.167],
        26: [ -0.738, -14.155],
        27: [  7.986,   1.622],
        #Column 7
        28: [  0.448,  -2.170],
        29: [  1.037,   0.851],
        30: [ -1.086,  12.766],
        31: [ -5.915, -15.286],
    },
    'gain': n.array([ 0.0041775,0.0040569]),
    'dB': {
        #Column 0
         0: n.array([ -0.411, -0.269]),
         1: n.array([ -0.515, -0.007]),
         2: n.array([ -0.610, -0.068]),
         3: n.array([ -0.276,  0.054]),
        #Column 4
         4: n.array([ -0.029,  0.094]),
         5: n.array([ -0.029,  0.117]),
         6: n.array([  0.361,  0.124]),
         7: n.array([  0.044,  0.279]),
        #Column 2
         8: n.array([  0.095, -0.014]),
         9: n.array([  0.034,  0.242]),
        10: n.array([  0.245,  0.066]),
        11: n.array([ -0.232, -0.151]),
        #Column 6
        12: n.array([ -0.199, -0.025]),
        13: n.array([  0.447, -0.083]),
        14: n.array([ -0.018,  0.313]),
        15: n.array([ -0.166,  0.229]),
        #Column 1
        16: n.array([  0.034,  0.161]),
        17: n.array([  0.475, -2.918]),
        18: n.array([  0.156, -0.403]),
        19: n.array([  0.235, -0.127]),
        #Column 5
        20: n.array([ -2.282, -0.261]),
        21: n.array([  0.585,  0.787]), 
        22: n.array([  0.156,  0.258]),
        23: n.array([  0.136,  0.225]),
        #Column 3
        24: n.array([  0.176, -0.351]),
        25: n.array([ -0.029,  0.271]),
        26: n.array([  0.418,  0.036]),
        27: n.array([  0.917,  0.735]),
        #Column 7
        28: n.array([ -0.793, -0.304]),
        29: n.array([  0.136,  0.316]),
        30: n.array([ -0.199,  0.245]),
        31: n.array([  0.185, -0.715]),
    },
    'bp_r': n.array([-546778459030.53168, 664643788581.23596, -352000715429.32422, 106069000024.00294, -19886868672.0816, 2375187771.2150121, -176441928.4305163, 7452103.7565970663, -136981.43950786022]) * 1.0178**0.5, # from J2214-170 in Helmboldt, with adjustment for tempgain.py gain adjustment
    'bp_i': n.array([[0., 0., 0.]]),
    'beam': a.fit.BeamAlm,
}

def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    for i in range(nants):
        beam = bm.prms['beam'](freqs, nside=32, lmax=20, mmax=20, deg=7)
        try: beam.set_params(bm.prms['bm_prms'])
        except(AttributeError): pass
        pos = prms['antpos'][i]
        try: dly = prms['delays'][i]
        except(KeyError): dly = [0.,0.]
        bp_r = [prms['bp_r'],prms['bp_r']]
        try: amp = prms['dB'][i]
        except(KeyError): amp = n.array([0.003,0.003])
        phsoff = [[dly[0],0.],[dly[1],0.]]
        twist = 0.
        antennas.append(
            Antenna(pos[0], pos[1], pos[2], beam, dp=True, phsoff=phsoff, dB=amp, bp_r=bp_r)
            )
    aa = AntennaArray(prms['loc'],antennas,gain=prms['gain'])
    return aa

src_prms = {
    'cen':{ 'jys':10**3.282102, 'index':   0.235166 , },
    'cyg':{ 'jys':10**3.566410, 'index':  -0.266315 , },
    'hyd':{ 'jys':10**2.448816, 'index':  -0.866462 , },
    'pic':{ 'jys':10**2.714456, 'index':  -0.436361 , },
    'vir':{ 'jys':10**2.200725, 'index':   0.202425 , },
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

