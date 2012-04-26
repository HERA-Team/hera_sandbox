''' This is a calibration file for data collected at PAPER in Karoo, South Africa
on JD 2455819 '''

import aipy as a, numpy as n,glob,ephem
import bm_prms as bm
import generic_catalog
import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('psa743_v004_gc')

class Antenna(a.pol.Antenna):
    def __init__(self,*args,**kwargs):
        a.pol.Antenna.__init__(self,*args,**kwargs)
        self.dpos = kwargs['dpos']
        self._pos = self.pos.copy()
        self.update()
    def update(self):
        a.pol.Antenna.update(self)
        self.pos = self._pos + self.dpos
        self._bm_gain = a.pol.Antenna.bm_response(self,(0,0,1),pol=self.pol)[0,0]
    def bm_response(self,*args,**kwargs):
        return a.pol.Antenna.bm_response(self,*args,**kwargs)/self._bm_gain
    def get_params(self, prm_list=['*']):
        """Return all fitable parameters in a dictionary."""
        x,y,z = self._pos
        aprms = {'x':x, 'y':y, 'z':z, 'dly':self._phsoff[-2],
            'off':self._phsoff[-1], 'phsoff':self._phsoff}
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

prms = {
    'loc': ('-30:43:17.5', '21:25:41.9'), # KAT, SA (GPS)
    'antpos':{
        0:[147.659407413, 336.269469733, 264.566180759],
        1:[-120.566931266, -270.142735412, -201.208899961],
        2:[175.483874,-282.046474,309.593714],
        3:[-24.5939776529, -369.022493234, -35.9049669793],
        #--------------------------------------------------------
        4:[-135.977107,-65.373043,-223.715356],
        5:[-184.222167454,  60.9256169154, -307.675312464],
        6:[84.568610,-397.906007,151.703088],
        7:[60.9037241018, 422.222408268, 116.124879563],
        #--------------------------------------------------------
        8:[148.405177,-231.475974,263.305593],
        9:[-121.15655,-329.37685,-197.06224],
        10:[-28.800063,-420.849441,-43.564604],
        11:[-180.112865,-190.297251,-301.062917],
        #--------------------------------------------------------
        12:[161.032208592, 207.530151484, 286.703007713],
        13:[-79.830818,266.416356,-122.828467],
        14:[90.491568,406.666552,171.303074],
        15:[136.833937217,-349.10409, 256.16691],
        #========================================================
	    16:[75.008275,-366.663944,135.807286],
        17:[-170.082246,113.392564,-280.090332],
        18:[-173.308984588, -52.5844630491, -289.495946419],
        19:[35.6156894023, -76.4247822222, 68.8003235664],
        #-------------------------------------------------------
        20:[ 223.405495506, -111.371927382, 391.352958368],
        21:[ 211.984088554, -181.820834933, 372.672243377],
        22:[-52.9274701935, -409.284993158, -84.1268196632],
        23:[-75.327786,379.129646,-113.829018],
        #--------------------------------------------------------
        24:[-90.374808,3.548977,-144.207995],
        25:[-23.653561,-153.921245,-31.289596],
        26:[208.418197,189.287085,370.725255],
        27:[-22.2282015089, 311.926612877, -26.8228657991],
        #--------------------------------------------------------
        28:[-18.1453146192, 166.083642242, -21.2052534495],
        29:[89.6597220746, -22.1294190136, 162.698139384],
        30:[-139.053365,312.917932,-223.870462],
        31:[229.945829,48.161862,406.414507],
        #--------------------------------------------------------
        32:[-112.893563,109.228967,-182.880941],
        33:[121.355347,-319.429590,209.575748],
        34:[-1.186004,298.790781,-1.572735],
        35:[-150.754218,-224.460782,-258.594058],
        #--------------------------------------------------------
        36:[-148.166345,285.390149,-254.152706],
        37:[73.704070,-378.161280,127.753480],
        38:[183.238623,145.046381,314.997386],
        39:[201.110057,270.608943,345.388038],
        #--------------------------------------------------------
        40:[-187.753175,101.634584,-322.330703],
        41:[32.859445,-311.361270,57.492402],
        42:[111.791791,-360.752264,193.124569],
        43:[185.296482,12.473870,318.948404],
        #--------------------------------------------------------
        44:[66.840886,269.989165,115.139909],
        45:[208.327549,-181.024029,358.713760],
        46:[222.401981,114.559981,382.329808],
        47:[82.998742,-157.005822,143.375763],
        #-------------------------------------------------------
        48:[-123.364050,7.568406,-211.391982],
        49:[42.324815,-394.596554,73.800150],
        50:[155.428104,103.981800,267.545140],
        51:[4.002712,454.858259,7.086482],
        #-------------------------------------------------------
        52:[40.840441,382.998141,70.689703],
        53:[228.948582,78.038958,393.662509],
        54:[208.232148,171.396294,357.761446],
        55:[22.162702,221.120016,38.449461],
        #--------------------------------------------------------
        56:[-85.962903,360.456826,-147.018238],
        57:[-22.182170,447.517664,-37.585541],
        58:[-40.132905,-349.207661,-68.174661],
        59:[-38.864384,362.866457,-66.270033],
        #--------------------------------------------------------
        60:[134.062901,401.074665,230.468279],
        61:[-81.496611,-277.174777,-139.301327],
        62:[-161.608043,226.512058,-277.243397],
        63:[170.275635,-299.764724,293.554481],
    }, 
    'dpos':{
        },
    'delays': {
     33 : {'x':  4.13077454567, 'y':  4.13764266968},
     34 : {'x':0.0027533769607, 'y':  1.19679474831},
     35 : {'x':  -1.5557454139, 'y': -1.31143268347},
     36 : {'x': -0.30007629469, 'y': -11.6215662467},
     37 : {'x':  -1.5211947009, 'y': -4.44083583355},
     38 : {'x':  17.8601805056, 'y': -11.8484287703},
     39 : {'x': -0.66980201751, 'y':  2.17936239168}, 
     40 : {'x':  6.49049766175, 'y': -13.4600170244},
     41 : {'x':  6.69653853768, 'y': -2.09064245224},
     42 : {'x':  5.79218932521, 'y': -4.63502132893},
     43 : {'x': -3.79584121704, 'y':  -9.7678035412},
     44 : {'x':-0.476948642731, 'y': -4.75042368418},
     45 : {'x':  4.91747128064, 'y':-0.277057933585},
     46 : {'x': -0.10716642797, 'y':-0.870703578671},
     47 : {'x':  1.75211724526, 'y':  -9.1290052938},
     48 : {'x':  4.85273297876, 'y':  7.26587415506},
     49 : {'x':  4.12397879362, 'y':   4.0667406589},
     50 : {'x':  1.54909890747, 'y': -12.6659270277},
     51 : {'x':  4.34739637189, 'y':  2.06444499756},
     52 : {'x':-0.458397638984, 'y': -6.60696714581},
     53 : {'x':  18.6865768433, 'y': -10.0768639594},
     54 : {'x':  4.54112693295, 'y':  4.91485251223},
     55 : {'x':  1.03537361526, 'y':  1.46194950109},
     56 : {'x':  3.46370260434, 'y':   14.108922866},
     57 : {'x':  5.31316812504, 'y':  5.39908059973},
     58 : {'x':  2.19500834905, 'y': -7.99523816478},
     59 : {'x':  3.04371864229, 'y':  17.1613820952},
     60 : {'x':  4.77957011759, 'y': 0.117690386636},
     61 : {'x': -1.97948118889, 'y':  11.3830283351},
     62 : {'x': -2.78464973342, 'y': -15.3278012643},
     63 : {'x': 0., 'y': 0.}, #REF
    },
    'amps': {
    33 : {'x': 9.72516325217, 'y': 9.73097880943}, 
    34 : {'x': 9.21297362356, 'y': 9.26477957328}, 
    35 : {'x': 7.51702254594, 'y': 7.55784919956}, 
    36 : {'x': 8.68938134663, 'y': 8.71941885817}, 
    37 : {'x': 9.21283616263, 'y': 9.22503959952},
    38 : {'x': 8.92463809164, 'y': 8.90650370966}, 
    39 : {'x': 9.23197678741, 'y': 9.19515292546}, 
    40 : {'x': 7.67449130184, 'y': 7.69111819703}, 
    41 : {'x': 8.96865066716, 'y': 8.98578994032},
    42 : {'x': 9.64412535757, 'y': 9.65422123645}, 
    43 : {'x': 9.08093315475, 'y': 9.05341998541}, 
    44 : {'x': 10.2476146848, 'y': 10.2581869408}, 
    45 : {'x': 9.67664888385, 'y': 9.67421618915}, 
    46 : {'x': 8.93627427798, 'y': 8.95537481592}, 
    47 : {'x': 9.01951289077, 'y': 9.01005555578}, 
    48 : {'x': 9.24325045912, 'y': 9.26984677324},
    49 : {'x': 8.51880738617, 'y': 8.52651156673}, 
    50 : {'x': 8.82618366191, 'y': 8.79400049108},  
    51 : {'x': 9.43840154501, 'y': 9.41583869828},  
    52 : {'x':  9.3264037881, 'y': 9.27877384586}, 
    53 : {'x': 8.22367511605, 'y': 8.19652450092}, 
    54 : {'x': 8.52174281478, 'y': 8.50642350972},  
    55 : {'x': 10.7374271282, 'y': 10.7591473633}, 
    56 : {'x': 9.80584440267, 'y': 9.82587855102},
    57 : {'x': 5.66536042137, 'y': 5.61411067325}, 
    58 : {'x': 7.52875330119, 'y': 7.51075861386}, 
    59 : {'x': 9.74585275578, 'y': 9.77128267254}, 
    60 : {'x': 9.43308836029, 'y': 9.44136460781}, 
    61 : {'x': 7.53215954956, 'y': 7.50974073631}, 
    62 : {'x': 7.50142373412, 'y':  7.5120721268}, 
    63 : {'x':  9.6841489388, 'y': 9.67969912107},
    },
    'off':{
        },
    'bp_r':  n.array([-167333390752.98276, 198581623581.65594, -102487141227.4993, 30027423590.548084, -5459067124.669095, 630132740.98792362, -45056600.848056234, 1822654.0034047314, -31892.9279846797]) * 1.0178**0.5,
    '32_64_pos_offset': n.array([-13.44311,-21.21483,-2.31634])
}

def get_aa(freqs):
    '''Return the AntennaArray to be used for simulation.'''
    location = prms['loc']
    antennas = []
    nants = len(prms['antpos'])
    for pi in ('x','y'):
        for i in prms['antpos'].keys():
            beam = bm.prms['beam'](freqs,nside=32,lmax=20,mmax=20,deg=7)
            try: beam.set_params(bm.prms['bm_prms'])
            except(AttributeError): pass
            pos = prms['antpos'][i]
            if i > 32: pos -= prms['32_64_pos_offset'] 
            try: dly = prms['delays'][i][pi]
            except(KeyError): dly = 0.
            try: off = prms['off'][i][pi]
            except(KeyError): off = 0.
            bp_r = prms['bp_r']
            try: amp = prms['amps'][i][pi]
            except(KeyError): amp = 10.
            try: dpos = prms['dpos'][i]
            except(KeyError): dpos = n.array([0.,0.,0.,])
            antennas.append(
                Antenna(pos[0],pos[1],pos[2], beam, dpos=dpos, num=i, pol=pi, phsoff=[dly,off], amp=amp, bp_r = bp_r, lat=prms['loc'][0] )
                )
    aa = a.pol.AntennaArray(prms['loc'], antennas)
    return aa

src_prms = {
    'pic':{'jys':600.,'ra':'5:19:49.70','dec':'-45:46.45.0','mfreq':0.150,'index':-1.},
    'forA':{'jys':210.,'ra':'3:22:37','dec':'-37:07:36','mfreq':0.150,'index':-1.},
    'forB':{'jys':96.,'ra':'3:25:01','dec':'-37:14:06','mfreq':0.150,'index':-1.},
    'J0445-23':{'jys':174.,'ra':'4:45:43','dec':'-28:05:56','mfreq':0.150,'index':-1.},
    'J0522-36':{'jys':96.,'ra':'5:23:43','dec':'-36:22:48','mfreq':0.150,'index':-1.},
    'J0625-53':{'jys':60.,'ra':'6:25:11','dec':'-53:39:31','mfreq':0.150,'index':-1.},
}

def get_catalog(srcs=None, cutoff=None, catalogs=['helm','misc']):
    '''Return a catalog containing the listed sources.'''
    log.info("psa743_v002_gc")
    specials = ['pic','forA','forB','J0445-23','J0522-36','J0625-53']
    srclist =[]
    for c in catalogs:     
        log.info("looking for %s in a local file"%(c,))
        this_srcs = generic_catalog.get_srcs(srcs=srcs,
              cutoff=cutoff,catalogs=[c])
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
                if src in specials:
                    cat[src] = a.fit.RadioFixedBody(**src_prms[src])
    return cat

if __name__=='__main__':
    import sys, numpy as n
    if len(sys.argv)>1:
        print "loading catalog: ",sys.argv[1]
        logging.basicConfig(level=logging.DEBUG)
        cat = get_catalog(catalogs=[sys.argv[1]])
        names = [cat[src].src_name for src in cat]
        print "loaded",len(names)," sources"
        flx = [cat[src]._jys for src in cat]
        print names
        print "brightest source in catalog"
        print names[flx.index(n.max(flx))],n.max(flx)
        log.info("loaded %d items from %s"%(len(cat),sys.argv[1]))
        try: assert([cat[src].e_S_nu for src in cat])
        except(AttributeError): print "this catalog does not have flux errors"
