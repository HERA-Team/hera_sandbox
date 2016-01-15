import numpy as np, omnical
import numpy.linalg as la
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=DeprecationWarning)
    import scipy.sparse as sps
    
POL_TYPES = 'xylrab'
#XXX this can't support restarts or changing # pols between runs
POLNUM = {} # factor to multiply ant index for internal ordering, 
NUMPOL = {}
#    'x':0, 
#    'y':1,
#    'l':2,
#    'r':3,
#    'a':4,
#    'b':5,
#}

def add_pol(p):
    global NUMPOL
    assert(p in POL_TYPES)
    POLNUM[p] = len(POLNUM)
    NUMPOL = dict(zip(POLNUM.values(), POLNUM.keys()))
    
class Antpol:
    def __init__(self, *args):
        try:
            ant,pol,nant = args
            if not POLNUM.has_key(pol): add_pol(pol)
            self.val, self.nant = POLNUM[pol] * nant + ant, nant
        except(ValueError): self.val, self.nant = args
    def antpol(self): return self.val % self.nant, NUMPOL[self.val / self.nant]
    def ant(self): return self.antpol()[0]
    def pol(self): return self.antpol()[1]
    def __int__(self): return self.val
    def __hash__(self): return self.ant()
    def __str__(self): return ''.join(map(str, self.antpol()))
    def __eq__(self, v): return self.ant() == v
    def __repr__(self): return str(self)
        
## XXX filter_reds w/ pol support should probably be in omnical
def filter_reds(reds, bls=None, ex_bls=None, ants=None, ex_ants=None, ubls=None, ex_ubls=None, crosspols=None, ex_crosspols=None):
    '''Filter redundancies to include/exclude the specified bls, antennas, and unique bl groups and polarizations.
    Assumes reds indices are Antpol objects.'''
    def pol(bl): return bl[0].pol() + bl[1].pol()
    if crosspols: reds = [r for r in reds if pol(r[0]) in crosspols]
    if ex_crosspols: reds = [r for r in reds if not pol(r[0]) in ex_crosspols]
    return omnical.arrayinfo.filter_reds(reds, bls=bls, ex_bls=ex_bls, ants=ants, ex_ants=ex_ants, ubls=ubls, ex_ubls=ex_ubls)

class RedundantInfo(omnical.info.RedundantInfo):
    def __init__(self, nant, filename=None):
        omnical.info.RedundantInfo.__init__(self, filename=filename)
        self.nant = nant
    def bl_order(self):
        '''Return (i,j) baseline tuples in the order that they should appear in data.  Antenna indicies
        are in real-world order (as opposed to the internal ordering used in subsetant).'''
        return [(Antpol(self.subsetant[i],self.nant),Antpol(self.subsetant[j],self.nant)) for (i,j) in self.bl2d]
    def order_data(self, dd):
        '''Create a data array ordered for use in _omnical.redcal.  'dd' is
        a dict whose keys are (i,j) antenna tuples; antennas i,j should be ordered to reflect
        the conjugation convention of the provided data.  'dd' values are 2D arrays
        of (time,freq) data.'''
        d = []
        for i,j in self.bl_order():
            bl = (i.ant(),j.ant())
            pol = i.pol() + j.pol()
            try: d.append(dd[bl][pol])
            except(KeyError): d.append(dd[bl[::-1]][pol[::-1]].conj())
        return np.array(d).transpose((1,2,0))

def compute_reds(nant, *args, **kwargs):
    reds = omnical.arrayinfo.compute_reds(*args, **kwargs)
    return [map(lambda bl: (Antpol(bl[0],nant),Antpol(bl[1],nant)), gp) for gp in reds]
    
def aa_to_info(aa, pols=['x'], **kwargs):
    '''Use aa.ant_layout to generate redundances based on ideal placement.
    The remaining arguments are passed to omnical.arrayinfo.filter_reds()'''
    layout = aa.ant_layout
    nant = len(aa)
    antpos = -np.ones((nant*len(pols),3)) # -1 to flag unused antennas
    xs,ys = np.indices(layout.shape)
    for ant,x,y in zip(layout.flatten(), xs.flatten(), ys.flatten()):
        for z,pol in enumerate(pols):
            z = 2**z # exponential ensures diff xpols aren't redundant w/ each other
            i = Antpol(ant,pol,len(aa)) # creates index in POLNUM/NUMPOL for pol
            antpos[i,0],antpos[i,1],antpos[i,2] = x,y,z
    reds = compute_reds(nant, antpos,tol=.1)
    # XXX haven't enforced xy = yx yet.  need to conjoin red groups for that
    ex_ants = [Antpol(i,nant).ant() for i in range(antpos.shape[0]) if antpos[i,0] < 0]
    kwargs['ex_ants'] = kwargs.get('ex_ants',[]) + ex_ants
    reds = filter_reds(reds, **kwargs)
    info = RedundantInfo(nant)
    info.init_from_reds(reds,antpos)
    return info


def redcal(data, info, xtalk=None, gains=None, vis=None,removedegen=False, uselogcal=True, maxiter=50, conv=1e-3, stepsize=.3, computeUBLFit=True, trust_period=1):
    # XXX add layer to support new gains format
    if gains:
        _gains = {}
        for pol in gains:
            for i in gains[pol]:
                ai = Antpol(i,pol,info.nant)
                _gains[int(ai)] = gains[pol][i].conj()
    else: _gains = gains
    if vis:
        _vis = {}
        for pol in vis:
            for i,j in vis[pol]:
                ai,aj = Antpol(i,pol[0],info.nant), Antpol(j,pol[1],info.nant)
                _vis[(int(ai),int(aj))] = vis[pol][(i,j)]
    else: _vis = vis
    meta, gains, vis = omnical.calib.redcal(data, info, xtalk=xtalk, gains=_gains, vis=_vis, removedegen=removedegen, uselogcal=uselogcal, maxiter=maxiter, conv=conv, stepsize=stepsize, computeUBLFit=computeUBLFit, trust_period=trust_period)    
    # rewrap to new format
    def mk_ap(a): return Antpol(a, info.nant)
    for i,j in meta['res'].keys():
        api,apj = mk_ap(i),mk_ap(j)
        pol = api.pol() + apj.pol()
        bl = (api.ant(), apj.ant())
        if not meta['res'].has_key(pol): meta['res'][pol] = {}
        meta['res'][pol][bl] = meta['res'].pop((i,j))
    #XXX make chisq a nested dict, with individual antpol keys?
    for k in [k for k in meta.keys() if k.startswith('chisq')]:
        try:
            ant = int(k.split('chisq')[1])
            meta['chisq'+str(mk_ap(ant))] = meta.pop(k)
        except(ValueError): pass
    for i in gains.keys():
        ap = mk_ap(i)
        if not gains.has_key(ap.pol()): gains[ap.pol()] = {}
        gains[ap.pol()][ap.ant()] = gains.pop(i).conj()
    for i,j in vis.keys():
        api,apj = mk_ap(i),mk_ap(j)
        pol = api.pol() + apj.pol()
        bl = (api.ant(), apj.ant())
        if not vis.has_key(pol): vis[pol] = {}
        vis[pol][bl] = vis.pop((i,j))
    return meta, gains, vis

def compute_xtalk(res, wgts):
    '''Estimate xtalk as time-average of omnical residuals.'''
    xtalk = {}
    for key in res:
        r,w = np.where(wgts[key] > 0, res[key], 0), wgts[key].sum(axis=0)
        w = np.where(w == 0, 1, w)
        xtalk[key] = (r.sum(axis=0) / w).astype(res[key].dtype) # avg over time
    return xtalk

def to_npz(filename, meta, gains, vismdl, xtalk, jds, lsts, freqs, conj=True):
    '''Write results from omnical.calib.redcal (meta,gains,vismdl) to npz file.
    Each of these is assumed to be a dict keyed by pol, and then by bl/ant/keyword'''
    d = {}
    for pol in meta:
        for k in meta[pol]:
            if k.startswith('chisq'):
                d[k + ' %s' % pol] = meta[pol][k]
            if k.startswith('history'):
                d['history' + ' %s' % pol] = meta[pol][k]
        # XXX chisq after xtalk removal
    for pol in gains:
        for ant in gains[pol]:
            if conj: d['%d%s' % (ant,pol)] = gains[pol][ant].conj() # conj to miriad 
            else: d['%d%s' % (ant,pol)] = gains[pol][ant]
    for pol in vismdl:
        for bl in vismdl[pol]:
            d['<%d,%d> %s' % (bl[0],bl[1],pol)] = vismdl[pol][bl]
    for pol in xtalk:
        for bl in xtalk[pol]: 
            d['(%d,%d) %s' % (bl[0],bl[1],pol)] = xtalk[pol][bl]
    d['jds'] =  jds
    d['lsts'] = lsts
    d['freqs'] = freqs
    np.savez(filename,**d)

def from_npz(filename, meta={}, gains={}, vismdl={}, xtalk={}, jds=[], lsts=[], freqs=[]):
    '''Reconstitute results from to_npz, returns meta, gains, vismdl, xtalk, each
    keyed first by polarization, and then by bl/ant/keyword.'''
    npz = np.load(filename)
    def parse_key(k):
        bl,pol = k.split()
        bl = tuple(map(int,bl[1:-1].split(',')))
        return pol,bl
    for k in [f for f in npz.files if f.startswith('(')]:
        pol,bl = parse_key(k)
        if not xtalk.has_key(pol): xtalk[pol] = {}
        xtalk[pol][bl] = npz[k]
    for k in [f for f in npz.files if f.startswith('<')]:
        pol,bl = parse_key(k)
        if not vismdl.has_key(pol): vismdl[pol] = {}
        vismdl[pol][bl] = npz[k]
    for k in [f for f in npz.files if f[0].isdigit()]:
        pol,ant = k[-1:],int(k[:-1])
        if not gains.has_key(pol): gains[pol] = {}
        gains[pol][ant] = npz[k]
    for k in [f for f in npz.files if f.startswith('c')]: #[0].isalpha()]:
        key,pol = k.split()
        if not meta.has_key(pol): meta[pol] = {}
        meta[pol][k.split()[0]] = npz[k]
    for k in [f for f in npz.files if f.startswith('h')]:
        key,pol = k.split()
        if not meta.has_key(pol): meta[pol] = {}
        meta[pol][k.split()[0]] = npz[k]
    for k in [f for f in npz.files if f.startswith('j')]: 
        jds = npz[k] 
    for k in [f for f in npz.files if f.startswith('l')]: 
        lsts = npz[k]
    for k in [f for f in npz.files if f.startswith('f')]:
        freqs = npz[k]
    return meta, gains, vismdl, xtalk, jds, lsts, freqs
