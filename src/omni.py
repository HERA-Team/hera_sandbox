import numpy as np, omnical

POLNUM = {
    'x':0, # factor to multiply ant index for internal ordering
    'y':1,
    'l':2,
    'r':3,
    'a':4,
    'b':5,
}

NUMPOL = dict(zip(POLNUM.values(), POLNUM.keys()))
    
class Antpol:
    def __init__(self, *args):
        try:
            ant,pol,nant = args
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
            i = Antpol(ant,pol,len(aa))
            antpos[int(i),0],antpos[int(i),1],antpos[int(i),2] = x,y,z
    reds = compute_reds(nant, antpos,tol=.1)
    # XXX haven't enforced xy = yx yet.  need to conjoin red groups for that
    ex_ants = [Antpol(i,nant).ant() for i in range(antpos.shape[0]) if antpos[i,0] < 0]
    kwargs['ex_ants'] = kwargs.get('ex_ants',[]) + ex_ants
    reds = filter_reds(reds, **kwargs)
    info = RedundantInfo(nant)
    info.init_from_reds(reds,antpos)
    return info


def redcal_Apkeys(data, info, xtalk=None, gains=None, vis=None,removedegen=False, uselogcal=True, maxiter=50, conv=1e-3, stepsize=.3, computeUBLFit=True, trust_period=1):
    _meta, _gain, _vis = omnical.calib.redcal(data, info, xtalk=None, gains=None, vis=None, removedegen=False, uselogcal=True, maxiter=50, conv=1e-3, stepsize=.3, computeUBLFit=True, trust_period=1)    
    meta, gain, vis, res = {}, {}, {}, {}
    mk_ap = lambda a: Antpol(a,NUMPOL[ant / info.nant], info.nant)
    for key in _meta.keys():
		if key == 'iter': 
			meta[key] = _meta[key]
			continue
		if key == 'res':
			for bl in _meta[key].keys():
				i,j = bl
				api = mk_ap(i)
				apj = mk_ap(j)
				res[(api,apj)] = _meta[key][bl]
			meta['res'] = res
			continue
		try: ant = int(key.split('chisq')[1])
		except(ValueError): meta[key] = _meta[key] #XXX this is due to a single array with key "chisq" i.e. no antnum associated
		ap = mk_ap(ant)
		meta['chisq'+str(ap)] = _meta[key] #XXX it might be worth making chisq a nested dictionary, with individual antpol keys
    for ant in _gain.keys(): gain[mk_ap(ant)] = _gain[ant]
    for bl in _vis.keys():
        i,j = bl
        api = mk_ap(i)
        apj = mk_ap(j)
        vis[(api,apj)] = _vis[bl]
    return meta, gain, vis

def redcal(data, info, xtalk=None, gains=None, vis=None,removedegen=False, uselogcal=True, maxiter=50, conv=1e-3, stepsize=.3, computeUBLFit=True, trust_period=1):
    _meta, _gain, _vis = omnical.calib.redcal(data, info, xtalk=None, gains=None, vis=None, removedegen=False, uselogcal=True, maxiter=50, conv=1e-3, stepsize=.3, computeUBLFit=True, trust_period=1)    
    meta, gain, vis, res = {}, {}, {}, {}
    mk_ap = lambda a: Antpol(a,NUMPOL[ant / info.nant], info.nant)
    Npols = info.nant/max(_gain.keys())
    for pi in range(Npols): #gains and chisqare per Antpol
    	gains[NUMPOL[pi]] = {} 
    	meta[NUMPOL[pi]] = {}
    	for pj in range(Npols): #res and vis per Antpol pair
    		vis[NUMPOL[pi]+NUMPOL[pj]] = {}
    		res[NUMPOL[pi]+NUMPOL[pj]] = {} 
    for key in _meta.keys():
    	if key=='iter': #XXX what is iter there for?
    		meta[key] = _meta[key]
    		continue
    	if key=='chisq': #XXX what is this there for (no ant index)?
    		meta[key] = _meta[key]
    		continue
    	if key=='res':
    		for i,j in _meta[key].keys():
    			api,apj = mk_ap(i),mk_ap(j)
    			res[api.pol()+apj.pol()][(api.ant(),apj.ant())] = _meta[key][(i,j)]
    		meta['res'] = res
    		continue
    	ant = int(key.split('chisq')[1])
    	ap = mk_ap(ant)
    	meta[ap.pol()]['chisq'+str(ap.ant())] = _meta[key]
    for ant in _gain.keys():
    	ap = mk_ap(ant)
    	gain[ap.pol()][ap.ant()] = _gain[key]
    for i,j in _vis.keys(): #XXX this could happen in the "res" loop above
    	api,apj = mk_ap(i),mk_ap(j)
    	vis[api.pol()+apj.pol()][(api.ant(),apj.ant())] = _vis[(i,j)]
    return meta, gain, vis

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
