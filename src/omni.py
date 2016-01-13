import numpy as n, omnical

POLNUM = {
    'x':0, # factor to multiply ant index for internal ordering
    'y':1,
    'l':2,
    'r':3,
    'a':4,
    'b':5,
}

NUMPOL = dict(zip(POLNUM.values(), POLNUM.keys()))
    

def antpol2ind(ant, pol, nant):
    return POLNUM[pol] * nant + ant

def ind2antpol(i, nant):
    return i % nant, NUMPOL[i / nant]

# XXX filter_reds w/ pol support should probably be in omnical
def filter_reds(nant, reds, bls=None, ex_bls=None, ants=None, ex_ants=None, ubls=None, ex_ubls=None, crosspols=None, ex_crosspols=None):
    '''Filter redundancies to include/exclude the specified bls, antennas, and unique bl groups and polarizations.'''
    def m(i): return i % nant
    pols = [NUMPOL[r[0][0]/nant]+NUMPOL[r[0][1]/nant] for r in reds]
    if crosspols:
        reds = [r for r,p in zip(reds,pols) if p in crosspols]
        pols = [p for p in pols if p in crosspols]
    if ex_crosspols:
        reds = [r for r,p in zip(reds,pols) if not p in ex_crosspols]
        pols = [p for p in pols if not p in ex_crosspols]
    upols = {}
    for p in pols: upols[p] = None
    upols = upols.keys()
    if ubls or ex_ubls:
        bl2gp = {}
        for i,gp in enumerate(reds):
            for bl in gp:
                blmod = m(bl[0]), m(bl[1])
                bl2gp[blmod] = bl2gp[blmod[::-1]] = bl2gp.get(blmod,[]) + [i]
        if ubls:
            ubls = [bl2gp[bl] for bl in ubls if bl2gp.has_key(bl)]
            ubls = reduce(lambda x,y: x+y, ubls)
        else: ubls = range(len(reds))
        if ex_ubls:
            ex_ubls = [bl2gp[bl] for bl in ex_ubls if bl2gp.has_key(bl)]
            ex_ubls = reduce(lambda x,y: x+y, ex_ubls)
        else: ex_ubls = []
        reds = [gp for i,gp in enumerate(reds) if i in ubls and i not in ex_ubls]
    if bls is None: bls = [bl for gp in reds for bl in gp]
    else:
        _bls = []
        for pi,pj in upols:
            di,dj = POLNUM[pi]*nant, POLNUM[pj]*nant
            _bls += [(i+di,j+dj) for i,j in bls]
        bls = _bls
    if ex_bls: bls = [(i,j) for i,j in bls 
            if (m(i),m(j)) not in ex_bls and (m(j),m(i)) not in ex_bls]
    if ants: bls = [(i,j) for i,j in bls if m(i) in ants and m(j) in ants]
    if ex_ants: bls = [(i,j) for i,j in bls 
            if m(i) not in ex_ants and m(j) not in ex_ants]
    bld = {}
    for bl in bls: bld[bl] = bld[bl[::-1]] = None
    reds = [[bl for bl in gp if bld.has_key(bl)] for gp in reds]
    return [gp for gp in reds if len(gp) > 1]
        
def aa_to_info(aa, pols=['x'], **kwargs):
    '''Use aa.ant_layout to generate redundances based on ideal placement.
    The remaining arguments are passed to omnical.arrayinfo.filter_reds()'''
    layout = aa.ant_layout
    antpos = -n.ones((len(aa)*len(pols),3)) # -1 to flag unused antennas
    xs,ys = n.indices(layout.shape)
    for ant,x,y in zip(layout.flatten(), xs.flatten(), ys.flatten()):
        for z,pol in enumerate(pols):
            z = 2**z # exponential ensures diff xpols aren't redundant w/ each other
            ind = antpol2ind(ant,pol,len(aa))
            antpos[i,0],antpos[i,1],antpos[i,2] = x,y,z
    reds = omnical.arrayinfo.compute_reds(antpos,tol=.1)
    # XXX haven't enforced xy = yx yet.  need to conjoin red groups for that
    ex_ants = [i for i in range(antpos.shape[0]) if antpos[i,0] < 0]
    kwargs['ex_ants'] = kwargs.get('ex_ants',[]) + ex_ants
    #reds = omnical.arrayinfo.filter_reds(reds, **kwargs)
    reds = filter_reds(reds, **kwargs)
    info = omnical.info.RedundantInfo()
    info.init_from_reds(reds,antpos)
    return info

def compute_xtalk(res, wgts):
    '''Estimate xtalk as time-average of omnical residuals.'''
    xtalk = {}
    for key in res:
        r,w = n.where(wgts[key] > 0, res[key], 0), wgts[key].sum(axis=0)
        w = n.where(w == 0, 1, w)
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
    n.savez(filename,**d)

def from_npz(filename, meta={}, gains={}, vismdl={}, xtalk={}, jds=[], lsts=[], freqs=[]):
    '''Reconstitute results from to_npz, returns meta, gains, vismdl, xtalk, each
    keyed first by polarization, and then by bl/ant/keyword.'''
    npz = n.load(filename)
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
