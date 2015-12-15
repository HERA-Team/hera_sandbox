import numpy as n, omnical

def aa_to_info(aa, **kwargs):
    '''Use aa.ant_layout to generate redundances based on ideal placement.
    The remaining arguments are passed to omnical.arrayinfo.filter_reds()'''
    layout = aa.ant_layout
    antpos = -n.ones((len(aa),3)) # -1 to flag unused antennas
    xs,ys = n.indices(layout.shape)
    for i,x,y in zip(layout.flatten(), xs.flatten(), ys.flatten()):
        antpos[i,0],antpos[i,1] = x,y
    reds = omnical.arrayinfo.compute_reds(antpos,tol=.1)
    ex_ants = [i for i in range(antpos.shape[0]) if antpos[i,0] < 0]
    kwargs['ex_ants'] = kwargs.get('ex_ants',[]) + ex_ants
    reds = omnical.arrayinfo.filter_reds(reds, **kwargs)
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

def to_npz(filename, meta, gains, vismdl, xtalk, jds, lsts, freqs):
    '''Write results from omnical.calib.redcal (meta,gains,vismdl) to npz file.
    Each of these is assumed to be a dict keyed by pol, and then by bl/ant/keyword'''
    d = {}
    for pol in meta:
        for k in meta[pol]:
            if k.startswith('chisq'):
                d[k + ' %s' % pol] = meta[pol][k]
        # XXX chisq after xtalk removal
    for pol in gains:
        for ant in gains[pol]:
            d['%d%s' % (ant,pol)] = gains[pol][ant].conj() # conj to miriad 
    for pol in vismdl:
        for bl in vismdl[pol]:
            d['<%d,%d> %s' % (bl[0],bl[1],pol)] = vismdl[pol][bl].conj() # conj to miriad 
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
    for k in [f for f in npz.files if f.startswith('j')]: 
        jds = npz[k] 
    for k in [f for f in npz.files if f.startswith('l')]: 
        lsts = npz[k]
    for k in [f for f in npz.files if f.startswith('f')]:
        freqs = npz[k]
    return meta, gains, vismdl, xtalk, jds, lsts, freqs
