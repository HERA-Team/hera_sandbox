'''
Tools for dealing with redundant array configurations.
'''

import numpy as n
from aipy.miriad import ij2bl, bl2ij

def group_redundant_bls(antpos):
    '''Return 2 dicts: bls contains baselines grouped by separation ('drow,dcol'), conj indicates for each
    baseline whether it must be conjugated to be redundant with the rest of the baselines in its redundancy group.'''
    bls,conj = {}, {}
    for ri in xrange(antpos.shape[0]):
        for ci in xrange(antpos.shape[1]):
            for rj in xrange(antpos.shape[0]):
                for cj in xrange(ci,antpos.shape[1]):
                    if ri >= rj and ci == cj: continue # exclude repeat +/- listings of certain bls
                    sep = '%d,%d' % (rj-ri, cj-ci)
                    i,j = antpos[ri,ci], antpos[rj,cj]
                    if i > j: i,j,c = j,i,True
                    else: c = False
                    bl = ij2bl(i,j)
                    bls[sep] = bls.get(sep,[]) + [bl]
                    conj[bl] = c
    return bls, conj

class LogCalMatrix:
    def __init__(self, nants, npols, nseps):
        self.nants = nants
        self.npols = npols
        self.nseps = nseps
        self.nprms = nants * npols + nseps # total number of parameters being solved for
        self.M_phs = []
        self.M_amp = []
        self.M_phs_inv = None
        self.M_amp_inv = None
        self.seps = {}
        self.pols = {}
        self.ants = {}
        self.meas_order = []
    def sep_index(self, sep):
        try: return self.nants * self.npols + self.seps[sep]
        except:
            assert(len(self.seps) < self.nseps)
            self.seps[sep] = len(self.seps)
            return self.nants * self.npols + self.seps[sep]
    def antpol_index(self, i, pol):
        # Add this pol to the list of pols we're solving for
        try: p = self.pols[pol]
        except:
            assert(len(self.pols) < self.npols)
            self.pols[pol] = len(self.pols)
            p = self.pols[pol]
        # Add this ant to the list of pols we're solving for
        try: a = self.ants[i]
        except:
            assert(len(self.ants) < self.nants)
            self.ants[i] = len(self.ants)
            a = self.ants[i]
        return self.npols * a + p
    def add_meas_record(self, bl, pol, sep, conj):
        amp_line = n.zeros((self.nprms,), dtype=n.double)
        phs_line = n.zeros((self.nprms,), dtype=n.double)
        i,j = bl2ij(bl)
        s = self.sep_index(sep)
        i,j = self.antpol_index(i,pol), self.antpol_index(j,pol)
        amp_line[i], amp_line[j], amp_line[s] = 1, 1, 1
        self.M_amp.append(amp_line)
        if not conj: phs_line[i], phs_line[j], phs_line[s] = -1, 1, 1
        else: phs_line[i], phs_line[j], phs_line[s] = 1, -1, 1
        self.M_phs.append(phs_line)
        self.meas_order.append((bl,pol,sep,conj))
    #def invert(self, amp_meas, phs_meas):
    def invert(self, vis, einstr='pm,mq'):
        '''Return the antenna gain solutions and sky/sep solutions for the provided set of amp/phs measurements.
        phs_meas needs to have had any necessary conjugation applied already'''
        # Invert the matrices recording which parameters are involved in each measurement
        # M_phs is (nmeas,nprm), so M_phs_inv is (nprm,nmeas)
        if self.M_phs_inv is None: self.M_phs_inv = n.linalg.pinv(n.array(self.M_phs))
        if self.M_amp_inv is None: self.M_amp_inv = n.linalg.pinv(n.array(self.M_amp))
        logvis = n.log(vis)
        phs = n.einsum(einstr, self.M_phs_inv, logvis.imag)
        amp = n.einsum(einstr, self.M_amp_inv, logvis.real)
        g = n.exp(amp + 1j * phs)
        i = self.nants * self.npols
        g_avg = n.sqrt(n.average(n.abs(g[:i])**2, axis=0))
        g[:i] /= n.where(g_avg > 0, g_avg, 1)
        g[i:] *= n.abs(g_avg)**2
        ant_sol = {}
        for ant in self.ants:
            ant_sol[ant] = {}
            for pol in self.pols:
                i = self.antpol_index(ant,pol)
                ant_sol[ant][pol] = g[i]
        sep_sol = {}
        for sep in self.seps:
            s = self.sep_index(sep)
            sep_sol[sep] = g[s]
        return ant_sol, sep_sol

def estimate_xtalk(ant_sol, vis, meas_order):
    dsum,dwgt = {},{}
    for (bl,pol,sep,cnj),d in zip(meas_order,vis):
        i,j = bl2ij(bl)
        gi,gj = ant_sol[i][pol], ant_sol[j][pol]
        if cnj: gij = gi * n.conj(gj)
        else: gij = gj * n.conj(gi)
        d /= gij
        dsum[sep] = dsum.get(sep,0) + d
        dwgt[sep] = dwgt.get(sep,0) + 1
    xtalk = []
    for (bl,pol,sep,cnj),d in zip(meas_order,vis):
        i,j = bl2ij(bl)
        gi,gj = ant_sol[i][pol], ant_sol[j][pol]
        if cnj: gij = gi * n.conj(gj)
        else: gij = gj * n.conj(gi)
        davg = dsum[sep] / dwgt[sep]
        davg *= gij
        xtalk.append(d - davg)
    return n.array(xtalk)

def xtalk_to_dict(xtalk, meas_order):
    xd = {}
    for (bl,pol,sep,cnj),x in zip(meas_order,xtalk):
        if cnj: x = x.conj()
        if not xd.has_key(bl): xd[bl] = {}
        xd[bl][pol] = x
    return xd
        
