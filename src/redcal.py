import capo
import capo.linsolve as ls
import numpy as np
from copy import deepcopy

def sim_red_data(reds, pols, gains=None, stokes_v=True, shape=(10,10), gain_scatter=.1):
    data = {}
    if not stokes_v: assert('xy' in pols and 'yx' in pols)
    bls = reduce(lambda x,y: x+y, reds)
    ants = set(reduce(lambda x,y: x+y, bls))
    if gains is None: gains = {}
    else: gains = deepcopy(gains)
    for ai in ants:
        for pol in pols:
            gains[(ai,pol[0])] = gains.get((ai,pol[0]), 1+gain_scatter*capo.oqe.noise((1,))) * np.ones(shape,dtype=np.complex)
            gains[(ai,pol[1])] = gains.get((ai,pol[1]), 1+gain_scatter*capo.oqe.noise((1,))) * np.ones(shape,dtype=np.complex)
    for bls in reds:
        for pol in pols:
            vis = capo.oqe.noise(shape)
            for bl in bls:
                data[bl+(pol,)] = vis * gains[(bl[0],pol[0])] * gains[(bl[1],pol[1])].conj()
                if not stokes_v:
                    data[bl+(pol[::-1],)] = vis * gains[(bl[0],pol[1])] * gains[(bl[1],pol[0])].conj()
    return gains, data

class RedundantInfo:
    def __init__(self, reds):
        self.reds = reds
    def build_eqs(self, bls, pols, stokes_v=True):
        eqs = {}
        for ubl, blgrp in enumerate(self.reds):
            blgrp = set(blgrp)
            for i,j in blgrp.intersection(bls):
                for pi,pj in pols:
                    eqs[self.pack_eqs_key(i,pi,j,pj,ubl,stokes_v=stokes_v)] = (i,j,pi+pj)
        return eqs
    def _solver(self, solver, data, wgts={}, stokes_v=True, detrend_phs=False, sparse=False, **kwargs):
        dc = capo.metrics.DataContainer(data)
        eqs = self.build_eqs(dc.bls(), dc.pols(), stokes_v=stokes_v)
        self.phs_avg = {} # detrend phases within redundant group, used for logcal to avoid phase wraps
        if detrend_phs:
            for ubl, blgrp in enumerate(self.reds):
                for pol in dc.pols():
                    self.phs_avg[(blgrp[0],pol)] = np.exp(-1j*np.median(np.unwrap([np.log(dc[bl+(pol,)]).imag for bl in blgrp],axis=0), axis=0))
                    for bl in blgrp: self.phs_avg[bl+(pol,)] = self.phs_avg[(blgrp[0],pol)]
        d_ls,w_ls = {}, {}
        for eq,key in eqs.items():
            d_ls[eq] = dc[key] * self.phs_avg.get(key,1)
        if len(wgts) > 0:
            wc = capo.metrics.DataContainer(wgts)
            for eq,key in eqs.items(): w_ls[eq] = wc[key]
        return solver(data=d_ls, wgts=w_ls, sparse=sparse, **kwargs)
    def pack_eqs_key(self, ant_i, pol_i, ant_j, pol_j, ubl_num, stokes_v=True):
        if stokes_v: pol = pol_i + pol_j
        else: pol = ''.join(sorted([pol_i,pol_j])) # make xy and yx the same
        return 'g%d%s * g%d%s_ * u%d%s' % (ant_i,pol_i,ant_j,pol_j,ubl_num,pol)
    def unpack_sol_key(self, k):
        if k.startswith('g'): # 'g' = gain solution
            return (int(k[1:-1]),k[-1])
        else: # 'u' = unique baseline solution
            return self.reds[int(k[1:-2])][0]+ (k[-2:],) 
    def pack_sol_key(self, k):
        if len(k) == 2: # 'g' = gain solution
            return 'g%d%s' % k
        else: # 'u' = unique baseline solution
            ubl_num = [cnt for cnt,blgrp in enumerate(self.reds) if blgrp[0] == k[:2]][0]
            return 'u%d%s' % (ubl_num, k[-1])
    def compute_ubls(self, data, gain_sols):
        dc = capo.metrics.DataContainer(data)
        ubl_sols = {}
        for ubl, blgrp in enumerate(self.reds):
            for pol in dc.pols():
                #d_gp = [dc[bl+(pol,)] / (gain_sols[(bl[0],pol[0])] * gain_sols[(bl[1],pol[1])].conj()) for bl in blgrp]
                d_gp = [dc[bl+(pol,)] for bl in blgrp]
                ubl_sols[blgrp[0]+(pol,)] = np.average(d_gp, axis=0) # XXX add option for median here?
        return ubl_sols
    def logcal(self, data, wgts={}, stokes_v=True, sparse=False):
        ls = self._solver(capo.linsolve.LogProductSolver, data, wgts=wgts, stokes_v=stokes_v, detrend_phs=True, sparse=sparse)
        sol = ls.solve()
        sol = {self.unpack_sol_key(k): sol[k] for k in sol.keys()}
        for ubl_key in [k for k in sol.keys() if len(k) == 3]:
            sol[ubl_key] = sol[ubl_key] * self.phs_avg[ubl_key].conj()
        return sol
    def lincal(self, data, sol0, wgts={}, stokes_v=True, sparse=False, conv_crit=1e-10, maxiter=50): # XXX use itersolve eventually
        #sol0 = dict(zip([self.pack_sol_key(k) for k in sol0.keys()],sol0.values()))
        sol0 = {self.pack_sol_key(k):sol0[k] for k in sol0.keys()}
        ls = self._solver(capo.linsolve.LinProductSolver, data, sol0=sol0, wgts=wgts, 
                stokes_v=stokes_v, sparse=sparse)
        meta, sol = ls.solve_iteratively(conv_crit=conv_crit, maxiter=maxiter)
        return meta, {self.unpack_sol_key(k):sol[k] for k in sol.keys()}

        
