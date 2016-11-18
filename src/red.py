'''
Tools for dealing with redundant array configurations.
'''

import numpy as n
from aipy.miriad import ij2bl, bl2ij
import aipy as a
import pylab as p
import time 
import multiprocessing as mpr

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
                    #bl = ij2bl(i,j)
                    bl = (i,j)
                    bls[sep] = bls.get(sep,[]) + [bl]
                    conj[bl] = c
    return bls, conj

def redundant_bl_cal(d1, w1, d2, w2, fqs, use_offset=False, maxiter=10, window='blackman-harris',
        clean=1e-4, verbose=False, tau=0., off=0.):
    '''Return gain and phase difference between two redundant measurements
    d1,d2 with respective weights w1,w2.'''
    # Compute measured values
    dtau,doff,mx = 0,0,0
    d12 = d2 * n.conj(d1)
    # For 2D arrays, assume first axis is time and integrate over it
    if d12.ndim > 1: d12_sum,d12_wgt = n.sum(d12,axis=0), n.sum(w1*w2,axis=0)
    else: d12_sum,d12_wgt = d12, w1*w2
    if n.all(d12_wgt == 0): return n.zeros_like(d12_sum), 0.
    d11 = d1 * n.conj(d1)
    if d11.ndim > 1: d11_sum,d11_wgt = n.sum(d11,axis=0), n.sum(w1*w1,axis=0)
    else: d11_sum,d11_wgt = d11, w1*w1
    window = a.dsp.gen_window(d12_sum.size, window=window)
    dlys = n.fft.fftfreq(fqs.size, fqs[1]-fqs[0])
    # Begin at the beginning
    d12_sum *= n.exp(-2j*n.pi*(fqs*tau+off))
    for j in range(maxiter):
        d12_sum *= n.exp(-2j*n.pi*(fqs*dtau+doff))
        tau += dtau; off += doff
        _phs = n.fft.fft(window*d12_sum)
        _wgt = n.fft.fft(window*d12_wgt)
        _phs,info = a.deconv.clean(_phs, _wgt, tol=clean)
        #_phs += info['res'] / a.img.beam_gain(_wgt)
        _phs = n.abs(_phs)
        mx = n.argmax(_phs)
        if j > maxiter/2 and mx == 0: # Fine-tune calibration with linear fit
            valid = n.where(d12_wgt > d12_wgt.max()/2, 1, 0)
            valid *= n.where(n.abs(d12_sum) > 0, 1, 0) # Throw out zeros, which NaN in the log below
            fqs_val = fqs.compress(valid)
            dly = n.real(n.log(d12_sum.compress(valid))/(2j*n.pi)) # This doesn't weight data
            wgt = d12_wgt.compress(valid); wgt.shape = (wgt.size,1)
            B = n.zeros((fqs_val.size,1)); B[:,0] = dly
            if use_offset: # allow for an offset component
                A = n.zeros((fqs_val.size,2)); A[:,0] = fqs_val; A[:,1] = 1
                dtau,doff = n.linalg.lstsq(A*wgt,B*wgt)[0].flatten()
            else:
                #A = n.zeros((fqs_val.size,1)); A[:,0] = fqs_val
                #dtau = n.linalg.lstsq(A*wgt,B*wgt)[0].flatten()[0]
                dtau = n.sum(wgt.flatten()*dly/fqs_val) / n.sum(wgt.flatten())
        else: # Pull out an integral number of phase wraps
            if mx > _phs.size/2: mx -= _phs.size
            dtau,doff = mx / (fqs[-1] - fqs[0]), 0
            mxs = mx + n.array([-1,0,1])
            dtau = n.sum(_phs[mxs] * dlys[mxs]) / n.sum(_phs[mxs])
            #dtau = n.sum(_phs**2 * dlys) / n.sum(_phs**2)
            #dtau = n.sum(_phs * dlys) / n.sum(_phs)
        if verbose: print j, dtau, doff, (tau, off), mx
        #P.subplot(211); P.plot(n.fft.fftshift(dlys), n.fft.fftshift(_phs)); P.xlim(-200,200)
        #P.subplot(212); P.plot(fqs, n.angle(d12_sum))
        #P.show()
    #import pylab as P
    #P.subplot(211); P.plot(n.fft.fftshift(dlys), n.fft.fftshift(_phs)); P.xlim(-200,200)
    #P.subplot(212); P.plot(fqs, n.angle(d12_sum))
    #P.show()
    off %= 1
    info = {'dtau':dtau, 'doff':doff, 'mx':mx} # Some information about last step, useful for detecting screwups
    g12 = d12_sum / d12_wgt.clip(1,n.Inf)
    g11 = d11_sum / d11_wgt.clip(1,n.Inf)
    gain = n.where(g11 != 0, g12/g11, 0)
    if use_offset: return gain, (tau,off), info
    else: return gain, tau, info

def fit_line(phs, fqs, valid, offset=False):
    fqs = fqs.compress(valid)
    dly = phs.compress(valid)
    B = n.zeros((fqs.size,1)); B[:,0] = dly
    if offset:
        A = n.zeros((fqs.size,2)); A[:,0] = fqs*2*n.pi; A[:,1] = 1
    else:
        A = n.zeros((fqs.size,1)); A[:,0] = fqs*2*n.pi#; A[:,1] = 1
    try:
        if offset:
            dt,off = n.linalg.lstsq(A,B)[0].flatten()
        else:
            dt = n.linalg.lstsq(A,B)[0][0][0]
            off = 0.0
        return dt,off
    except ValueError:
        import IPython 
        IPython.embed()

def mpr_clean(args):
    return a.deconv.clean(*args,tol=1e-4)[0]

def redundant_bl_cal_simple(d1,w1,d2,w2,fqs, cleantol=1e-4, window='none', finetune=True, verbose=False, plot=False, noclean=True, average=False, offset=False):
    '''Gets the phase differnce between two baselines by using the fourier transform and a linear fit to that residual slop. 
        Parameters
        ----------
        d1,d2 : NXM numpy arrays.
            Data arrays to find phase difference between. First axis is time, second axis is frequency. 
        w1,w2 : NXM numpy arrays. 
            corrsponding data weight arrays.
        fqs   : 1XM numpy array
            Array of frequencies in GHz.
        window: str
            Name of window function to use in fourier transform. Default is 'none'.
        finetune  : boolean 
            Flag if you want to fine tune the phase fit with a linear fit. Default is true.
        verbose: boolean
            Be verobse. Default is False.
        plot: boolean
            Turn on low level plotting of phase ratios. Default is False.
        cleantol: float
            Clean tolerance level. Default is 1e-4.
        noclean: boolean
            Don't apply clean deconvolution to remove the sampling function (weights).
        average: boolean
            Average the data in time before applying analysis. collapses NXM -> 1XM.
        offset: boolean
            Solve for a offset component in addition to a delay component.

        Returns
        -------
        delays : ndarrays
            Array of delays (if average == False), or single delay.    
        offsets: ndarrays
            Array of offsets (if average == False), or single delay.
    
'''
    d12 = d2 * n.conj(d1)
    # For 2D arrays, assume first axis is time. 
    if average:
        if d12.ndim > 1:
            d12_sum = n.sum(d12,axis=0).reshape(1,-1)
            d12_wgt = n.sum(w1*w1,axis=0).reshape(1,-1)
        else:
            d12_sum = d12.reshape(1,-1)
            d12_wgt = w1.reshape(1,-1)*w2.reshape(1,-1)
    else:
        d12_sum = d12
        d12_wgt = w1*w2
    #normalize data to maximum so that we minimize fft articats from RFI
    d12_sum *= d12_wgt
    d12_sum = d12_sum/n.where(n.abs(d12_sum)==0., 1., n.abs(d12_sum)) 
    #import IPython; IPython.embed()
    #exit()
    window = a.dsp.gen_window(d12_sum[0,:].size, window=window)
    dlys = n.fft.fftfreq(fqs.size, fqs[1]-fqs[0])
    # FFT and deconvolve the weights to get the phs
    _phs = n.fft.fft(window*d12_sum,axis=-1)
    _wgt = n.fft.fft(window*d12_wgt,axis=-1)
    _phss = n.zeros_like(_phs)
    if not noclean and average:
        _phss = a.deconv.clean(_phs, _wgt, tol=cleantol)[0]
    elif not noclean:
        pool=mpr.Pool(processes=4)
        _phss = pool.map(mpr_clean, zip(_phs,_wgt))
    else:
        _phss = _phs
    #print('Cleaning is taking %f seconds'%(t2-t1))
    _phss = n.abs(_phss)
    #get bin of phase
    mxs = n.argmax(_phss, axis=-1)
    #Fine tune with linear fit.
    mxs[mxs>_phss.shape[-1]/2] -= _phss.shape[-1]
    dtau = mxs / (fqs[-1] - fqs[0])
    mxs = n.dot(mxs.reshape(len(mxs),1),n.ones((1,3),dtype=int)) + n.array([-1,0,1])
    #import IPython; IPython.embed()
    taus = n.sum(_phss[n.arange(mxs.shape[0],dtype=int).reshape(-1,1),mxs] * dlys[mxs],axis=-1) / n.sum(_phss[n.arange(mxs.shape[0]).reshape(-1,1),mxs],axis=-1)
    dts,offs = [],[]
    if finetune:
    #loop over the linear fits
        t1 = time.time()
        for ii,(tau,d) in enumerate(zip(taus,d12_sum)):
            valid = n.where(d != 0, 1, 0) # Throw out zeros, which NaN in the log below
            #print valid.any()
            valid = n.logical_and(valid, n.logical_and(fqs>.11,fqs<.19))
            dly = n.angle(d*n.exp(-2j*n.pi*tau*fqs))
            dt,off = fit_line(dly,fqs,valid,offset=offset)
            dts.append(dt), offs.append(off)
            if plot:
                p.subplot(411)
                p.plot(fqs,n.angle(d12_sum[ii]), linewidth=2)
                p.plot(fqs,d12_sum[ii], linewidth=2)
                p.plot(fqs, n.exp((2j*n.pi*fqs*(tau+dt))+off))
                p.hlines(n.pi, .1,.2,linestyles='--',colors='k')
                p.hlines(-n.pi, .1,.2,linestyles='--',colors='k')
                p.subplot(412)
                p.plot(fqs,n.unwrap(dly)+2*n.pi*tau*fqs, linewidth=2)
                p.plot(fqs,dly+2*n.pi*tau*fqs, linewidth=2,ls='--')
                p.plot(fqs,2*n.pi*tau*fqs, linewidth=2,ls='-.')
                p.plot(fqs,2*n.pi*(tau+dt)*fqs + off, linewidth=2,ls=':')
                p.subplot(413)
                p.plot(dlys, n.abs(_phs[ii]),'-.')
                p.xlim(-400,400)
                p.subplot(414)
                p.plot(fqs,dly, linewidth=2)
                p.plot(fqs,off+dt*fqs*2*n.pi, '--')
                p.hlines(n.pi, .1,.2,linestyles='--',colors='k')
                p.hlines(-n.pi, .1,.2,linestyles='--',colors='k')
                print 'tau=', tau
                print 'tau + dt=', tau+dt
                p.xlabel('Frequency (GHz)', fontsize='large')
                p.ylabel('Phase (radians)', fontsize='large')
        p.show()

        #print('Linear Fitting is taking %f seconds'%(time.time()-t1))
        dts = n.array(dts)
        offs = n.array(offs)

           # Pull out an integral number of phase wraps
    #if verbose: print tau, dtau, mxs, dt, off
    info = {'dtau':dts, 'off':offs, 'mx':mxs} # Some information about last step, useful for detecting screwups
    if verbose: print info, taus, taus+dts, offs
    return taus+dts,offs




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
        if conj: i,j = j,i
        s = self.sep_index(sep)
        i,j = self.antpol_index(i,pol), self.antpol_index(j,pol)
        amp_line[i], amp_line[j], amp_line[s] = 1, 1, 1
        self.M_amp.append(amp_line)
        phs_line[i], phs_line[j], phs_line[s] = -1, 1, 1
        self.M_phs.append(phs_line)
        self.meas_order.append((bl,pol,sep,conj))
    #def invert(self, amp_meas, phs_meas):
    def invert(self, logvis, einstr='pm,mq'):
        '''Return the antenna gain solutions and sky/sep solutions for the provided set of amp/phs measurements.
        phs_meas needs to have had any necessary conjugation applied already'''
        # Invert the matrices recording which parameters are involved in each measurement
        # M_phs is (nmeas,nprm), so M_phs_inv is (nprm,nmeas)
        if self.M_phs_inv is None: self.M_phs_inv = n.linalg.pinv(n.array(self.M_phs))
        if self.M_amp_inv is None: self.M_amp_inv = n.linalg.pinv(n.array(self.M_amp))
        #logvis = n.log(vis)
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

class RelativeLogCalMatrix(LogCalMatrix):
    def __init__(self, nants, npols):
        LogCalMatrix.__init__(self, nants, npols, 0)
        self.constraint = []
        self.constraint_val = []
        self.constraint_wgt = []
    # xXX sep_index shouldn't exist
    def add_meas_record(self, bl1, pol1, conj1, bl2, pol2, conj2):
        #amp_line = n.zeros((self.nprms,), dtype=n.double)
        phs_line = n.zeros((self.nprms,), dtype=n.double)
        i1,j1 = bl2ij(bl1)
        i2,j2 = bl2ij(bl2)
        if conj1: i1,j1 = j1,i1
        if conj2: i2,j2 = j2,i2
        i1,j1 = self.antpol_index(i1,pol1), self.antpol_index(j1,pol1)
        i2,j2 = self.antpol_index(i2,pol2), self.antpol_index(j2,pol2)
        #amp_line[j2], amp_line[i2], amp_line[j1], amp_line[i1] = 1, 1, -1, -1
        #self.M_amp.append(amp_line)
        #phs_line[j2], phs_line[i2], phs_line[j1], phs_line[i1] = 1, -1, -1, 1
        phs_line[j2] += 1; phs_line[i2] += -1; phs_line[j1] += -1; phs_line[i1] += 1
        self.M_phs.append(phs_line)
        self.meas_order.append((bl1,pol1,conj1,bl2,pol2,conj2))
    def add_constraint(self, ant, pol, val, wgt=1000.):
        '''Additional constraints is a list of (ant,pol,val) constraints that are priced into the fit with the
        provided weight.'''
        line = n.zeros((self.nprms,), dtype=n.double)
        line[self.antpol_index(ant,pol)] = wgt
        self.constraint.append(line)
        self.constraint_val.append(val)
        self.constraint_wgt.append(wgt)
    def invert(self, dly, wgts=None, einstr='pm,m'):
        '''Return the antenna gain solutions and sky/sep solutions for the provided set of amp/phs measurements.
        phs_meas needs to have had any necessary conjugation applied already.'''
        # Invert the matrices recording which parameters are involved in each measurement
        # M_phs is (nmeas,nprm), so M_phs_inv is (nprm,nmeas)
        #if self.M_phs_inv is None:
        #    self.M_phs_inv = n.linalg.pinv(n.array(self.constraint + self.M_phs))
        if wgts is None: wgts = n.ones_like(dly)
        wgts = n.concatenate([self.constraint_wgt, wgts])
        #wgts.shape = (wgts.size, 1)
        M = n.array(self.constraint + self.M_phs)
        val = n.concatenate([self.constraint_val, dly])
        
        self.M_phs_inv = n.linalg.pinv(M * wgts.reshape((wgts.size,1)))
        #if self.M_amp_inv is None: self.M_amp_inv = n.linalg.pinv(n.array(self.M_amp))
        dly = n.einsum(einstr, self.M_phs_inv, val * wgts)
        #dly = n.linalg.lstsq(n.array([phs_line1, phs_line2, phs_line3] + self.M_phs), dly)[0]
        #gain = n.einsum(einstr, self.M_amp_inv, gain)
        dly_sol = {}
        #gain_sol = {}
        for ant in self.ants:
            #dly_sol[ant], gain_sol[ant] = {}, {}
            dly_sol[ant] = {}
            for pol in self.pols:
                i = self.antpol_index(ant,pol)
                #gain_sol[ant][pol] = gain[i]
                dly_sol[ant][pol] = dly[i]
        return dly_sol#, gain_sol

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
        
