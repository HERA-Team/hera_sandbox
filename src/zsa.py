import aipy as a, numpy as n, pylab as p
import sys, scipy
import capo.omni as omni


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
    #p.plot(d12_sum)
    #p.show()
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

def noise(size):
    #generates a complex random gaussian noise with std=1 and mean=0.
    sig = 1./n.sqrt(2)
    return n.random.normal(scale=sig, size=size) + 1j*n.random.normal(scale=sig, size=size)


def grid2ij(GRID):
    '''
        bl_str = given sep, returns bls in string format.
        bl_conj = given a baseline (miriad bl) gives separation.
        bl2sep_str = given baseline (miriad) return its separation.    
    '''
    bls, conj = {}, {}
    for ri in range(GRID.shape[0]):
        for ci in range(GRID.shape[1]):
            for rj in range(GRID.shape[0]):
                for cj in range(GRID.shape[1]):
                    if ci > cj: continue
#                    if ri > rj and ci == cj: continue
#                    if ci > cj and ri == rj: continue
                    sep = (rj-ri, cj-ci)
                    sep = '%d,%d'%sep
                    i,j = GRID[ri, ci], GRID[rj,cj]
                    bls[sep] = bls.get(sep,[]) + [(i,j)]
    for sep in bls.keys():
        if sep == '0,0' or len(bls[sep]) < 2 or (sep[-1] == '0' and sep[0] == '-'): del(bls[sep])
    for sep in bls:
        conj[sep] = [i>j for i,j in bls[sep]]

    bl_str,bl_conj,bl2sep_str = {}, {}, {}
    for sep in bls:
        bl_str[sep],bl_list = [], []
        for (i,j),c in zip(bls[sep],conj[sep]):
            if c: i,j = j,i
            bl_list.append(a.miriad.ij2bl(i,j))
            bl_str[sep].append('%d_%d'%(i,j))
            bl2sep_str[a.miriad.ij2bl(i,j)] = bl2sep_str.get(a.miriad.ij2bl(i,j),'') + sep
            bl_conj[a.miriad.ij2bl(i,j)] = c
        bls[sep] = bl_list
        bl_str[sep] = ','.join(bl_str[sep])
    return bl_str,bl_conj,bl2sep_str

def get_dict_of_uv_data(filenames, antstr, polstr, decimate=1, decphs=0, verbose=False, recast_as_array=True, return_lsts=False):
    lsts, times, dat, flg = [], [], {}, {}
    if type(filenames) == 'str': filenames = [filenames]
    for filename in filenames:
        if verbose: print '   Reading', filename
        uv = a.miriad.UV(filename)
        a.scripting.uv_selector(uv, antstr, polstr)
        if decimate > 1: uv.select('decimate', decimate, decphs)
        for (crd,t,(i,j)),d,f in uv.all(raw=True):
            if len(times) == 0 or t != times[-1]:
                times.append(t)
                lsts.append(uv['lst'])
            bl = a.miriad.ij2bl(i,j)
            if not dat.has_key(bl): dat[bl],flg[bl] = {},{}
            pol = a.miriad.pol2str[uv['pol']]
            if not dat[bl].has_key(pol):
                dat[bl][pol],flg[bl][pol] = [],[]
            dat[bl][pol].append(d)
            flg[bl][pol].append(f)
    if recast_as_array:
        # This option helps reduce memory footprint, but it shouldn't
        # be necessary: the replace below should free RAM as quickly
        # as it is allocated.  Unfortunately, it doesn't seem to...
        for bl in dat.keys():
          for pol in dat[bl].keys():
            dat[bl][pol] = n.array(dat[bl][pol])
            flg[bl][pol] = n.array(flg[bl][pol])
    if return_lsts: times = lsts
    return n.array(times), dat, flg

def list2str(li):
    '''Take list of baselines and convert to string format for plot_uv'''
    s = ''
    for i in li:
        s += '_'.join(map(str,i)) + ','
    return s[:-1]

def flag_by_chisq(filenames, nsig=12, deg=8, outfile=False):
    '''Use the omnical global chisq to flag the model visibilities.'''
    m,g,v,x = omni.from_npz(filenames)
    f = {}
    for pol in g.keys():
        if not pol in f: f[pol] = {}
        for k in g[pol].keys():        
            chisq = m['chisq'+str(k)+pol[0]]
            #iterate twice on flattening bandpass to find rfi
            mask = n.zeros_like(chisq, dtype=n.bool)
            #Run loop twice to get better fit after removing large rfi
            for i in range(2):
                wgts = n.logical_not(mask)
                chisq *= wgts
                med_chisq = n.median(chisq, axis=0) 
                w = n.median(chisq,axis=0)!=0. 
                fit = n.polyfit(n.arange(len(med_chisq))[w], n.log10(med_chisq[w]), deg=deg)
                flat_chisq = chisq/10**n.polyval(fit, n.arange(len(med_chisq)))
                med = n.median(flat_chisq)
                sig = n.sqrt(n.median(n.abs(flat_chisq-med)**2))
                mask |= n.where(flat_chisq > (med + nsig*sig), True, False)
                #import IPython; IPython.embed()
            f[pol][k] = n.logical_not(mask)#weights for the data
    if outfile:
        pass 
    return m,g,v,x,f

def flatten_reds(reds):
    'Take a list of lists and flattens it'
    freds = []
    for r in reds:
        freds += r
    return freds

def order_data(dd, info):
    '''Order data in dict, where pol is first key and bl tuple is second key, the same way an info object is oriented'''
    d = {}
    for bl in dd.keys():
        for pol in dd[bl].keys():
            if bl in info.bl_order(): 
                if not d.has_key(bl): d[bl] = {}
                d[bl][pol] = dd[bl][pol]
            else:
                if not d.has_key(bl[::-1]): d[bl[::-1]] = {}
                d[bl[::-1]][pol] = n.conj(dd[bl][pol])
    return d

def run_nb(workdir, fileroot, basenotebook, agdir=''):

    from shutil import copy
    from subprocess import call
    import os

    if agdir: #git directory
        os.environ['AGDIR'] = agdir

    nb = basenotebook
    os.chdir(workdir)
    copy(nb, '{0}/{1}.ipynb'.format(workdir, fileroot))

    cmd = 'jupyter nbconvert {0}.ipynb --inplace --execute --to notebook --allow-errors --ExecutePreprocessor.timeout=15000'.format(fileroot).split(' ')
    status = call(cmd)

    cmd = 'jupyter trust {0}.ipynb'.format(fileroot).split(' ')
    status = call(cmd)

    return True
