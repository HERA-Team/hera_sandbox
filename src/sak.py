import numpy as np, pandas as pd, matplotlib.pyplot as plt, aipy as a
from mpl_toolkits.axes_grid1 import make_axes_locatable


def get_ants(mfile, pol='xx', search='f1', columns=['fengine', 'plate', 'rx', 'stations', 'ant']):
    '''Get antennas corresponding to search string. Given mfile.
       Search string can be 
            f{i} for fengine i : returns all antennas on a specific fengine. 
            r{i} for receiverator i : returns all antennas on a specific recieverator.
            f{i} and r{j} : overlap of fengine i and recieverator j.
            etc.. 
        Logicals are read from left to right.
            e.g. f2 and r3 or r3 = ((f2 and r3) or r3)
        
        returns x and y pol strings as a dictionary. 
    '''
    #get dataframe from map file. Make sure columns are in order below.
    mapping = {'f':'fengine', 'r': 'rx'}
    d = pd.read_csv(mfile, sep=' -> ', names=columns)

    splits = np.array(search.split(' '))
    merge_these = []
    for i in splits[::2]:
        #loop throguh f's and r's. Get where True.
        map_key = mapping[i[0]]
        merge_these.append( d[map_key].str.contains(i) )
    final_truth = merge_these[0]
    for i,k in enumerate(splits[1::2]):
        #loop through and's and or's. 
        if k=='and':
            final_truth = np.logical_and(final_truth, merge_these[i+1])
        elif k=='or':
            final_truth = np.logical_or(final_truth, merge_these[i+1])
        else:
            continue

    final_ants = d['ant'][final_truth].as_matrix()
    final_strings = {'x': [], 'y': []}
    for f in final_ants:
        if f.endswith('X'): final_strings['x'].append(f[1:-1])
        elif f.endswith('Y'): final_strings['y'].append(f[1:-1])
    for k in final_strings:
        final_strings[k] = ','.join(final_strings[k])

    return final_strings[pol[0]]

def data_mode(data, mode='abs'):
    if mode.startswith('phs'): data = np.angle(data)
    elif mode.startswith('lin'):
        data = np.absolute(data)
        #data = np.masked_less_equal(data, 0)
    elif mode.startswith('real'): data = data.real
    elif mode.startswith('imag'): data = data.imag
    elif mode.startswith('log'):
        data = np.absolute(data)
        data = np.log10(data)
    else: raise ValueError('Unrecognized plot mode.')
    return data

def waterfall(d, mode='log', mx=None, drng=None, recenter=False, **kwargs):
    if np.ma.isMaskedArray(d): d = d.filled(0)
    if recenter: d = a.img.recenter(d, np.array(d.shape)/2)
    d = data_mode(d, mode=mode)
    if mx is None: mx = d.max()
    if drng is None: drng = mx - d.min()
    mn = mx - drng
    return plt.imshow(d, vmax=mx, vmin=mn, aspect='auto', interpolation='nearest', **kwargs)

def full_data_view(d, f, mx=None, drng=None, save=None, clean=1e-3, verbose=False, **kwargs):
    flags = np.logical_not(f).astype(np.float)
    if np.ma.isMaskedArray(d): d = d.filled(0)
    if mx is None: mx = d.max()
    if drng is None: drng = mx - d.min()
    mn = mx - drng
    
    f,axarr = plt.subplots(1,4)
    plt.figure(num=1,figsize=(16,12))
    #absolute
    if verbose: print 'Plotting abs waterfall'
    im0 = axarr[0].imshow(np.log10(np.absolute(d)), aspect='auto', interpolation='nearest', vmax=mx, vmin=mn, **kwargs)
    axarr[0].set_xlabel('Frequency bin')
    axarr[0].set_ylabel('Integration')
    
    #angle
    if verbose: print 'Plotting phase waterfall'
    im1 = axarr[1].imshow(np.angle(d), aspect='auto', interpolation='nearest', vmax=np.pi, vmin=-1*np.pi, **kwargs)
    axarr[1].set_xlabel('Frequency bin')
    axarr[1].set_ylabel('Integration')
    
    #delay
    if verbose: print 'Plotting delay waterfall'
    w = a.dsp.gen_window(d.shape[-1], window='blackman-harris')
    _dw = np.fft.ifft(d*w)
    _ker= np.fft.ifft(flags*w)
    gain = a.img.beam_gain(_ker)
    for time in range(_dw.shape[0]):
        _dw[time,:],info = a.deconv.clean(_dw[time,:], _ker[time,:], tol=clean)
        _dw[time,:] += info['res']/gain
    dd = np.fft.fftshift(np.ma.array(_dw),axes=1)
    im2 = axarr[2].imshow(np.log10(np.absolute(dd)), aspect='auto', interpolation='nearest', vmax=mx, vmin=mn, **kwargs)
    axarr[2].set_xlabel('Delay bin')
    axarr[2].set_ylabel('Integration')
    
    #delay-rate
    if verbose: print 'Plotting delay-rate waterfall'
    w = a.dsp.gen_window(d.shape[0], window='blackman-harris')
    w.shape += (1,)
    wgts = np.where(d!= 0, 1., 0.) * w
    gain = np.sqrt(np.average(wgts**2, axis=0))
    _ker = np.fft.ifft(wgts, axis=0) # w already put in 2 lines above
    _dw = np.fft.ifft(d*w, axis=0)
    
    for chan in range(_dw.shape[-1]):
        if gain[chan] == 0: continue
        _dw[:,chan],info = a.deconv.clean(_dw[:,chan],_ker[:,chan],tol=clean)
        _dw[:,chan] += info['res'] / gain[chan]
    dr = np.ma.array(np.fft.fftshift(_dw, axes=0))
    im3 = axarr[3].imshow(np.log10(np.absolute(dr)), aspect='auto', interpolation='nearest', vmax=mx, vmin=mn, **kwargs)
    axarr[3].set_xlabel('Frequency')
    axarr[3].set_ylabel('Delay-rate bin')
        
    #attach colorbars
    imshows = [im0,im1,im2,im3]
    for i,ax in enumerate(axarr):
        div = make_axes_locatable(ax)
        cax = div.append_axes("right", size="10%", pad=0.05)
        cbar = plt.colorbar(imshows[i],cax=cax)
    plt.subplots_adjust(wspace=0.5,left=0.08,right=0.93,top=0.95,bottom=0.13)
    if not save is None: plt.savefig(save)
    
def plot_phase_ratios(data,ref=None):
    '''
    Plots ratios of baselines given in data. 
    Data is a nested dictionary. First key is baseline, second key is pol. 
    "ref" specifies an index  in the array of baselines to use as a reference.
    '''
    bls = data.keys()
    nbls = len(bls)
    pol = data[bls[0]].keys()[0]

    nratios = (nbls * (nbls-1))/2
    r = int(divmod(nratios,3)[0] + np.ceil(divmod(nratios,3)[1]/3.))
    c = 3
    ncross = []
    for k in range(nbls): 
        for i in range(k+1,nbls): 
            if ref is None: ncross.append((bls[k],bls[i]))
            else: ncross.append((bls[k],bls[ref]))
    
    fig = plt.figure(figsize=(16,12))
    for i,k in enumerate(ncross):
        ax = plt.subplot(r,c,i+1)
        plt.title(str(k),color='magenta')
        g = 1.0
        waterfall(data[k[0]][pol]*np.conj(data[k[-1]][pol])*g, mode='phs', cmap='jet', mx=np.pi, drng=2*np.pi)
        plt.grid(0)
        if divmod(i,c)[-1] != 0:  ax.yaxis.set_visible(False) 
        if divmod(i,c)[0] != r-1: ax.xaxis.set_visible(False)
    cax = fig.add_axes([0.2, 0.06, 0.6, 0.01])
    try: plt.colorbar(cax=cax, orientation='horizontal')
    except RuntimeError: pass

def order_data(dd, info):
    '''
    Order data in dict, where pol is first key and bl tuple is second key, the same way an info object is oriented
    '''
    d = {}
    for bl in dd.keys():
        for pol in dd[bl].keys():
            if bl in info.bl_order(): 
                if not d.has_key(bl): d[bl] = {}
                d[bl][pol] = dd[bl][pol]
            else:
                if not d.has_key(bl[::-1]): d[bl[::-1]] = {}
                d[bl[::-1]][pol] = np.conj(dd[bl][pol])
    return d
