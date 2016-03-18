import sys, capo, aipy, numpy as np, pylab, glob

def get_delayfall(filenamearray,antstr,pol,blo=0,bhi=203):
    timeinfo,d,f = capo.arp.get_dict_of_uv_data(filenamearray,antstr=antstr,polstr=pol)
    sp = np.array(antstr.split('_'))
    bl = (int(sp[0]),int(sp[1]))
    
    data = d[bl][pol][:,blo:bhi]
    flag = f[bl][pol][:,blo:bhi]
    
    def b2f(b): return 0.1 + 0.1*b/203.
    
    freqs = np.linspace(b2f(blo),b2f(bhi),bhi-blo) #XXX this is for full bandwidth
    
    assert(data.shape[1] == freqs.shape[0])
    
    delays = np.fft.fftfreq(freqs.size,freqs[1]-freqs[0])
    delays = np.fft.fftshift(delays)
    dtrans = []
    
    for t in range(data.shape[0]):
        weights = aipy.dsp.gen_window(data.shape[1],'blackman-harris')
        d = np.fft.ifft(data[t,:]*np.logical_not(flag[t,:]).astype(np.float)*weights)
        ker = np.fft.ifft(np.logical_not(flag[t,:])*weights)
        gain = aipy.img.beam_gain(ker)
        if not np.all(d == 0):
            _d,info = aipy.deconv.clean(d,ker,tol=1e-4)
            _d += info['res']/gain
        _d = np.fft.fftshift(np.ma.array(_d),axes=0)
        dtrans.append(_d)
    
    dtrans = np.array(dtrans)
    return {'waterfall':data, 'delayfall':dtrans, 'timeinfo':timeinfo,'delays':delays,'freqs':freqs}
"""
pol = 'I'

old = sorted(glob.glob('Stokes%s/zen.2455819.5*P'%pol)+glob.glob('Stokes%s/zen.2455819.6*P'%pol))
new = sorted(glob.glob('Stokes%s_newcal/zen.2455819.5*P'%pol)+glob.glob('Stokes%s_newcal/zen.2455819.6*P'%pol))

old = get_delayfall(old,'6_14',pol,blo=90,bhi=127)
new = get_delayfall(new,'6_14',pol,blo=90,bhi=127)

f, (ax1, ax2) = pylab.subplots(1, 2, sharey=True)
figold = ax1.imshow(np.abs(old['delayfall']),aspect='auto',interpolation='None',extent=(old['delays'][0],old['delays'][-1],old['delayfall'].shape[0],0))

fignew = ax2.imshow(np.abs(new['delayfall']),aspect='auto',interpolation='None',extent=(new['delays'][0],new['delays'][-1],new['delayfall'].shape[0],0))

cbar_ax = f.add_axes([0.92, 0.15, 0.01, 0.7])
cb = f.colorbar(figold, cax=cbar_ax)
pylab.show()

print np.nanmean(np.abs(old['delayfall'])),np.nanmean(np.abs(new['delayfall']))
"""
import IPython; IPython.embed()
