#! /usr/bin/env python
import capo.omni as omni, numpy as np, pylab as p, aipy as a, capo.zsa as zsa
import optparse, sys

o=optparse.OptionParser()
a.scripting.add_standard_options(o,cal=True)
opts,files = o.parse_args(sys.argv[1:]) #files are arguments.

def gen_delay_limits(aa, freqs, max_bl_add=0.):
    dly_per_bin = 1/((freqs[1] - freqs[0])*len(freqs))
    filters = {}
    for i in range(len(aa.ants)):
        for j in range(len(aa.ants)):
            if j<i: continue
            max_bl = aa.get_baseline(i,j)
            max_bl = max_bl_add + np.sqrt(np.dot(max_bl, max_bl))
            uthresh, lthresh = max_bl/dly_per_bin + 1.5, -max_bl/dly_per_bin - 0.5
            uthresh, lthresh = int(np.ceil(uthresh)), int(np.floor(lthresh))
            filters[(i,j)] = (uthresh, lthresh)

    return filters

#def delay_filter(d,w,uthresh,lthresh):
    

freqs = np.arange(.1,.2,.1/1024)
aa = a.cal.get_aa('hsa7458_v000', freqs)
filters = gen_delay_limits(aa, freqs) 
window = a.dsp.gen_window(len(freqs), window="blackman-harris")




#get data
#m,g,v,x = omni.from_npz(files, verbose=True)
m,g,v,x,f = zsa.flag_by_chisq(files)
pol='xx'
bandpass = np.ones_like(v['xx'][(9,88)][0])
#loop over pol, baseline, and time.
niters=0
while niters!=10:
    betas = np.zeros_like(v['xx'][(9,88)][0])
    count = 0
    for bl in v[pol].keys():
        for vis,flags in zip(v[pol][bl],f):
            w = flags
            d = vis*w
            _d = np.fft.ifft(d*window)
            _w = np.fft.ifft(w*window)
            ut,lt = filters[bl]
            area = np.ones(_d.size, dtype=np.int)
            _d_cl, info = a.deconv.clean(_d, _w, tol=1e-5, area=area, stop_if_div=False, maxiter=100)
            d_cl = np.fft.fft(_d_cl)
            betas += d/d_cl
            count+=1
             
    bandpass*=betas/count                
    niters+=1
    p.plot(freqs, bandpass)        
    p.title(str(niters))
    p.show() 
