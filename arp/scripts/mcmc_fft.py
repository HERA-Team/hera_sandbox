#! /usr/bin/env python
import numpy as n

def mcmc_fft(x, wgts, maxiter=1000, noise_lev=.0001, callback=None, maxhistory=500):
    ker = n.fft.fft(wgts)
    gain = n.sqrt(n.average(wgts**2))
    mdl = n.fft.fft(x) / gain
    res = x - wgts * n.fft.ifft(mdl)
    score = n.sqrt(n.sum(n.abs(res)**2) / n.sum(n.abs(wgts**2)))
    history = [(mdl,score)]
    if callback: callback(history)
    for i in xrange(maxiter):
        #_mdl = mdl + n.random.normal(scale=n.abs(n.fft.fft(res)*.1/gain).clip(noise_lev,n.Inf)) * n.exp(1j*n.random.uniform(0,2*n.pi,size=res.size).astype(n.complex))
        _mdl = mdl + n.fft.fft(res)*.1/gain 
        _mdl += n.random.normal(scale=noise_lev, size=_mdl.size) * n.exp(1j*n.random.uniform(0,2*n.pi,size=_mdl.size).astype(n.complex))
        _res = x - wgts * n.fft.ifft(_mdl)
        _score = n.sqrt(n.sum(n.abs(_res)**2) / n.sum(n.abs(wgts**2)))
        if n.random.uniform() < n.exp((score - _score) / noise_lev):
            mdl,res,score = _mdl,_res,_score
            history.append((mdl,score))
            if len(history) > maxhistory: history = history[-maxhistory:]
            if callback: callback(history)
        #else: print 'Rejecting score of', _score
    best = sum([h[0] for h in history]) / len(history)
    return best
   

if __name__ == '__main__':
    import pylab as p
    p.ion()
    spec = 10.**-n.arange(0,1,1./256, dtype=n.complex) * n.exp(1j*n.random.uniform(0,2*n.pi,size=256).astype(n.complex))
    print spec
    x = n.fft.ifft(spec)
    #wgts = n.random.uniform(0,1,size=256)
    #wgts = n.where(wgts < .25, 0., 1)
    wgts = n.ones_like(x); wgts[10:20] = 0
    _x = wgts * x
    line1, = p.plot(n.abs(spec), 'k')
    line2, = p.plot(n.zeros_like(x), 'r')
    def callback(history):
        mdl,score = history[-1]
        print 'Accepting jump with score', score
        line2.set_ydata(n.abs(mdl))
        p.draw()
    best = mcmc_fft(_x, wgts, callback=callback, noise_lev=1e-3, maxiter=100000, maxhistory=1000) 
    p.plot(n.abs(best), 'b')
    p.show()
