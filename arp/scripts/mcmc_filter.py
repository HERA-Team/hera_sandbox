#! /usr/bin/env python
import aipy as a, numpy as n, pylab as P
import sys

maxiter = 1000
P.ion()

uv = a.miriad.UV(sys.argv[-1])
uv.select('antennae', 2, 9)
for p,d,f in uv.all(raw=True):
    mdl, score = 0., n.Inf
    P.subplot(121)
    line0, = P.semilogy(n.abs(d),'k')
    line1, = P.semilogy(n.abs(d),'b')
    line2, = P.semilogy(n.abs(d),'r')
    P.subplot(122)
    line0_, = P.semilogy(n.abs(n.fft.fft(d)),'k')
    line1_, = P.semilogy(n.abs(n.fft.fft(d)),'b')
    line2_, = P.semilogy(n.abs(n.fft.fft(d)),'r')
    P.ylim(1e-2,1e3)
    for iter in xrange(maxiter):
        ns = n.sqrt(n.median(n.abs(d-mdl)**2))
        ns_prob = 1 - n.exp(-n.abs(d-mdl)**2 / (3*ns)**2)
        _v = n.where(n.random.uniform(size=ns_prob.size) < ns_prob, 0., 1)
        _d = n.where(_v, d-mdl, 0)
        #_v_ = n.fft.fft(_v)
        _d_ = n.fft.fft(_d)
        ns_ = n.sqrt(n.median(n.abs(_d_)**2))
        ns_prob_ = 1 - n.exp(-n.abs(_d_)**2 / (2*ns_)**2)
        _v_ = n.where(n.random.uniform(size=ns_prob_.size) < ns_prob_, 1., 0)
        if iter > 0: _mdl_ = n.where(_v_, .1 * _d_ + n.fft.fft(mdl), 0)
        else: _mdl_ = n.where(_v_, .1 * _d_, 0)
        #_mdl_ = n.where(_v_, .1 * _d_, 0)
        #_mdl_ = _d_.copy()
        #_mdl_[10:-10] = 0
        #_mdl = mdl + n.fft.ifft(_mdl_)
        _mdl = n.fft.ifft(_mdl_)
        _score = n.sqrt(n.sum(n.abs((_d - _mdl)*_v)**2) / n.sum(_v**2))
        #flagging_penalty = n.log(n.sum(1-_v)).clip(1,n.Inf)
        #flagging_penalty = n.sqrt(n.sum(1-_v)).clip(1,n.Inf)
        flagging_penalty = (1-_v).sum() + _v_.sum()
        _score *= flagging_penalty
        if n.random.uniform() < n.exp((score - _score) / .1):
            mdl,score = _mdl,_score
            line1.set_ydata(n.abs(_mdl))
            line2.set_ydata(n.abs(_d - _mdl)*_v)
            line1_.set_ydata(n.abs(n.fft.fft(_mdl) + _mdl_))
            line2_.set_ydata(n.abs(_d_))
            print iter, score, ns, _v.sum(), _v_.sum(), flagging_penalty
            P.draw()
    P.show()
    sys.exit(0)
