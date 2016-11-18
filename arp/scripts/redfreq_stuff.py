#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p, capo as C
import sys

seps = [
    '<0,103> xx', #  1
    '<1,4> xx',   #  2
    '<0,101> xx', #  3
    '<0,62> xx',  #  4
    '<0,100> xx', #  5
    '<1,13> xx',  #  6
    '<1,70> xx',  #  7
    '<1,56> xx',  #  8
    '<1,71> xx',  #  9
    '<1,59> xx',  # 10
    '<0,97> xx',  # 11
    '<12,43> xx', # 12
    '<9,71> xx',  # 13
    '<9,59> xx',  # 14
    '<57,64> xx', # 15
]

conj = [
    '<0,103> xx', #  1
    '<1,4> xx',   #  2
    '<0,101> xx', #  3
    '<0,62> xx',  #  4
    '<0,100> xx', #  5
    '<0,97> xx',
    '<12,43> xx',
    '<57,64> xx']

def get_data(filelist, seps, conj):
    d,w = {}, {}
    for filename in filelist:
        print 'Reading', filename
        npz = n.load(filename)
        for sep in seps: d[sep] = d.get(sep,[]) + [npz[sep]]
    for sep in seps:
        d[sep] = n.concatenate(d[sep])
        w[sep] = n.where(d[sep] == 0, 1., 0) # XXX add in noise estimate here
    for sep in conj:
        try: d[sep] = d[sep].conj()
        except(KeyError): pass
    return d,w

def gen_mset_dict(d, seps, fqi, BL_MAX=15, umin=0):
    chi = n.arange(fqi.size)
    sdf = fqi[1]-fqi[0]
    M = {}
    for i,si in enumerate(seps):
        for j,sj in enumerate(seps):
            if j <= i: continue
            Cij = n.dot(d[si].T, d[sj].conj())
            bi,bj = float(i+1),float(j+1) # proxy for baseline length
            if bi < umin: continue
            #print bi,si, bj, sj
            fqj = fqi * bi/bj
            chj = n.around((fqj - fqi[0]) / sdf).astype(n.int)
            inds = n.where(n.logical_and(chj >= 0, chj < fqi.size))[0]
            if len(inds) < 4: continue # XXX why 4?
            iis,jjs = chi[inds],chj[inds]
            uui,uuj = n.around(fqi[iis] * bi / (BL_MAX*sdf)), n.around(fqi[jjs] * bj / (BL_MAX*sdf))
            rrs = Cij[iis,jjs]
            for ci,cj,ui,uj,rij in zip(iis,jjs,uui,uuj,rrs):
                k = (ci,cj,int(ui))
                if ui != uj or n.isnan(rij): continue
                #if wij == 0: continue
                M[k] = M.get(k, []) + [rij]
    return M

def logcal(mset, iAtA_At=None, iBtB_Bt=None):
    chs, us = {}, {}
    for ci,cj,u in mset.keys():
        chs[ci] = chs[cj] = us[u] = None
    chs, us = chs.keys(), us.keys()
    chs.sort(); us.sort()
    ch2ind, u2ind = {}, {}
    for i,ch in enumerate(chs): ch2ind[ch] = i
    for i,u in enumerate(us): u2ind[u] = i
    if iAtA_At is None:
        A = n.zeros((len(mset), len(chs)+len(us)), dtype=n.float32)
        B = n.zeros((len(mset), len(chs)+len(us)), dtype=n.float32)
        for i,(ci,cj,u) in enumerate(mset.keys()):
            #d /= w.clip(1,n.inf)
            A[i,ch2ind[ci]] = 1
            A[i,ch2ind[cj]] = 1
            A[i,len(chs)+u2ind[u]] = 1
            B[i,ch2ind[ci]] = 1
            B[i,ch2ind[cj]] = -1
            B[i,len(chs)+u2ind[u]] = 1
        AtA = n.dot(A.T,A)
        BtB = n.dot(B.T,B)
        iAtA = n.linalg.pinv(AtA, rcond=1e-10)
        iBtB = n.linalg.pinv(BtB, rcond=1e-10)
        iAtA_At = n.dot(iAtA, A.T)
        iBtB_Bt = n.dot(iAtA, A.T)
    ma = n.zeros((len(mset),), dtype=n.float32)
    mb = n.zeros((len(mset),), dtype=n.float32)
    for i,(ci,cj,u) in enumerate(mset.keys()):
        d = mset[(ci,cj,u)][0]
        ma[i] = n.log10(n.abs(d).clip(1e-4,n.Inf))
        mb[i] = n.angle(d)
    a = n.dot(iAtA_At, ma)
    b = n.dot(iBtB_Bt, mb)
    bp_sol, u_sol = {}, {}
    for ch in chs: bp_sol[ch] = 10**a[ch2ind[ch]] * n.exp(1j*b[ch2ind[ch]])
    for u in us: u_sol[u] = 10**a[len(chs) + u2ind[u]] * n.exp(1j*b[len(chs) + u2ind[u]])
    return bp_sol, u_sol, iAtA_At, iBtB_Bt

def lincal(mset, bp_in, u_in, rcond=1e-5):
    chs, us = {}, {}
    for ci,cj,u in mset.keys():
        chs[ci] = chs[cj] = us[u] = None
    chs, us = chs.keys(), us.keys()
    chs.sort(); us.sort()
    ch2ind, u2ind = {}, {}
    for i,ch in enumerate(chs): ch2ind[ch] = i
    for i,u in enumerate(us): u2ind[u] = i
    A = n.zeros((len(mset), 2*len(chs)+2*len(us)), dtype=n.complex64)
    ma = n.zeros((len(mset),), dtype=n.complex64)
    for i,(ci,cj,u) in enumerate(mset.keys()):
        d = mset[(ci,cj,u)][0]
        #if w > 0: d,w = d/w,1.
        BBu = bp_in[ci] * bp_in[cj].conj() * u_in[u]
        A[i,ch2ind[ci]] = BBu
        A[i,ch2ind[cj]] = BBu
        A[i,len(chs)+ch2ind[ci]] = 1j*BBu
        A[i,len(chs)+ch2ind[cj]] = -1j*BBu
        A[i,2*len(chs)+u2ind[u]] = BBu
        A[i,2*len(chs)+len(us)+u2ind[u]] = 1j*BBu
        ma[i] = d - BBu
    #AtA = n.dot(A.T,A)
    AtA = n.dot(A.T.conj(),A)
    iAtA = n.linalg.pinv(AtA, rcond=rcond)
    iAtA_At = n.dot(iAtA, A.T.conj())
    #a = n.dot(iAtA, n.dot(A.T, ma))
    a = n.dot(iAtA_At, ma)
    bp_sol, u_sol = {}, {}
    for ch in chs: bp_sol[ch] = bp_in[ch] * n.exp(a[ch2ind[ch]]+1j*a[len(chs)+ch2ind[ch]])
    for u in us: u_sol[u] = u_in[u] * n.exp(a[2*len(chs) + u2ind[u]] + 1j*a[2*len(chs) + len(us) + u2ind[u]])
    return bp_sol, u_sol

import glob, pylab as plt
#days = ['979','980','981','982']
#for day in days:
#for f in sys.argv[1:]:
if True:
    #print f
    #files = [f]
    #files = glob.glob('*%s.*npz' % day)
    files = sys.argv[1:]
    d,w = get_data(files, seps, conj)
    iAtA_At, iBtB_Bt = None, None
    for i in xrange(d.values()[0].shape[0]):
        print i
        if i > 3: break
        _d = {}
        for k in d: _d[k] = d[k][i:i+1]
        fqi = n.linspace(.1,.2,_d.values()[0].shape[1])
        mset = gen_mset_dict(_d, seps, fqi)
        if i == 0: bp,us,iAtA_At,iBtB_Bt = logcal(mset, iAtA_At=iAtA_At, iBtB_Bt=iBtB_Bt)
        bp,us = lincal(mset, bp, us)
        chs = bp.keys(); chs.sort()
        bp_sol = n.array([bp[ch] for ch in chs])

        plt.subplot(211)
        plt.plot(n.abs(bp_sol))
        plt.subplot(212)
        plt.plot(n.angle(bp_sol))
plt.show()

import IPython; IPython.embed()
