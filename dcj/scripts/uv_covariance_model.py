#! /usr/bin/env python

import aipy as a, numpy as n, pylab as p, capo, capo.frf_conv as fringe
import glob, optparse, sys, random
from capo.dcj import condition_goodness
from matplotlib.pyplot import *
def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = n.array(m, ndmin=2, dtype=n.complex128)
    X -= X.mean(axis=1)[(slice(None),n.newaxis)]
    N = X.shape[1]
    fact = float(N - 1) #normalization
    return (n.dot(X, X.T.conj()) / fact).squeeze()
def cov2(M,N):
    M = np.array(M,ndmin=2,dtype=np.complex128)
    M -= M.mean(axis=1)[(slice(None),np.newaxis)]
    N = np.array(N,ndmin=2,dtype=np.complex128)
    N -= N.mean(axis=1)[(slice(None),np.newaxis)]
    fact = np.sqrt(M.shape[1]-1)*np.sqrt(N.shape[1]-1)
    return (np.dot(M,N.T.conj())/fact).squeeze()
def diffdata(data1,data2):
    outdata = {}
    for bl in data:
        outdata[bl] = outdata.get(bl,{})
        for pol in data[bl]:
            outdata[bl][pol] = data1[bl][pol] - data2[bl][pol]
    return outdata

o = optparse.OptionParser()
a.scripting.add_standard_options(o, ant=True, pol=True, chan=True, cal=True,dec=True)
o.add_option('--plot',action='store_true',
    help='Plot plots')
opts,args = o.parse_args(sys.argv[1:])

print "split data into into even/odd based on path name"
evens = np.sort([filename for filename in args if filename.find('even')>-1])
odds = np.sort([filename for filename  in args if filename.find('odd')>-1])
assert(len(evens)==len(odds))
print "found ",len(evens), "data files"

print "loading evens..",
oddinfo, evendata, flag = capo.miriad.read_files(evens,opts.ant,opts.pol,decimate=opts.decimate,recast_as_array=True)
print "loading odds..",
eveninfo, odddata, flag = capo.miriad.read_files(odds,opts.ant,opts.pol,decimate=opts.decimate,recast_as_array=True)


print "found {n} times".format(n=len(oddinfo['times']))
print "comparing even/odd times to make sure we're congruent"
print "the largest lst discrpancy is ",np.max(np.abs(oddinfo['lsts'] - eveninfo['lsts']))*12/np.pi*3600,"s"

Cs = []
conds = []
goods = []
chans = None
for bl in evendata:
    for pol in evendata[bl]:
        if chans is None:
            chans = a.scripting.parse_chans(opts.chan,evendata[bl][pol].shape[1])
        #D = data[bl][pol][:,chans].T.astype(np.complex128)
        D_even = evendata[bl][pol][:,chans].T.astype(np.complex128)
        D_odd = odddata[bl][pol][:,chans].T.astype(np.complex128)

        #difference
        C_diff = cov(D_even-D_odd)
        S_diff = np.linalg.eigvals(C_diff)
        _C_diff = np.linalg.inv(C_diff)


        #individual
        C_even = cov(D_even)
        S_even,V_even = np.linalg.eig(C_even)
        _C_even = np.linalg.inv(C_even)
        #delay filter the eigenvectors
        #for now the horizon is hardcoded at 10 bins wide 
        dV_even = np.fft.ifft(V_even,axis=0)
        dV_even[:,5:-5] = 0 # mask out stuff outside the horizon
        V_even_filtered = np.fft.fft(dV_even,axis=0)
        S_even[4:] = 0 #kill all but the top 4 eigs. cut chosen by inspection...
        #cov with only top 4 eigs and those are delay filtered at the horizon
        C_even_robusto = np.einsum('ij,j,jk',V_even_filtered,S_even,V_even_filtered.T)
        #the projection matrix version (with the nonzero eigs replaced by 1s
        P_even = np.einsum('ij,j,jk',V_even_filtered,S_even!=0,V_even_filtered.T)

        #make a projected noise PNP
        PNP = np.einsum('ij,jk,kl',P_even,C_diff,P_even.T)
        #subtract the projected noise
        C_FG = C_even_robusto - PNP
        S_FG = np.linalg.eigvals(C_FG)

        _C_FG = np.linalg.inv(C_FG)
        #C_N should be the "ideal" noise 
        C_N = np.identity(C_FG.shape[0])*np.median(S_diff)
        # or the real noise 
        #C_N = C_diff
        C = C_FG + C_N
        _C = np.linalg.inv(C)
        S = np.linalg.eigvals(C)
        #weighted data
        D_even_weighted = D_even.T.dot(_C).T


        C_odd = cov(D_odd)
        S_odd = np.linalg.eigvals(C_odd)
        _C_odd = np.linalg.inv(C_odd)

        #cross
        C_cross = cov2(D_even,D_odd.conj())
        U_cross,S_cross,V_cross = n.linalg.svd(C_cross.conj())
        _C_cross = n.einsum('ij,j,jk', V_cross.T, 1./S_cross, U_cross.T)

        Cs.append(C)
        conds.append(n.log(n.abs(n.linalg.cond(C))/n.log(2)))
        goods.append(condition_goodness(C))
        print pol,bl,n.round(conds[-1],2),goods[-1]
        if opts.plot:
            #plot the dspecs
            figure()
            subplot(221)
            imshow(np.log(np.abs(n.fft.ifft(V_even,axis=0))));title('C_FG in delay space')
            subplot(222)
            plot(np.abs(np.fft.ifft(V_even,axis=0)[:,:6]));title('first 5 eigenvectors')
            xlabel('delay bins')
            subplot(212)
            plot(V_even_filtered[:,:6]);title('delay filtered eigenvectors')
            xlabel('channel')

            

            #plot the covariances
            mx,mn = np.log(np.abs(C_even)).max(),np.log(np.abs(C_even)).min()
            figure()
            subplot(241)
            imshow(np.log(np.abs(C_even)),vmin=mn,vmax=mx,interpolation='nearest')
            title('even C')
            subplot(242)
            imshow(np.log(np.abs(_C_even)),interpolation='nearest')
            title('C even inv')


            mxd,mnd = np.log(np.abs(C_diff)).max(),np.log(np.abs(C_diff)).min()
            subplot(243)
            imshow(np.log(np.abs(C_diff)),vmin=mnd,vmax=mxd,interpolation='nearest') 
            title('diff_C')
            subplot(244)
            imshow(np.log(np.abs(C_N)),vmin=mnd,vmax=mxd,interpolation='nearest')
            title('C_N')



            subplot(245)
            imshow(np.log(np.abs(C_FG)),vmin=mn,vmax=mx,interpolation='nearest')
            title('C_FG')
            subplot(246)
            imshow(np.log(np.abs(_C_FG)),interpolation='nearest')
            title('C_FG inv')

            subplot(247)
            imshow(np.log(np.abs(C)),vmin=mn,vmax=mx,interpolation='nearest')
            title('C = C_FG + C_N')
            subplot(248)
            imshow(n.log(n.abs(_C)),interpolation='nearest')
            title('C inv')



            tight_layout()

            #plot the eigenvectors
            figure()
            semilogy(S,label='C')
            semilogy(S_diff,label='diff')
            semilogy(S_even,label='even')
            #semilogy(S_odd,label='odd')
            semilogy(np.abs(S_FG),label='FG')
            semilogy(np.abs(S_odd - S_diff), label='odd - diff')
            legend(loc='best')
            title('eigenvalues')
            xlabel('eigenmode')
            tight_layout()

            #plot the data
            figure()
            subplot(311)
            imshow(np.abs(D_even),interpolation='nearest',aspect='auto');title('even data')
            subplot(312)
            imshow(np.abs(D_odd-D_even),interpolation='nearest',aspect='auto');title('even odd difference')
            subplot(313)
            imshow(np.abs(D_even_weighted),interpolation='nearest',aspect='auto');title('weighted data')
            xlabel('time')
            
            tight_layout()

            show()
            
            


