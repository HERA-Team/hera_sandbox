#! /usr/bin/env python
import aipy as ap, numpy as np, pylab as plt, capo
import sys

def cov(m):
    '''Because numpy.cov is stupid and casts as float.'''
    #return n.cov(m)
    X = np.array(m, ndmin=2, dtype=np.complex)
    X -= X.mean(axis=1)[(slice(None),np.newaxis)]
    N = X.shape[1]
    fact = float(N - 1)
    return (np.dot(X, X.T.conj()) / fact).squeeze()

aa = ap.cal.get_aa('psa6622_v003', np.array([.150]))
bls,conj = capo.red.group_redundant_bls(aa.ant_layout)
POL = 'xx'
#antstr = '49_64,65_66,72_73,80_81,88_89,96_97,104_105,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_49,66_67,73_74,81_82,89_90,97_98,105_106,3_25,1_58,35_61,42_63,2_33,15_21,8_45,19_41,47_67,74_75,82_83,90_91,98_99,106_107,25_48,1_4,18_35,37_42,6_33,15_16,8_11,19_29,47_68,75_76,83_84,91_92,99_100,107_108,24_48,4_17,5_18,37_40,6_52,16_62,11_36,28_29,68_69,76_77,84_85,92_93,100_101,108_109,24_55,13_17,5_32,14_40,7_52,44_62,36_60,28_34,69_70,77_78,85_86,93_94,101_102,109_110,27_55,13_56,30_32,14_54,7_12,0_44,39_60,34_51,70_71,78_79,86_87,94_95,102_103,110_111,27_57,56_59,23_30,50_54,12_38,0_26,39_46'
#antstr = '1_4,49_64,65_66,72_73,80_81,88_89,96_97,104_105,3_10'
antstr = '49_64,65_66,72_73,80_81,88_89,96_97,104_105,3_10,9_58,22_61,20_63,2_43,21_53,31_45,41_49,66_67,73_74,81_82,89_90,97_98,105_106,3_25,1_58,35_61,42_63,2_33,15_21,8_45,19_41,47_67,74_75,82_83,90_91,98_99,106_107,25_48,1_4'

t,d,f = capo.arp.get_dict_of_uv_data(sys.argv[1:], antstr=antstr, polstr=POL, verbose=True)

for bl in d:
    blm = ap.miriad.ij2bl(*bl)
    for pol in d[bl]:
        if conj[blm]: d[bl][pol] = d[bl][pol].conj()

C = {}
for i,bl in enumerate(d.keys()):
    if not C.has_key(bl): C[bl] = {}
    for pol in d[bl]:
        di = d[bl][pol]
        wi = np.where(f[bl][pol], 0., 1.)
        #Ci = np.dot(di.T, di.conj())
        Ci = cov(di.T).T
        Wi = np.dot(wi.T, wi.conj())
        C[bl][pol] = np.where(Wi > 0, Ci / Wi, 0)

davg = np.zeros_like(di)
times = range(di.shape[0])

if True:
    plt.ion()
    plt.subplot(121)
    #plt1 = capo.arp.waterfall(davg, mx=3, drng=3)
    plt1 = capo.arp.waterfall(davg, mx=0, drng=4)
    plt.subplot(122)
    plt2 = capo.arp.waterfall(d[(1,4)][POL], mx=0, drng=4)
    plt.draw()

for t in times:
    print t, times[-1]
    _Cdt = 0.
    _Csum = 0.
    for bl in d:
        for pol in d[bl]:
            Cbl = C[bl][pol]
            #I = 1e-8 * np.identity(Cbl.shape[0])
            I = 0
            dt, wt = d[bl][pol][t:t+1], np.where(f[bl][pol][t:t+1], 0., 1.)
            Wt = np.dot(wt.T, wt.conj())
            #_Ct = np.linalg.pinv(Cbl * Wt)
            #_Ct = np.linalg.pinv((Cbl+I) * Wt)
            if False:
                U,S,V = np.linalg.svd((Cbl+I) * Wt)
                cutoff = 0.1 * np.median(S)
                print S[0], np.median(S)
                _Ct = np.dot(V.T.conj(), np.dot(np.diag(1/(S+cutoff)), U.T.conj()))
            else:
                _Ct = np.identity(Cbl.shape[0])
            _Ct *= Wt
            #S = np.ones_like(S)
            #_Ct = np.dot(V.T.conj(), np.dot(np.diag(1/(S+2e-2)), U.T.conj()))
            #_Ct = np.dot(V.T.conj(), np.dot(np.diag(1/(S)), U.T.conj()))
            #_Ct = np.identity(S.shape[0])
            if False:
                plt.subplot(221); capo.arp.waterfall(Cbl, drng=5)
                plt.subplot(222); capo.arp.waterfall(Cbl*Wt, drng=5)
                plt.subplot(223); capo.arp.waterfall(np.linalg.pinv(Cbl), drng=5)
                plt.subplot(224); capo.arp.waterfall(_Ct, drng=5)
                plt.show()
            _Cdt += np.dot(_Ct, dt.T).T
            _Csum += _Ct
    #Cavg = np.linalg.inv(_Csum)
    try: Cavg = np.linalg.pinv(_Csum) * Wt
    except(np.linalg.LinAlgError): continue
    #Cavg = np.identity(Cavg.shape[0])
    davg[t] = np.dot(Cavg, (_Cdt).T).T
    plt1.set_data(np.log10(np.abs(davg)))
    plt.draw()
            
plt.show()
import IPython; IPython.embed()
Cov1 = dot(d1.T, d1.conj())
#Cov2 = dot(d1.T, d1.conj(), med=True)

#p.subplot(121); C.arp.waterfall(Cov1, mx=0, drng=3); p.colorbar()
#p.subplot(122); C.arp.waterfall(Cov2, mx=0, drng=3); p.colorbar(); p.show()

#U2,S2,V2 = n.linalg.svd(Cov2)
#iCov2 = n.dot(V2.T.conj(), n.dot(n.diag(1/(S2+.01)), U2.T.conj()))

ad1 = n.abs(d1)
d1C1 = n.where(ad1 > 0, n.dot(iCov1,d1.T).T, 0)
#d1C2 = n.where(ad1 > 0, n.dot(iCov2,d1.T).T, 0)

p.subplot(131); C.arp.waterfall(d1, drng=3); p.colorbar()
p.subplot(132); C.arp.waterfall(d1C1, drng=3); p.colorbar()
#p.subplot(133); C.arp.waterfall(d1C2, drng=3); p.colorbar()
p.subplot(133); C.arp.waterfall(chisq, drng=5); p.colorbar()
p.show()

import IPython; IPython.embed()

