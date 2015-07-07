__author__ = 'yunfanzhang'
import pylab as p, numpy as n
def P_v_Eta(k,P):
    fig = p.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('kz')
    ax.set_ylabel(r'$P(k) mK^{2} (h^{-1} Mpc)^{3}$')
    ax.set_yscale('log')
    p.plot(k,P,'bo')
    p.show()
