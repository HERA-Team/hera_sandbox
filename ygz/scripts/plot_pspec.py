__author__ = 'yunfanzhang'
import pylab as p, numpy as n
def P_v_Eta(ax,k,P):
    ax.set_xlabel('kz')
    ax.set_ylabel(r'$P(k) K^{2} (h^{-1} Mpc)^{3}$')
    #ax.set_yscale('log')
    p.plot(k,P,'bo')

def P_v_Eta_log(ax,k,P):
    ax.set_xlabel('kz')
    ax.set_ylabel(r'$P(k) K^{2} (h^{-1} Mpc)^{3}$')
    ax.set_yscale('log')
    p.plot(k,P,'bo')

def Del_v_Eta(ax,k,P,sgn=False):
    ax.set_xlabel('kz')
    ax.set_ylabel(r'$\Delta^{2}(k) K^{2}$')
    k, P = n.abs(n.array(k)), n.array(P)
    #Del = k*k*k/P/2/(n.pi**2)
    Del = k*k*k*P/2/(n.pi**2)
    p.plot(k,Del,'bo')
