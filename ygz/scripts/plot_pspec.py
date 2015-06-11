__author__ = 'yunfanzhang'
import pylab as p, numpy as n
def P_v_Eta(k,P):
    fig = p.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel('kz')
    ax.set_yscale('log')
    p.plot(k,P)
    p.show()