__author__ = 'yunfanzhang'
import pylab as p, numpy as n
#plot all the closest approach points
def plot_closapp(clos_app,corr_tol,figname):
    fig = p.figure()
    ax = fig.add_subplot(111)
    U,V = [],[]
    for key in clos_app.keys():
        U.append(clos_app[key][3][0])
        V.append(clos_app[key][3][1])
    p.plot(U,V,'.',ms=2)
    p.grid()
    p.xlabel('u',size=14)
    p.ylabel('v',size=14)
    p.savefig(figname)

#plot the pairings near one of the crossings
def plot_pair_xampl(pair_xampl,figname='PairSample.png'):
    fig2 = p.figure()
    ax2 = fig2.add_subplot(111)
    N = len(pair_xampl)
    for i in n.arange(len(pair_xampl)):
        p.plot([pair_xampl[i][0][0], pair_xampl[i][1][0]], [pair_xampl[i][0][1], pair_xampl[i][1][1]], marker = 'o', color = (i/1.0001/N, 0.5, 0.5))
    p.grid()
    p.xlabel('u',size=14)
    p.ylabel('v',size=14)
    p.savefig(figname)

