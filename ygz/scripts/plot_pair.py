__author__ = 'yunfanzhang'
import pylab as p
#plot all the closest approach points
def plot_closapp(clos_app,corr_tol,figname):
    fig = p.figure()
    ax = fig.add_subplot(111)
    U,V = [],[]
    for key in clos_app.keys():
        for item in clos_app[key]:
            U.append(item[3][0])
            V.append(item[3][1])
    p.plot(U,V,'.',ms=5)
    p.grid()
    p.xlabel('u',size=14)
    p.ylabel('v',size=14)
    p.savefig(figname)

#plot the pairings near one of the crossings
def plot_pair_xampl(pair_xampl,figname='PairSample.png'):
    fig2 = p.figure()
    ax2 = fig2.add_subplot(111)
    for i in n.arange(len(pair_xampl)):
        p.plot([pair_xampl[i][0][0], pair_xampl[i][1][0]], [pair_xampl[i][0][1], pair_xampl[i][1][1]],color = 'brown', marker = 'o')
    p.grid()
    p.xlabel('u',size=14)
    p.ylabel('v',size=14)
    p.savefig(figname)
