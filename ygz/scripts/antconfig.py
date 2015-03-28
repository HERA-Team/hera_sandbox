#plots configuration of antenna array
import aipy as a, numpy as n, pylab as p, ephem as e

aa=a.cal.get_aa('psa6622_v001',n.array([.15]))
#aa=a.cal.get_aa('paper128',n.array([.15]))
nants=128
rad2deg=180/n.pi
src= a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
#src=a.fit.RadioSpecial("Sun")
U = []
V = []
I = []
J = []
aa.set_jultime(2456240.3)
src.compute(aa)
for i in range(nants):
    #print i
    u,v,w = aa.gen_uvw(i,0,src=src)
    u,v = u.flatten(), v.flatten()
    U.append(u)
    V.append(v)
    I.append(i)

        #p.plot(time, aa.sidereal_time(),'k.')

fig = p.figure()
ax = fig.add_subplot(111)
p.plot(U,V,'.',ms=5)
for u,v,i in zip(U, V,I):
    #print u,v,i# <--
    ax.annotate('%s' %i, xy=(u,v), textcoords='data') # <--
p.grid()
#p.xlim(-200,200)
#p.ylim(-200,200)
p.xlabel('u',size=14)
p.ylabel('v',size=14)
p.show()
