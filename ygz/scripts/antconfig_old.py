#plots configuration of antenna array
import aipy as a, numpy as n, pylab as p, ephem as e

aa = a.cal.get_aa('psa6240_v003',n.array([.15]))
#aa=a.cal.get_aa('paper128',n.array([.15]))
nants = 64
rad2deg = 180/n.pi
ltsec = 299792458.  #meters
U,V,I,J,X,Y = [],[],[],[],[],[]

aa.set_jultime(2456240.365)
src = a.fit.RadioFixedBody(aa.long, aa.lat, janskies=0., mfreq=.15, name='test')
src.compute(aa)
for i in range(nants):
    #print i
    u,v,w = aa.gen_uvw(i,0,src=src)
    u,v = u.flatten(), v.flatten()
    U.append(u)
    V.append(v)
    I.append(i)

#pos0 = aa.ants[0].pos
#i = 0
#for ant in aa.ants:
#    Y.append(ant.pos[0]-pos0[0])
#    X.append(ant.pos[1]-pos0[1])
#    if i<112: J.append(aa.ant_layout.flatten()[i])
#    i = i+1
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
