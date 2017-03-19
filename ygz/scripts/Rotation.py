import aipy as a, numpy as n, pylab as p, ephem as e
#Plot tracks of the entire array as the earth rotates
#aa=a.cal.get_aa('psa6622_v001',n.array([.15]))
aa = a.cal.get_aa('psa6622_v003',n.array([.15]))
nants = len(aa)
nants = 112
rad2deg=180/n.pi
src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15)
#src=a.fit.RadioSpecial("Sun")

#aa.set_jultime(2456240.2)
dt = 0.01
TIME = n.arange(2456249.3,2456249.5, dt)
U, V = {}, {}

for time in TIME:
    aa.set_jultime(time)
    src.compute(aa)
    for i in range(nants):
        for j in range(i+1,nants):
            if src.alt>0:
                u,v,w = aa.gen_uvw(i,j,src=src)
                u,v = u.flatten().flatten(), v.flatten().flatten()
                U[(i,j)] = U.get((i,j),[])+[u]
                V[(i,j)] = V.get((i,j),[])+[v]
                #p.plot(u,v,'.',ms=2,color = (((j*17)%127)/127., (32-i)*17%127/127., (i+j)/357., 1))
                #print u,v
        
            #p.plot(-u,-v,'ko')
    

#rs = 10**n.arange(1,2.5,rstep)
#rs = 2**(n.arange(3,8,1) +.5)
#for r in rs:
#    th = n.arange(0, 2*n.pi+.02, .01)
#    x,y = r*n.cos(th), r*n.sin(th)
#    p.plot(x,y,'r-')
import seaborn as sns
sns.set_context("paper")
sns.set(style="ticks", color_codes=True,font='DejaVu Serif', font_scale=2)
p.rc('axes', linewidth=2.5)
p.figure()
#import IPython; IPython.embed()
for k, u in U.iteritems(): 
    v = V[k]
    p.plot(u,v)
#p.xlim(-200,200)
#p.ylim(-200,200)
p.gcf().subplots_adjust(bottom=0.15, left=0.15)

p.grid()
p.xlabel('u')
p.ylabel('v')
p.show()
