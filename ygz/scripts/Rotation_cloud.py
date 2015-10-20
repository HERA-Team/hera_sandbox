import aipy as a, numpy as n, pylab as p, ephem as e
#Plot tracks of the entire array as the earth rotates
#aa=a.cal.get_aa('psa6622_v001',n.array([.15]))

def plot_c(src, TIME, i=0,J=[26,38]):
    for j in J:
        U,V=[],[]
        for time in TIME:
            aa.set_jultime(time)
            src.compute(aa)
            if src.alt>0:
                u,v,w = aa.gen_uvw(i,j,src=src)
                u,v = u.flatten(), v.flatten()
                U.append(u); V.append(v)
        p.plot(U,V,label='0_'+str(j),color=str(float(j)/64))
    #for xy in zip(U, V):                                                # <--
    #    ax.annotate('(%s, %s)' % xy, xy=xy, textcoords='offset points') # <--
    #p.legend()
    return

aa = a.cal.get_aa('psa6240_v003',n.array([.15]))
nants = len(aa)
rad2deg=180/n.pi

#src=a.fit.RadioSpecial("Sun")

fig = p.figure()
ax = fig.add_subplot(111)
dt = 0.001
TIME = n.arange(2456249.25,2456249.35, dt)
DEC = aa.lat+n.arange(-0.3,0.3,0.1)
for dec in DEC:
    src = a.fit.RadioFixedBody(0, dec, janskies=0., mfreq=.15)
    plot_c(src,TIME)


            #print u,v

            #p.plot(-u,-v,'ko')


#rs = 10**n.arange(1,2.5,rstep)
#rs = 2**(n.arange(3,8,1) +.5)
#for r in rs:
#    th = n.arange(0, 2*n.pi+.02, .01)
#    x,y = r*n.cos(th), r*n.sin(th)
#    p.plot(x,y,'r-')

p.grid()
#p.xlim(-200,200)
#p.ylim(-200,200)
p.xlabel('u',size=14)
p.ylabel('v',size=14)
p.show()
