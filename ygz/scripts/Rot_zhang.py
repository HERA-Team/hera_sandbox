import aipy as a, numpy as n, pylab as p, ephem as e

aa=a.cal.get_aa('psa6622_v001',n.array([.15]))
#aa=a.cal.get_aa('paper128',n.array([.15]))
nants=32
rad2deg=180/n.pi
src= a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
#src=a.fit.RadioSpecial("Sun")

p.figure()
#aa.set_jultime(2456240.2)
for i in range(nants):
  for j in range(i+1,nants):
     # for time in n.arange(2456240.2,2456241,.01):
          #  aa.set_jultime(time)
            #print aa.get_jultime()
            src.compute(aa)
            if src.alt>0:
                u,v,w = aa.gen_uvw(i,j,src=src)
                u,v = u.flatten(), v.flatten()
                #p.plot(time, aa.sidereal_time(),'k.')
                p.plot(u,v,'.',ms=2,color = (((j*17)%33)/33., (32-i)*17%33/33., (i+j)/65., 1)) 
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
