import aipy as a, numpy as n, pylab as p, ephem as e
#Plot tracks of the entire array as the earth rotates
#aa=a.cal.get_aa('psa6622_v001',n.array([.15]))

def plot_c(src, TIME, ann=False, C=None,i=0,J=[26,23], ia=0):
    if C==None: C=['b','g','y','c','m']
    for j in J:
        c=C[j%len(C)]
        U,V=[],[]
        for time in TIME:
            aa.set_jultime(time)
            src.compute(aa)
            if src.alt>0:
                u,v,w = aa.gen_uvw(i,j,src=src)
                u,v = u.flatten(), v.flatten()
                U.append(u); V.append(v)
        p.plot(U,V,label='0_'+str(j),color=c,alpha=0.3)
        if ann == 'spread':
            #annotate times
            step = len(U)/10
            UA,VA,TA = U[1::step],V[1::step],TIME[1::step]-2456249
            ind=0
            for xy in zip(UA,VA):
                #ax.annotate('%s,%s' % xyt[2], xy=xyt[:1], textcoords='offset points') # <--
                ax.annotate('%.3f' % TA[ind], xy=xy, textcoords='data') # <--
                ind=ind+1
            p.scatter(UA,VA)
        if ann == 'T':
            UA,VA,TA = U[ia],V[ia],TIME[ia]-2456249
            p.scatter(UA,VA)
            #ax.annotate('%.3f' % TA, xy=(U,), textcoords='data')
    #p.legend()
    return

aa = a.cal.get_aa('psa6240_v003',n.array([.15]))
nants = len(aa)
rad2deg=180/n.pi

#src=a.fit.RadioSpecial("Sun")

fig = p.figure()
ax = fig.add_subplot(111)
dt = 0.001
#TIME = n.arange(2456249.20,2456249.35, dt)
TIME = n.arange(2456249.0,2456250.0, dt)
dd = 0.05
DEC = aa.lat+n.arange(-1,0.41,dd)
rad2deg = 180./n.pi
deg2rad = n.pi/180

####################################
for dec in DEC:
   src = a.fit.RadioFixedBody(0, dec, janskies=0., mfreq=.15)
   C = None
   if abs(dec-aa.lat)<dd/2:
       ann='spread',
       C = ['r']
   elif abs(dec-60*deg2rad)<dd:
       ann='spread',
       C = ['y']
   else: ann=False
   plot_c(src,TIME,ann=ann,C=C)

#############################################
# for dec in DEC:
#     src = a.fit.RadioFixedBody(0, dec, janskies=0., mfreq=.15)
#     plot_c(src,TIME,ann='T',ia=200)
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
#import IPython; IPython.embed()
p.show()
