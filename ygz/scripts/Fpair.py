__author__ = 'yunfanzhang'
import aipy as a, numpy as n

def Fpair_coarse(aa,nants, src,TIME,dt,dist):


    d={}
    time=TIME[0]+10*dt
    aa.set_jultime(time)
    src.compute(aa)
    dist=10*dist
    for i in range(nants):
        for j in range(i+1,nants):
            u0,v0,w0 = aa.gen_uvw(i,j,src=src)
            u0r,v0r = n.around(u0[0,0]/dist,decimals=0)*dist, n.around(v0[0,0]/dist,decimals=0)*dist
            d[u0r,v0r]=d.get((u0r,v0r),[])+[([i,j],time,[u0[0,0],v0[0,0]])]

    repbl=[]
    for key in d.keys():
        repbl.append(d[key][0][0])

    #print repbl
   # for i,j in repbl:
   #     print i,j, aa.gen_uvw(i,j,src=src)[0][0][0],aa.gen_uvw(i,j,src=src)[1][0][0]

    d={}

    print 'd=', d

    dist=dist/10
    for time in TIME:
        aa.set_jultime(time)
        src.compute(aa)
        print time, n.round(time,decimals=2)
        for i,j in repbl:
            u1,v1,w1 = aa.gen_uvw(i,j,src=src)
            u1r,v1r = n.around(u1.flatten()[0]/dist,decimals=0)*dist, n.around(v1.flatten()[0]/dist,decimals=0)*dist
            d[u1r,v1r]=d.get((u1r,v1r),[])+[([i,j],n.round(time,decimals=3),[u1[0,0],v1[0,0]])]
            #if n.sqrt(u**2+v**2)<distlim and delt>1.5*dt:
            #    antdict[key].append([i,j,time])
            #    #print key
            #else:
            #    antdict[u1,v1]=[[i,j,time]]


    for key in d.keys():
        if len(d[key])<2:
            del d[key]

    return d