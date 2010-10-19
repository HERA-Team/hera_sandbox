#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse

o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
o.add_option('-o','--origin',dest='origin',type='str',
    help="""Origin of coordinate system.  Necessary for outputing relative
    azimuth and distance""")
o.add_option('-z',dest='z',default=10,type='float',
    help="""Plot height as circles multiplied by this factor""")
o.add_option('-p',dest='p',action='store_true',
    help="""Plot phase as circles multiplied by this facter (override height).""")

opts, args = o.parse_args(sys.argv[1:])

if not opts.origin is None: origin = n.array(map(float,opts.origin.split(',')))
else: origin = n.array([-1,-1,0])
print "origin = ",origin
th = n.arange(0, 2*n.pi, .01)
fmts = (('k.','k-'), ('r.','r-.'))

aas = [a.cal.get_aa(opts.cal, .1, .1, 1)]
if len(args)>0: aas.append(a.cal.get_aa(args[0], .1, .1, 1))
print aas
all_antpos = n.zeros((len(aas[0].ants),3,len(aas)) )
if len(aas)>1:
    for cnt, aa in enumerate(aas):
        antpos = [aa.get_baseline(0,i,src='z') for i in range(len(aa.ants))]
        antpos = n.array(antpos) * a.const.len_ns / 100.
        x,y,z = antpos[:,0], antpos[:,1], antpos[:,2]
        all_antpos[:,:,cnt]=antpos
        for l in antpos: print "[",l[0],",",l[1],",",l[2],"]"
        def fitfunc(vec):
            r,cx, cy = vec
            return ((n.sqrt((x-cx)**2 + (y-cy)**2) - r)**2).sum()
        if cnt==0: r,cx,cy = a.optimize.fmin(fitfunc, n.array([200., 0, 0]))
        print r, cx, cy
        circ_x, circ_y = r*n.cos(th), r*n.sin(th)
        if cnt==1: 
           x_c =all_antpos[:,0,0] - cx
           y_c =all_antpos[:,1,0] - cy
           x_s =all_antpos[:,0,1] - cx
           y_s =all_antpos[:,1,1] - cy

           dr_a = n.zeros(len(x_c))
           dt_a = n.zeros(len(x_c))
           for i in range(len(x_c)):
               dr_a[i]=m.sqrt(x_c[i]**2+y_c[i]**2) - r
               dt_a[i]=(m.atan(y_c[i]/x_c[i]) - m.atan(y_s[i]/x_s[i]))*r
           #p.plot(circ_x,circ_y,fmt2)
        fmt1,fmt2 = fmts[cnt % len(fmts)]
        p.plot(x,y, fmt1)
        for ant,(xa,ya,za) in enumerate(zip(x,y,z)):
            hx,hy = 10*za*n.cos(th)+xa, 10*za*n.sin(th)+ya
            if za > 0: fmt3 = '#eeeeee'
            else: fmt3 = '#a0a0a0'
            p.fill(hx,hy, fmt3,alpha=0.5)
        
    p.grid()
    delta_coords = n.diff(all_antpos,axis=2)*100 #deltas in cm
    #print x_c

    for d in delta_coords:
        print d[0][0],',',d[1][0],',',d[2][0]
    #p.ylim((-25,25))
    #p.xlim(-150,150)
    #p.ylim(-150,150)
    a = p.gca()
    a.set_aspect('equal')
    #print p.axes()
#    p.axes([0.4,0.4,0.2,0.2])
    p.figure()
    p.plot(delta_coords[:,0],delta_coords[:,1],'.')
    p.xlim(-25,25)
    p.xlabel("delta x [cm]")
    p.ylim(-25,25)
    p.ylabel("delta y [cm]")
    p.figure()
    p.subplot(131)
    p.hist(delta_coords[:,0])
    #(c,bins)=n.histogram(delta_coords[:,0],bins=10,new=True)
    print '<dx>','<dy>','<dz> [cm]'
    print n.mean(delta_coords,axis=0).squeeze()
    print 'sigma(x,y,z) [cm]'
    print n.std(delta_coords,axis=0).squeeze()
    p.xlabel("delta x [cm]")
    p.subplot(132)
    p.hist(delta_coords[:,1])
    p.xlabel("delta y [cm]")
    p.subplot(133)
    p.hist(delta_coords[:,2])
    p.xlabel("delta z [cm]")
#    p.subplot(234)
#    p.hist(dr_a,20)
#    p.xlabel("delta r [m]")
#    p.subplot(235)
#    p.hist(dt_a,20)
#    p.xlabel("delta l (r*theta) [m]")
    #p.subplot(236)
    #p.plot(r_a,t_a)
    #print r_a,t_a
    p.show()
else:
    aa = aas[0]
    antpos = n.array([aa.get_baseline(0,i,src='z') for i in range(len(aa.ants))])
    antpos /= 3.3356
    antpos += origin
    print antpos
    fmt1,fmt2 = fmts[0 % len(fmts)]
    p.plot(antpos[:,0],antpos[:,1], fmt1)
    p.plot([0],[0],fmt1)
    p.text(0,0,"base [ref base2]",fontsize=14)
    if not opts.p is None: 
        prms = aa.get_params()
        print "Circles are %f * dly"%(opts.z)
    print "base1 relative sighting coords"
    print "ant #\t distance [m]\t azimuth [deg]"
    for ant,(xa,ya,za) in enumerate(antpos):
        if not opts.p is None: 
            za=prms[str(ant)]['dly']
            print "off = ",za,
        hx,hy = opts.z*za*n.cos(th)+xa, opts.z*za*n.sin(th)+ya
        if za > 0: fmt3 = '#eeeeee'
        else: fmt3 = '#a0a0a0'
        if not opts.z is None: p.fill(hx,hy, fmt3,alpha=0.5)
        p.text(xa,ya,str(ant),fontsize=14)
        az = n.arctan(xa/ya)*180/m.pi
        if xa<0 and ya<0: az += 180 
        if xa<0 and ya>0: az += 360
        if xa>0 and ya<0: az += 180
        if az==n.NaN: az=0
        print "%d\t%d\t%d"%(ant,n.round(n.sqrt(xa**2+ya**2)).astype(int), int(az))
    p.axis('equal')
    p.show()
