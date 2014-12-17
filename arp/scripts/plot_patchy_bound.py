#! /usr/bin/env python
import numpy as n, pylab as p

def Tcmb(z):
    return 2.7e3 * (1+z)

def Tgas(z):
    return Tcmb(z) * ((1+z)/150.)

def k3pk_patchy(T_b, x_i=0.5, kmax_kmin=2):
    x_H = 1-x_i
    return (x_H - x_H**2) * T_b**2 / n.log(kmax_kmin)

limits = {}
#limits['.05-.1'] = (n.arange(.05, .1, .001), 3100.,'c')
#limits['.1-.2'] = (n.arange(.1, .2, .001), 1450.,'m')
#limits['.2-.4'] = (n.arange(.2, .4, .001), 3460.,'b')
##limits['.4-.5'] = (n.arange(.4, .52, .001), 8210., 'y')

# Errors are 2sigma
if True: # PSA64 limits
    limits['0.130'] = (0.130, 1002.06079, 368.28219, 'k')
    limits['0.156'] = (0.156, 1387.49041, 711.60973, 'k')
    limits['0.182'] = (0.182, 2172.47976, 148.88083, 'k')
    limits['0.208'] = (0.208, 1411.65841, 177.60727, 'k')
    limits['0.234'] = (0.234,  822.62623, 673.71216, 'k')
    limits['0.259'] = (0.259,  842.20637, 195.00152, 'k')
    limits['0.285'] = (0.285,  710.55152, 483.71724, 'k')
    limits['0.311'] = (0.311,  819.47308, 399.65856, 'k')
    limits['0.337'] = (0.337, 1529.75429, 371.16516, 'k')
    limits['0.363'] = (0.363,    4.96123,1153.46259, 'k')
    limits['0.389'] = (0.389,   39.53246, 819.52705, 'k')
    limits['0.415'] = (0.415, -438.08045, 463.34828, 'k')
    limits['0.441'] = (0.441, 1024.38682,2163.02153, 'k')
    limits['0.466'] = (0.466, 2856.04752,4320.85202, 'k')
    limits['0.492'] = (0.492, 1574.46938,2494.22128, 'k')

#kwids = n.arange(1.1,3000.,.1)
#kwids = 10**n.arange(0.1, 6, .1)
#kcens = 10**n.arange(-2, 2, .1)
kmins = 10**n.arange(-3, 2.1, .03)
kmaxs = 10**n.arange(-3, 2.1, .03)
#for xi,sty in [(.5,'-'), (.3,'--'), (.1,'-.'), (.03,':')]:
for cnt,(xi,sty) in enumerate([(.1,'-'), (.3,'--'), (.5,'-.')]):
    T_b_lim = []
    xs = [kmins] * len(kmaxs); xs = n.array(xs).transpose()
    ys = [kmaxs] * len(kmins); ys = n.array(ys)
    print xs.shape, ys.shape
    for kmax in kmaxs:
        T_b_lim_kwid = []
        for kmin in kmins:
            tmp = []
            for L in limits:
                ks, k3pk, err, clr = limits[L]
                #if ks < kcen - kwid or ks > kcen + kwid:
                if ks < kmin or ks > kmax:
                    tmp.append(400)
                    continue
                lim = k3pk + err
                ks = n.linspace(ks-.014, ks+.14, 100)
                #avg = n.average(k3pk_patchy(ks, 1., xi, kwid))
                avg = n.average(k3pk_patchy(1., xi, kmax/kmin))
                tmp.append(n.sqrt(lim / avg))
            T_b_lim_kwid.append(min(tmp+[400]))
        T_b_lim.append(T_b_lim_kwid)
    T_b_lim = n.array(T_b_lim)
    print T_b_lim.shape
    #p.plot(Rs, T_b_lim, 'k'+sty)
    #p.fill_between(kwids, T_b_lim, 1000*n.ones_like(kwids), facecolor='k', alpha=0.2)
    p.subplot(1,3,cnt+1)
    p.imshow(T_b_lim, vmax=400, vmin=0, interpolation='nearest', aspect='auto', origin='lower',
        #extent=(-3,2,-3,2))
        extent=(1e-3,1e2,1e-3,1e2))
    p.xscale('log')
    p.yscale('log')
    p.xlim(1e-3,1e0)
    p.ylim(1e-1,1e2)
    #p.fill_between([-3,3],[-3,3],[-3,-3], facecolor='w', alpha=1)
    if xi == 0.5: p.title('$x_i=%3.1f$'%xi)
    else: p.title('$x_i=%3.1f,%3.1f$' % (xi,1-xi))
    #p.contourf(ys, xs, T_b_lim, [50,100,150,200,250,300,350,400])
    p.fill_between([1e-3,1e3],[1e-3,1e3],[1e-3,1e-3], facecolor='w', alpha=1)
    #p.colorbar()
    p.xlabel(r'$k_{\rm min}$ [$h^{-1}\ {\rm Mpc}$]', size=14)
    if cnt == 0: p.ylabel(r'$k_{\rm max}$ [$h^{-1}\ {\rm Mpc}$]', size=14)
    p.grid()
    if cnt in [1,2]:
        p.plot([.2], [30.], 'ko')
        p.plot([.1], [30.], 'k^')
        #p.plot([.3], [10.], 'ko')
        #p.plot([.15], [10.], 'k^')
    #if cnt == 2: p.colorbar()
    
p.show()
#p.semilogx([1., 1000], [400, 400], 'k--')
p.fill_between([1., 1e8], [400, 400], [1000,1000], facecolor='m', alpha=.5)
p.semilogx([1., 1e8], [30, 30], 'k--')
p.ylim(0,500)
p.xlim(1e1,1e6)
p.ylabel(r'$|\langle T_b\rangle|$ [${\rm mK}$]', size=14)
p.xlabel(r'$k_{\rm max}/k_{\rm min}$ [$h^{-1}\ {\rm Mpc}$]', size=14)
p.grid()

for L in limits:
    ks, k3pk, err, clr = limits[L]
    lim = k3pk + err
    kmax,kmin = 10.,.1
    for xi in [.3, .5]:
        avg = n.average(k3pk_patchy(1., xi, kmax/kmin))
        print L, xi, n.sqrt(lim / avg)

if True: # set true to get the scalings for the lidz curves
    import glob, scipy.interpolate, re, os, capo as C, scipy.special

    re_z = re.compile(r'power_21cm.*_z(\d+\.\d+).*\.dat')
    # xi = .54, .71
    files_norm = [glob.glob('lidz_mcquinn_k3pk/power*%s*.dat' % s)[0] for s in ['7.3','7.0']]
    # xi = .77, .63, .51, .41
    # xi = .51, .77
    #files_hmass = [glob.glob('lidz_mcquinn_k3pk/hmass/power*%s*.dat' % s)[0] for s in ['7.1','7.3','7.4','7.6']]
    files_hmass = [glob.glob('lidz_mcquinn_k3pk/hmass/power*%s*.dat' % s)[0] for s in ['7.4','7.1']]

    #for cnt,filename in enumerate(glob.glob('lidz_mcquinn_k3pk/*dat')):
    for cnt,files in enumerate([files_norm, files_hmass]): #glob.glob('lidz_mcquinn_k3pk/hmass/power_21cm_hmass_*dat')):
        clr = 'cm'[cnt]
        for cnt,filename in enumerate(files):
            sym = ['-','--','-.',':'][cnt]
            print 'Reading', filename, clr + sym
            d = n.array([map(float, L.split()) for L in open(filename).readlines()])
            ks, pk = d[:,0], d[:,1]
            z_file = float(re_z.match(os.path.basename(filename)).groups()[0])
            z = C.pspec.f2z(.160)
            k3pk = ks**3 / (2*n.pi**2) * pk
            mdl = scipy.interpolate.interp1d(ks, k3pk, kind='linear')
            Ls = limits.keys(); Ls.sort()
            ks = n.array([limits[L][0] for L in Ls])
            k3pk = n.array([limits[L][1] for L in Ls])
            err = n.array([limits[L][2] for L in Ls])
            lvl = mdl(ks)
            Ts = n.sqrt((k3pk+err)/lvl)
            print ks
            print n.around(Ts)
            Tmin = Ts.min()
            print '-'*10
            p.plot([1., 1000], [Tmin, Tmin], clr+sym)

p.show()
import sys; sys.exit(0)

def V_to_pk(V, L, nbins=200):
    dL = L / V.shape[0]
    kmin,kmax = 2*n.pi/L, 2*n.pi/dL
    P_k = n.abs(n.fft.ifftn(V)/dL**3)**2 * L**3
    print P_k[0,0,0], L, dL
    kx = 2*n.pi * n.fft.fftfreq(V.shape[0], d=dL); kx.shape = (kx.size,1,1)
    ky = 2*n.pi * n.fft.fftfreq(V.shape[1], d=dL); ky.shape = (1,ky.size,1)
    kz = 2*n.pi * n.fft.fftfreq(V.shape[2], d=dL); kz.shape = (1,1,kx.size)
    k = n.sqrt(kx**2 + ky**2 + kz**2)
    print k[0,0,0], k[1,0,0]
    hist_sum,bins = n.histogram(k, range=(kmin,kmax), bins=nbins, weights=P_k)
    hist_wgt,bins = n.histogram(k, range=(kmin,kmax), bins=nbins, weights=n.ones_like(P_k))
    hist = hist_sum / n.where(hist_wgt > 0, hist_wgt, 1)
    return 0.5*(bins[:-1] + bins[1:]), hist
    
kmin,kmax = .03, 3. # h Mpc^-1
L,dL = 2*n.pi/kmin, 2*n.pi / kmax
SZ = int(L / dL)

print L, dL
print SZ


# Add in bubbles
def ionize(V, L, xi=0.5, avg_log10R=n.log10(2*n.pi/0.25), sig_log10R=.1):
    x,y,z = n.indices(V.shape)
    print n.average(V)
    while n.average(V) > xi:
        r = 10**n.random.normal(loc=avg_log10R, scale=sig_log10R)
        r /= L / V.shape[0]
        print n.average(V), r
        x0,y0,z0 = n.random.randint(V.shape[0], size=len(V.shape))
        print x0, y0, z0
        V = n.where((x-x0)**2+(y-y0)**2+(z-z0)**2 <= r**2, 0, V)
    print n.average(V)
    return V

def ionize_square(V, L, xL):
    x,y,z = n.indices(V.shape, dtype=n.float32)
    xL = xL / (L / V.shape[0])
    x = n.where((x/xL) % 2 < 1, 0, 1)
    y = n.where((y/xL) % 2 < 1, 0, 1)
    z = n.where((z/xL) % 2 < 1, 0, 1)
    xyz = (x + y + z) % 2
    return V * xyz
        
#for xL in n.arange(n.log10(kmin), n.log10(kmax), .1):
#    V = n.ones((SZ,SZ,SZ), dtype=n.float32)
#    xL = 2*n.pi/(10**xL)
#    print xL, dL
#    V = ionize_square(V, L, xL)
#    print n.average(V)
#    ks, pk = V_to_pk(V, L)
#    p.plot(ks, ks**3/(2*n.pi**2) * pk)
#    print 'done'
V = n.ones((SZ,SZ,SZ), dtype=n.float32)
V = n.random.normal(scale=V)
V *= (Tcmb(8.) - Tgas(8.))
ks, pk = V_to_pk(V, L)
p.plot(ks, ks**3/(2*n.pi**2) * pk)

p.show()

