#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import optparse, sys,time
from numpy import mgrid
from scipy.interpolate import interp1d
from scipy.interpolate import interp2d
from smooth import smooth
from cosmo_units import *
def length(A,B):
    return n.sqrt(n.dot(A,B))

p.ion()
#p.figure()
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True)
opts,args = o.parse_args(sys.argv)

#aa = a.cal.get_aa(opts.cal, n.array([.150]))
nchan = 50 #bandwidth of 3.6MHz
fo = 0.1738 #band frequency
df = 73.24 * 1e-6 #frequency resolution
B = nchan * df
aa = a.cal.get_aa(opts.cal, n.array([fo]))
F = n.arange(fo,fo+B,df)
kpar = n.linspace(eta2kparr(1/(B*1e9),f212z(fo*1e9)),
    eta2kparr(len(F)/(B*1e9)/2.,f212z(fo*1e9)),num=len(F)/2.)
h = a.healpix.HealpixMap(nside=64)
#h = a.healpix.HealpixMap(nside=32)
h.map = n.random.normal(size=h.map.size)
h.set_interpol(True)

#im = a.img.Img(size=200, res=.5)
SZ = 200
im = a.img.Img(size=SZ, res=0.5)
tx,ty,tz = im.get_top()
invalid = tx.mask.copy()
tx = tx.filled(0).flatten()
ty = ty.filled(0).flatten()
tz = tz.filled(0).flatten()
#generate a uv coordinate grid (need dimensions)
umax = im.size
du  = im.res
#U,V = mgrid[-umax:umax:du,-umax:umax:du]
#urange = n.arange(-umax,umax,du)
#Uf,Vf = mgrid[-umax:umax:du/4,-umax:umax:du/4]
#interpolate on the uv grid
#delayish stuff
resp = aa[0].bm_response((tx,ty,tz), pol='y')**2
resp.shape = invalid.shape
resp = n.where(invalid, 0, resp/resp[0,0])


Np = 2 #number of free parameters in smooth
window_len = len(F) - Np
lsts = []
uvdata0,uvdata1,uvdata2,uvdata3 = [],[],[],[]
PLOT = None
S,Sm = [],[]
Si_all = []
Si2_all = []
S_all = []
Sd_all = []
uv_subs = []
Sm_all = []
for ha in n.arange(-n.pi/16, n.pi/16, 2*n.pi / (24*60)):
    print ha
    lsts.append(ha)
    ex,ey,ez = im.get_eq(ra=ha, dec=aa.lat)
    ex = ex.filled(1).flatten()
    ey = ey.filled(0).flatten()
    ez = ez.filled(0).flatten()
    Tsky = h[ex,ey,ez]
    Tsky.shape = resp.shape
    img = Tsky * resp
    img = a.img.recenter(img, (img.shape[0]/2, img.shape[1]/2))
    uv = n.fft.fft2(img)
#    UV = lambda U,V: interp2d(urange,urange,n.real(uv)) + interp2d(urange,urange,n.imag(uv))*1j
    
    plt_uv = a.img.recenter(n.abs(uv), n.array(uv.shape)/2)[SZ:SZ+64,SZ:SZ+64]
#    plt_uvf = a.img.recenter(n.abs(UV(Uf,Vf)), n.array(Uf.shape)/2)[SZ*4:SZ*4+64*4,SZ*4:SZ*4+64*4]
    i,j = 0,6
#    for i in range(len(aa)):
#        for j in range(i,len(aa)):
#            if i==j:continue
    u,v,w = aa.gen_uvw(i,j)
    Sd = []
    Sa = []
    Sp = []
    Si = []
    Si_2 = []
    S = []
    di = 1
    for f in F:
        U = u*f/fo/du
        V = v*f/fo/du
        ui = n.round(U).astype(n.int)
        vj = n.round(V).astype(n.int)
        dU,dV = n.meshgrid(n.array(n.arange(ui-di,ui+di+1)-U).squeeze(),
            n.array(n.arange(vj-di,vj+di+1).squeeze()-V).squeeze())
        dR = n.sqrt(dU**2 + dV**2)
        Sa.append(n.sum(n.abs(uv[ui-di:ui+di+1,vj-di:vj+di+1])/dR/(n.sum(1/dR))))
        Sp.append(n.sum(n.angle(uv[ui-di:ui+di+1,vj-di:vj+di+1])/dR/(n.sum(1/dR))))
        Si.append(Sa[-1]*n.exp(1j*Sp[-1]))
        Si_2.append(n.sum(uv[ui-di:ui+di+1,vj-di:vj+di+1]/dR/(n.sum(1/dR))))
        Sd.append(uv[ui-1:ui+2,vj-1:vj+2])
        S.append(uv[ui,vj].squeeze())
    Si2_all.append(n.array(Si_2))
    S_all.append(n.array(S))
    Si_all.append(n.array(Si))
    Sd_all.append(n.array(Sd))
    Sm_all.append(smooth(n.abs(S),window_len=window_len)*n.exp(1j*smooth(n.angle(S),window_len=window_len)))
    print(len(Si_all))
#    sys.exit()
#        UV = uv[
        

#    fi = n.sort(n.array(fi))
#    fi[0],fi[-1] = f.min(),f.max()
#    fi_u[0],fi_u[-1] = f.min(),f.max()
#    fi_v[0],fi_v[-1] = f.min(),f.max()
#    uvRi = interp2d(ui_i,vj_i,n.real(uv_sub))
#    uvIi = interp2d(ui_i,vj_i,n.imag(uv_sub))
#    uvi = lambda u,v: uvRi(u,v) + uvIi(u,v)*1j
    


#    S.append(uv[ui,vj])
#    Sm.append(smooth(n.real(S[-1]),window_len=10,window='hanning') + \
#        smooth(n.imag(S[-1]),window_len=10,window='hanning')*1j)
#    Sm.append(n.array([uvi(U,V) for (U,V) in zip(u.squeeze()/du,v.squeeze()/du)]))
#    uv_subs.append(uv_sub)
#    print uv_sub.shape
#    F.append(f)
#    Fi.append(fi)
    if PLOT is None:
        PLOT = {}
        PLOT['bl'] = p.plot(F,n.abs(S),'k',alpha=0.5)[0]
        PLOT['bl_i'] = p.plot(F,n.abs(Si),'.b')[0]
        PLOT['avg2'] = p.plot(F,n.abs(n.average(Si2_all,axis=0)))[0]
        PLOT['avg'] = p.plot(F,n.abs(n.average(Si_all,axis=0)))[0]
    else:
        PLOT['bl'].set_ydata(n.abs(S))
        PLOT['bl_i'].set_ydata(n.abs(Si))
        PLOT['avg'].set_ydata(n.abs(n.average(Si_all,axis=0)))
        PLOT['avg2'].set_ydata(n.abs(n.average(Si2_all,axis=0)))
    p.draw()
n.save('delay_noise',n.average(Sm_all,axis=0))
p.ioff()    
#            umag = n.sqrt(u[:inds.shape[0]]**2+v[:inds.shape[0]]**2)
#            umag = 
#            print u,v,umag
#    #            umagi = n.sqrt(u[ui]**2+v[vj]**2)
#    #            print umagi.shape
#    #            Si = interp1d(umagi,S,kind='nearest')
#            if PLOT is None:
#                PLOT=  {}
#                p.subplot(221)
#                PLOT['img'] = p.imshow(n.abs(img),vmax=img.max(),vmin=0,origin='lower')
#                p.subplot(222)
#                PLOT['uv'] = p.imshow(plt_uv,origin='lower')
#    #                p.subplot(223)
#    #                PLOT['uvf'] = p.imshow(plt_uvf,origin='lower')
#                p.subplot(224)
#                PLOT['bls'] = p.plot(fbl,n.abs(S),'.')[0]
#            else:
#                print len(umag),len(S)
#                PLOT['bls'].set_ydata(S)
#                PLOT['bls'].set_xdata(fbl)
#            p.show()
#            sys.exit()
#    p.draw()
#    PLOT['img'].set_data(n.abs(img))
#    PLOT['uv'].set_data(plt_uv)
#    p.draw()        
    #for each baseline, sample uv plane
    
#    plt_uv0 = a.img.recenter(n.abs(uv), n.array(uv.shape)/2)[SZ:SZ+64,SZ:SZ+64]
#    plt_uv1 = a.img.recenter(n.angle(uv), n.array(uv.shape)/2)[SZ:SZ+64,SZ:SZ+64]
#    if PLOT is None:
#        PLOT = {}
#        p.subplot(221)
#        PLOT['img0'] = p.imshow(img, vmax=img.max(), vmin=-img.max(), origin='lower')
#        p.subplot(222)
#        PLOT['img1'] = p.imshow(n.abs(img), vmax=img.max(), vmin=0, origin='lower')
#        p.subplot(223)
#        PLOT['uv0'] = p.imshow(plt_uv0, origin='lower')
#        p.subplot(224)
#        PLOT['uv1'] = p.imshow(plt_uv1, origin='lower')
#    else:
#        PLOT['img0'].set_data(img)
#        PLOT['img1'].set_data(n.abs(img))
#        PLOT['uv0'].set_data(plt_uv0)
#        PLOT['uv1'].set_data(plt_uv1)
#    uvdata0.append(uv[0,:uv.shape[1]/2])
#    uvdata1.append(uv[16,:uv.shape[1]/2])
#    uvdata2.append(uv[32,:uv.shape[1]/2])
#    uvdata3.append(uv[:uv.shape[0]/2,0])
#    p.draw()

#uvdata0 = n.array(uvdata0); uvdata0 /= n.abs(uvdata0)
#uvdata1 = n.array(uvdata1); uvdata1 /= n.abs(uvdata1)
#uvdata2 = n.array(uvdata2); uvdata2 /= n.abs(uvdata2)
#uvdata3 = n.array(uvdata3); uvdata3 /= n.abs(uvdata3)
#lsts = n.array(lsts)
#
#p.ioff()
#p.clf()
#nacc = n.arange(uvdata0.shape[0])+1
#for step in [8,16,32,64]:
#    p.subplot(2,2,1)
#    p.loglog(nacc*60, n.abs(n.cumsum(uvdata0[:,step]) / nacc), label='uv=%2d, 0'%(step*.5))
#    p.subplot(2,2,2)
#    p.loglog(nacc*60, n.abs(n.cumsum(uvdata1[:,step]) / nacc), label='uv=%2d, 8'%(step*.5))
#    p.subplot(2,2,3)
#    p.loglog(nacc*60, n.abs(n.cumsum(uvdata2[:,step]) / nacc), label='uv=%2d,16'%(step*.5))
#    p.subplot(2,2,4)
#    p.loglog(nacc*60, n.abs(n.cumsum(uvdata3[:,step]) / nacc), label='uv= 0,%2d'%(step*.5))
#
#for i in range(4):
#    p.subplot(2,2,i+1)
#    p.loglog(nacc*60, 3*nacc**-.5, 'k:')
#    p.xlim(1e2, 1e4)
#    p.ylim(1e-2, 1e1)
#    p.legend(loc='best')
#    p.grid()
#
#p.show()
