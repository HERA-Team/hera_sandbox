import aipy as a, numpy as n, pylab as p, ephem as e
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import export_beam
# plots the beam and the Fourier transform of the beam squared
sz = 200
d = 1./sz
img = a.img.Img(200,res=0.5)
X,Y,Z = img.get_top(center=(200,200))
shape0 = X.shape
X,Y,Z = X.flatten(),Y.flatten(),Z.flatten()
aa = a.cal.get_aa('psa6622_v001',n.array([.15, .18]))
aa.set_jultime(2456240.2)
peak = []
sample_ant = 1
ntop = n.array([X,Y,Z])  #note all beams are the same
print X
print Z
bmp_list = export_beam.beam_real(aa[sample_ant], ntop, shape0, 'x',sq=True)
bmp = bmp_list[0]
freq, fbmamp = export_beam.beam_fourier(bmp, d, 400)


rax = [-1,1,-1,1]
freq_pl = 5
mid = len(freq)/2
flim = (mid-freq_pl,mid+freq_pl)
freqax = [freq[flim[0]],freq[flim[1]], freq[flim[0]],freq[flim[1]]]
f_range = n.array(n.arange(flim[0],flim[1]))
m_range = n.arange(0, flim[1]-flim[0])


#peak.append(fbmamp[mid,mid])
fig = p.figure()
ax1 = fig.add_subplot(121)
im1 = ax1.imshow(bmp, interpolation='nearest', extent=rax)
ax1.set_xlabel('l')
ax1.set_ylabel('m')
cbar1 = p.colorbar(im1,fraction=0.046, pad=0.04)
ax2 = fig.add_subplot(122)
ax2.set_xlabel('u')
ax2.set_ylabel('v')
im2 = plt.imshow(fbmamp[f_range[:,None],f_range], interpolation='nearest', extent=freqax)
cbar2 = p.colorbar(im2,fraction=0.046, pad=0.04)

fig2 = p.figure()
u_range = n.linspace(-3,3,100)
f = export_beam.beam_interpol(freq,fbmamp)
z = export_beam.get_overlap(f, u_range,u_range)
#z = z/n.amax(n.abs(z))
U,V = n.meshgrid(u_range,u_range)
#z = abs(z)
ax3 = fig2.add_subplot(111,projection='3d')
ax3.set_xlabel('u')
ax3.set_ylabel('v')
surf = ax3.plot_surface(U,V,z,cmap=cm.coolwarm)
#p.colorbar(surf, shrink=0.5, aspect=5)


#ax4 = fig2.add_subplot(122)
#ax4.set_xlabel('u')
#ax4.set_ylabel('FT')
#for v in [0.,1.]:
#    z = export_beam.get_overlap(f,u_range,v)
#    curv, = p.plot(u_range,z,label='v= %f' % v)
#p.legend()
p.show()



#calculate Omp, Ompp
bm_list = export_beam.beam_real(aa[sample_ant], ntop, shape0, 'x',sq=False)
bm = bm_list[1]
dl,dm=2./400,2./400
Ompp = n.sum(bmp)*dl*dm
Omp = n.sum(bm)*dl*dm
Om = Omp*Omp/Ompp

print "Omp, Ompp, OmPrime=", Omp,Ompp, Om
