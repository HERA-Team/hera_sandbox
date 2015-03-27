import aipy as a, numpy as n, pylab as p, ephem as e
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import export_beam

sz=200
d=1./sz
img=a.img.Img(200,res=0.5)
X,Y,Z=img.get_top(center=(200,200))
shape0=X.shape
X,Y,Z=X.flatten(),Y.flatten(),Z.flatten()
aa=a.cal.get_aa('psa6622_v001',n.array([.15]))

fig = p.figure()
aa.set_jultime(2456240.2)
peak = []
sample_ant = 1
ntop = n.array([X,Y,Z])  #note all beams are the same
bmp = export_beam.beam_real(aa[sample_ant], ntop, shape0, 'x')
freq, fbmamp = export_beam.beam_fourier(bmp, d, 400)
rax = [-1,1,-1,1]
freq_pl = 5
mid = len(freq)/2
flim = (mid-freq_pl,mid+freq_pl)
freqax = [freq[flim[0]],freq[flim[1]], freq[flim[0]],freq[flim[1]]]
f_range = n.array(n.arange(flim[0],flim[1]))
m_range = n.arange(0, flim[1]-flim[0])
#peak.append(fbmamp[mid,mid])
ax1 = fig.add_subplot(1,2,1)
im1 = ax1.imshow(bmp, interpolation='nearest', extent=rax)
ax1.set_xlabel('l')
ax1.set_ylabel('m')
cbar1 = p.colorbar(im1,fraction=0.046, pad=0.04)
ax2 = fig.add_subplot(1,2,2)
ax2.set_xlabel('u')
ax2.set_ylabel('v')
im2 = plt.imshow(fbmamp[f_range[:,None],f_range], interpolation='nearest', extent=freqax)
cbar2 = p.colorbar(im2,fraction=0.046, pad=0.04)

#plt.colorbar()
#p.xlabel('u',size=14)
#p.ylabel('v',size=14)
p.show()
