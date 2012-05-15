import numpy as np
from aipy import dsp
import RotMeasTools as RMT
from pylab import *
import pspec

c = 0.3 #m/ns
twopi = 2.*np.pi

N = 512
fq = np.linspace(0.1,0.2,N)
L2 = RMT.better_guess_l2(fq)

Rms,W = RMT.RMTmat(fq)
plot_me = np.zeros((N,N),dtype=np.complex)
delays = np.fft.fftshift(np.fft.fftfreq(N,fq[1]-fq[0]))
wgt = dsp.gen_window(N,window='blackman-harris')
Rm0 = 5.
for i,Rm in enumerate(Rms): 
    QiUi = RMT.gen_rm_spec(fq,Rm+Rm0)
    plot_me[:,i] = np.fft.fftshift(np.fft.ifft(QiUi*wgt))
    
plot_me1 = plot_me*np.conjugate(plot_me)
plot_me1 = plot_me1.real
plot_me1 /= np.max(plot_me1)

figure(0)
imshow(10.*np.log10(plot_me1),aspect='auto',interpolation='nearest',
    extent=[delays[0],delays[-1],Rms[-1],Rms[0]],
    vmax=0,vmin=-100)
ylabel('Rotation Measure [$m^{-2}$]')
xlabel('Delay [ns]')
title('$|W|^2$ [dB]')
colorbar()
draw()

centroid = np.sum(wgt * delays * plot_me1,axis=0)/np.sum(wgt*plot_me1,axis=0)
width = np.sum(wgt * (delays**2) * plot_me1,axis=0)/np.sum(wgt*plot_me1,axis=0)
width -= centroid**2
width = np.sqrt(np.abs(width))

figure(1)
Nslice = 8 
for i in range(Nslice-1):
    if i < Nslice/2 - 1: continue
    ii = int(N*((i+1.)/Nslice))
    plot(delays,10.*np.log10(plot_me1[:,ii]),label='$\Phi=%2.2f\ m^{-2}$'%Rms[ii])
legend()
title('Slices through $|W|^2$ at constant $\Phi$')
xlabel('delay [ns]')
ylabel('Amplidude [dB]')
xlim([0,delays[-1]])
ylim([-100,0])
draw()

#k = pspec.dk_deta(8.)*delays
#k = k[N/2:]
#Pk = (plot_me * np.conjugate(plot_me)).real 
#Pk *= pspec.jy2T(0.150)**2
#Pk *= pspec.X2Y(8.) *((fq[1]-fq[0])* 1e-9) / 2.
#Pk *= 1e6
#
#figure(2)
#for i in range(Nslice-1):
#    if i < Nslice/2 - 1: continue
#    ii = int(N*((i+1.)/Nslice))
#    k3pk = Pk[N/2:,ii] * k**3 / (2.*np.pi**2)
#    loglog(k,k3pk,label='$\Phi=%2.2f\ m^{-2}$'%Rms[ii])
#legend(loc='lower right')
#title('$\Delta^2(k)$ for 1 Jy source [$mK^2$]')
#xlabel('$k\ h$Mpc$^{-1}$')
#draw()

#figure(2)
#suptitle('Centroid and width of $W_{ij}$')
#subplot(211)
#plot(Rms,centroid)
#ylabel('Centroid [ns]')
#subplot(212)
#plot(Rms,width)
#ylabel('Width [ns]')
#xlabel('Rotation Measure [$m^{-2}$]')
#draw()

show()
