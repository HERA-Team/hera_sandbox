import numpy as np
import RotMeasTools as RMT
from pylab import *

nfig = 0
def plot3(x,y,labels,xlabels,title):
    global nfig
    figure(nfig)
    for i,s in enumerate(y):
        subplot(3,1,i+1)
        try: plot(x[i],s,label=labels[i])
        except(ValueError): plot(x,s,label=labels[i])
        legend()
        if i < len(labels)-1: xticks([])
        #yticks([])
    xlabel(xlabels)
    suptitle(title)
    #subplots_adjust(hspace=0)
    draw()
    nfig += 1
    

#GENERATE A FEW SIMULATED SPECTRA

fq = np.linspace(0.1,0.2,512)
RMsim = [10.,np.sqrt(200.)]
RMlabel = ['$10\ m^{-2}$','10$\sqrt{2}\ m^{-2}$','Sum']
QiUsim = []
for r in RMsim:
    QiUsim.append(RMT.gen_rm_spec(fq,r))
QiUsim.append(np.sum(QiUsim,axis=0)/2.)

plot3(fq,QiUsim,RMlabel,'Frequency (GHz)','Simulated Spectra $\Re\{Q + iU\}$')

##TEST THE REBIN/FFT METHOD

l2,QiUsim_l2 = [],[]
for i,spec in enumerate(QiUsim):
    l2i,spec_l2 = RMT.rebin_nu2lam2(fq,spec)
    l2.append(l2i)
    QiUsim_l2.append(spec_l2)
#plot3(l2,QiUsim_l2,RMlabel,'$\lambda^2\ (m^2)$ ','Rebinned into $\lambda^2$')

RM_samp,QiUsim_rm = [],[]
for i,spec in enumerate(QiUsim_l2): 
    QiUsim_rm.append(np.fft.fft(QiUsim_l2[i]))
    RM_samp.append(np.fft.fftfreq(len(spec),l2[i][1]-l2[i][0])*RMT.twopi)
#plot3(RM_samp,np.abs(QiUsim_rm),RMlabel,'RM $(m^{-2})$','FFTed')

#TEST THE DFT METHOD

RM_samp,W = RMT.RMTmat(fq)
for i,spec in enumerate(QiUsim):
    QiUsim_rm[i] = RMT.RMT(spec,W)
plot3(RM_samp,np.abs(QiUsim_rm),RMlabel,'$RM (m^2)$ ','DFT')

there_and_back_again = []
fuck_ups = []
iW = RMT.iRMTmat(fq)
for i,spec in enumerate(QiUsim_rm):
    there_and_back_again.append(RMT.iRMT(spec,iW))
    fuck_ups.append(there_and_back_again[-1] / QiUsim[i])

plot3(fq,there_and_back_again,RMlabel,'Frequency (GHz)','iDFT')

plot3(fq,np.abs(fuck_ups),RMlabel,'Frequency (GHz)','Recovered / Original [amp]')
plot3(fq,np.angle(fuck_ups),RMlabel,'Frequency (GHz)','Recovered / Original [phase]')

show()
