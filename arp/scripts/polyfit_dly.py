#! /usr/bin/env python

import numpy as np, pylab as plt
dk_deta = 0.0004921824635812991 # h Mpc^-1 / ns

DN,UP = 3,4 # chan interval to fit over
DEG = 3 # degree of polynomial to fit
fq = np.arange(.1,.2,.000750) # GHz
NCHAN = fq.size
dly = np.fft.fftshift(np.fft.fftfreq(NCHAN-DN-UP, fq[1]-fq[0])) # ns

loss, wgt = 0., 0.
for cnt in xrange(1000):
    eor = np.random.normal(size=NCHAN) + 1j * np.random.normal(size=NCHAN)
    mdl = np.zeros_like(eor)
    for c in xrange(DN,NCHAN-UP):
        pr = np.polyfit(fq[c-DN:c+UP], eor[c-DN:c+UP].real, deg=DEG)
        pi = np.polyfit(fq[c-DN:c+UP], eor[c-DN:c+UP].imag, deg=DEG)
        mdl[c] = np.polyval(pr,fq[c]) + 1j * np.polyval(pi,fq[c])

    eor_dly = np.fft.fftshift(np.fft.ifft(eor[DN:-UP]))
    mdl_dly = np.fft.fftshift(np.fft.ifft(mdl[DN:-UP]))
    loss += np.abs(eor_dly-mdl_dly)**2 / np.abs(eor_dly)**2
    wgt += 1.
plt.plot(dly*dk_deta, loss/wgt)
plt.xlabel('k [h Mpc^-1]')
plt.ylabel('Signal Loss')
plt.show()

#plt.figure(1)
#plt.plot(fq, eor)
#plt.plot(fq, mdl)
#plt.plot(fq, eor - mdl)

#plt.figure(2)
#plt.plot(dly, eor_dly)
#plt.plot(dly, mdl_dly)
#plt.plot(dly, eor_dly - mdl_dly)

