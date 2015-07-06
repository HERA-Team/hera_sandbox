#! /usr/bin/env python

import numpy
import pylab


time = numpy.linspace(0,10,1024) #s

freq = 1.1 #Hz

signal = numpy.exp(2j*numpy.pi*freq*time)
#signal = numpy.concatenate([signal,numpy.zeros(2048)])

ft_signal = numpy.fft.ifft(signal)
ft_freq = numpy.fft.fftfreq(signal.size,time[1]-time[0])


pylab.subplot(121)
pylab.plot(time,signal,'k-')
pylab.subplot(122)
pylab.semilogy(ft_freq,numpy.abs(ft_signal),'k-')

pylab.show()

