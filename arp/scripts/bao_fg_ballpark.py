#! /usr/bin/env python
import numpy as n, pylab as p

fq1 = n.linspace(.1,.2,1024)
fq2 = n.linspace(.7,.8,1024)

spec1 = (fq1/.150)**-1
spec2 = (fq2/.150)**-1

dspec1 = n.fft.fft(spec1)
dspec2 = n.fft.fft(spec2)

dly1 = n.fft.fftfreq(fq1.size, fq1[1]-fq1[0])
dly2 = n.fft.fftfreq(fq2.size, fq2[1]-fq2[0])

p.semilogy(n.fft.fftshift(dly1),n.fft.fftshift(n.abs(dspec1)))
p.semilogy(n.fft.fftshift(dly2),n.fft.fftshift(n.abs(dspec2)))

p.show()
