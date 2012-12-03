import numpy as np
import RotMeasTools as RMT
from pylab import *

N = 15

fq = np.linspace(0.12,0.18,N)

W = RMT.gen_RMtau_ker(fq)
Wt= RMT.gen_RMtau_ker(fq,inv=True).T

print W.shape,Wt.shape

KER = np.dot(Wt,W)
TEST = np.dot(W,Wt)

figure()
subplot(311)
imshow(np.abs(KER))
subplot(312)
INV = np.linalg.inv(KER)
imshow(np.abs(INV))
subplot(313)
one = np.dot(INV,KER)
imshow(np.abs(one))
draw()

figure()
TEST2 = np.zeros((N,N))
TEST2[N/2,N/2] = 1.
TEST2 = TEST2.flatten()
CONV2 = np.dot(KER,TEST2)
dCONV = np.dot(INV,CONV2)
dCONV = dCONV.reshape((N,N))
CONV2 = CONV2.reshape((N,N))
TEST2 = TEST2.reshape((N,N))

subplot(311)
imshow(np.abs(TEST2))
subplot(312)
imshow(np.abs(CONV2))
subplot(313)
imshow(np.abs(dCONV))
draw()


new_spec = np.dot(np.dot(INV,Wt),np.ones(N))
figure()
plot(np.ones(N))
plot(new_spec.real)
draw()

show()
