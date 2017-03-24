import numpy as n, pylab as p
N = 1
T = n.arange(-100,100,1)
for n in range(N):
	arr = n.random.normal(size=1000)
	for t in T:
		arshift = n.append(arr[-t:],arr[:t])
