import numpy as np, pylab as plt

fq = np.linspace(0.1, 0.2, 1024) # GHz
tau = 23 #ns

data = np.exp(-2 * np.pi * 1j * tau * fq)

plt.plot(fq, np.angle(data))
plt.plot(fq,np.unwrap( np.angle(data)))
plt.show()
