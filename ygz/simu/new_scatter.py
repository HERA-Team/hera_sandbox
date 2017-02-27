import numpy as n, matplotlib.pyplot as p
from mpl_toolkits.axes_grid.axislines import SubplotZero
import pandas as pd


markers = matplotlib.markers.MarkerStyle.markers

df = pd.DataFrame.from_csv('corr_res.csv')
sep = str(df['sep'])+':'+str(df['sep2'])
Npts = df.shape[0]
fig, ax = p.subplots()
ax = fig.add_subplot(211)
for i in xrange(Npts):
	marker = markers[i]
	ax.scatter(n.abs(df['dT'][i]), df['peak'][i], label=sep[i], maker=marker)
plt.legend()

ax.grid()
ax.set_xlabel('Time Delay [Sidereal Day]')
ax.set_ylabel('Correlation [Normalized to 1]')

ax = fig.add_subplot(212)
for i in xrange(Npts):
	marker = markers[i]
	ax.scatter(n.abs(df['dT'][i]), df['mult']*df['peak'][i], label=sep[i], maker=marker)
plt.legend()

ax.grid()
ax.set_xlabel('Time Delay [Sidereal Day]')
ax.set_ylabel('Correlation [Normalized to 1]')
p.show()