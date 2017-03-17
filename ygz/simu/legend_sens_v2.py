import numpy as np, matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from mpl_toolkits.axes_grid.axislines import SubplotZero
import pandas as pd
from itertools import cycle
from matplotlib import colors as mcolors
import seaborn as sns 
sns.set_context("paper", font_scale=2.5)

FILE = 'corr_res.csv'
#FILE = 'HERA_350_core_pm300.csv'
LEGEND = True

#markers = matplotlib.markers.MarkerStyle.markers
markers = cycle(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'H', 'D', 'd', 'x'])
print markers
df = pd.DataFrame.from_csv(FILE)
Npts = df.shape[0]
#Npts = 1000
C = np.random.random(Npts)
f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
# f.set_figheight(10);  f.set_figwidth(12)
markerszip = [markers.next() for i in xrange(Npts)]
Cc = cycle(['r', 'g', 'b', 'y'])
Czip = [Cc.next() for i in xrange(Npts)]
peakmax = np.amax(df['peak'])
colors = mcolors.cnames.keys()
Cz = [np.random.choice(colors) for i in xrange(Npts)]
#marker = markers[i%len(markers)]
#label = str(df['sep'])+':'+str(df['sep2'])
# for dt, dp, m, c in zip(np.abs(df['dT']), df['peak'], markerszip, C):
# 	ax.scatter(dt, dp, marker=m, c=c, cmap=cm.jet)
for i in xrange(Npts):
	marker = markerszip[i]
	label = str(df['sep'][i])+':'+str(df['sep2'][i]) if i<Npts/2 else None
	#print C[i]
	ax.scatter(np.abs(df['dT'][i]), df['peak'][i]/peakmax, label=label, s=65, marker=marker, c=Cz[i])
plt.setp(ax.get_xticklabels(), visible=False)
if LEGEND:
	legend = ax.legend(ncol=4, scatterpoints=1, frameon=False,fontsize = 'x-small')
	legend.get_frame().set_facecolor('none')

#ax.grid()
#ax.set_xlabel('Time Delay [Sidereal Day]')
ax.set_ylabel('Correlation')
ax.set_xlim([-0.01, 0.21])

mult1010 = float(np.amax(df['mult']))
mult1010 = float(11025)
for i in xrange(Npts):
	marker = markerszip[i]
	label = str(df['sep'][i])+':'+str(df['sep2'][i]) if i>=Npts/2 else None
	ax2.scatter(np.abs(df['dT'][i]), (df['mult']*df['peak']/mult1010)[i]/peakmax, s=65, c=Cz[i], label=label, marker=marker)
#plt.legend()
if LEGEND:
	legend = ax2.legend(ncol=4, scatterpoints=1, frameon=False, fontsize = 'x-small')
	legend.get_frame().set_facecolor('none')

#ax2.grid()
ax2.set_xlabel('Time Delay [Sidereal Day]')
ax2.set_ylabel('Weighted Correlation')
#plt.tight_layout()
f.subplots_adjust(hspace=0.05)
plt.show()


print "========= Statistics of sensitibity contribution =========="
df['Theta'] = df['peak']*df['mult']
df['Theta'] /= np.amax(df['Theta'])
def get_imp(df, Theta_min=0.4):
	dft = df.where(df['Theta']>Theta_min)
	dfeq = dft.where(dft['sep']==dft['sep2'])
	dfnq = dft.where(dft['sep']!=dft['sep2'])
	
	totalsumsq = np.sum(dft['Theta']**2)
	eqsumsq = np.sum(dfeq['Theta']**2)
	totalsens = np.sum(dft['Theta'])/totalsumsq
	eqsens = np.sum(dfeq['Theta'])/eqsumsq
	improve = (totalsens-eqsens)/eqsens
	return improve


