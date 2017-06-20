import numpy as np, matplotlib.pyplot as plt
import matplotlib
from matplotlib import cm
from mpl_toolkits.axes_grid.axislines import SubplotZero
import pandas as pd
from itertools import cycle
from matplotlib import colors as mcolors
from matplotlib.font_manager import FontProperties
import seaborn as sns 
from pylab import MaxNLocator
sns.set_context("paper", font_scale=1.5)
sns.set(style="ticks", color_codes=False,font='DejaVu Serif', font_scale=1.5)
plt.rc('axes', linewidth=1.5)
font = FontProperties()
# font.set_weight('bold')
# font.set_size('large')
FILE = 'corr_res.csv'
#FILE = 'HERA_350_core_pm300.csv'
LEGEND = False

#markers = matplotlib.markers.MarkerStyle.markers
markers = cycle(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'd', 'x', 'H', 'D'])
print markers
df = pd.DataFrame.from_csv(FILE)
Npts = df.shape[0]
#Npts = 1000
C = np.random.random(Npts)
#f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
# f.set_figheight(10);  f.set_figwidth(12)
X = np.arange(0, 4, 1).astype(np.int)
Y = np.arange(0, 3, 1).astype(np.int)
#X, Y = np.meshgrid(X, Y)
f = plt.figure()
axleg1 = plt.subplot2grid((3, 2), (0, 0))
axleg2 = plt.subplot2grid((3, 2), (0, 1))
ax = plt.subplot2grid((3, 2), (1, 0), colspan=2)
ax2 = plt.subplot2grid((3, 2), (2, 0), colspan=2)

markerszip = [markers.next() for i in xrange(Npts)]
Cc = cycle(['r', 'g', 'b', 'y'])
#Czip = [Cc.next() for i in xrange(Npts)]
Czip = ['r', 'g', 'b', 'y']
peakmax = np.amax(df['peak'])
colors = mcolors.cnames.keys()
Cz = [np.random.choice(colors) for i in xrange(Npts)]
ix = 0
for x in X:
	iy = 0
	czip = Czip[ix]
	for y in Y:
		cedge = Czip[iy]
		axleg1.scatter(x,y, marker=markers.next(), c='b', s=85)
		axleg2.scatter(x,y, marker='o', c=czip, edgecolors=cedge, s=85, linewidths=2)
		iy+=1
	ix += 1
lega = [axleg1, axleg2]
for legax in lega:
	ya = legax.get_yaxis()
	xa = legax.get_xaxis()
	ya.set_major_locator(MaxNLocator(integer=True))
	xa.set_major_locator(MaxNLocator(integer=True))
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
# =======
# ax.set_ylabel('Correlation [Normalized to 1]')
# ax.set_xlim([-0.01, 0.05])
# >>>>>>> 4edb61abe312fcaf2c92b9b7d1a92034c9afd333

mult1010 = float(np.amax(df['mult']))
#mult1010 = float(11025)
for i in xrange(Npts):
	marker = markerszip[i]
	label = str(df['sep'][i])+':'+str(df['sep2'][i]) if i>=Npts/2 else None
	ax2.scatter(np.abs(df['dT'][i]), (np.sqrt(df['mult']/mult1010)*df['peak'])[i]/peakmax, s=65, c=Cz[i], label=label, marker=marker)
#plt.legend()
if LEGEND:
	legend = ax2.legend(ncol=4, scatterpoints=1, frameon=False, fontsize = 'x-small')
	legend.get_frame().set_facecolor('none')

#ax2.grid()
ax2.set_xlabel('Time Offset [Sidereal Day]')
ax2.set_ylabel('Weighted Correlation')
#plt.tight_layout()
#f.subplots_adjust(hspace=0.02,bottom=0.2)

plt.show()


print "========= Statistics of sensitibity contribution =========="
df['Theta'] = df['peak']*np.sqrt(df['mult'])
df['Theta'] /= np.amax(df['Theta'])
def get_imp(df, Theta_min=0.0):
	dft = df.loc[df['Theta']>Theta_min]
	dfeq = dft.loc[dft['sep']==dft['sep2']]
	dfnq = dft.loc[dft['sep']!=dft['sep2']]
	
	totalsum = np.sum(dft['Theta'])
	eqsum = np.sum(dfeq['Theta'])
	totalsens = totalsum/np.sqrt(len(dft.index))
	eqsens = eqsum/np.sqrt(len(dfeq.index))
	improve = (totalsens-eqsens)/eqsens
	return improve
L = np.arange(0,1,0.01)
IL = [get_imp(df, L) for L in l]
import IPython; IPython.embed()

