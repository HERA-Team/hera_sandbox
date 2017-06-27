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
FILE = 'corr_res_rephs.csv'
#FILE = 'HERA_350_core_pm300.csv'
LEGEND = False

#markers = matplotlib.markers.MarkerStyle.markers
markers = cycle(['o', 'v', '^', '<', '>', '8', 's', 'p', '*', 'h', 'd', 'H', 'D'])
print markers
df = pd.DataFrame.from_csv(FILE)
Npts = df.shape[0]
#Npts = 1000
C = np.random.random(Npts)
#f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
# f.set_figheight(10);  f.set_figwidth(12)
X = np.arange(1, 4, 1).astype(np.int)
Y = np.arange(0, 4, 1).astype(np.int)
#X, Y = np.meshgrid(X, Y)
f = plt.figure()
axleg1 = plt.subplot2grid((3, 2), (0, 0))
axleg2 = plt.subplot2grid((3, 2), (0, 1))
axleg1.set_title("E-W separation")
axleg2.set_title("E-W separation")
axleg1.set_ylabel("N-S separation")
axleg1.set_xlabel("E-W separation")
ax = plt.subplot2grid((3, 2), (1, 0), colspan=2)
ax2 = plt.subplot2grid((3, 2), (2, 0), colspan=2, sharex=ax)

def get_marker(ix, iy):
	if ix< 2:
		if iy == 0 and ix == 1:
			marker = (iy+3, ix, 90)
		else:
			marker = (iy+3, ix, 0)
	else:
		angle = 360/(iy+3)/2.
		marker = (iy+3, 0, angle)
	return marker
Cx = ['r', 'g', 'b', 'y']
peakmax = np.amax(df['peak'])
colors = mcolors.cnames.keys()
#Cy = [np.random.choice(colors) for i in xrange(Y.size)]
Cy = Cx
ix = 0
for x in X:
	iy = 0
	cx = Cx[ix]
	for y in Y:
		cedge = Cy[iy]
		# axleg1.scatter(x,y, marker=markers.next(), c='b', s=100)
		marker = get_marker(ix,iy)
		axleg1.scatter(x,y, marker=marker, c='b', s=100)
		axleg2.scatter(x,y, marker='o', c=cx, edgecolors=cedge, s=100, linewidths=2)
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
#import IPython; IPython.embed()



for i in xrange(Npts):
	label = str(df['sep'][i])+':'+str(df['sep2'][i])
	ix1, iy1 = int(str(df['sep'][i])[0])-1, int(str(df['sep'][i])[1])
	ix, iy = int(str(df['sep2'][i])[0])-1, int(str(df['sep2'][i])[1])
	marker = get_marker(ix1, iy1)

	ax.scatter(np.abs(df['dT'][i]), df['peak'][i]/peakmax, 
				s=50, 
				marker=marker, 
				c=Cx[ix], 
				edgecolors=Cy[iy], 
				linewidths=2)
plt.setp(ax.get_xticklabels(), visible=False)

ax.grid()
#ax.set_xlabel('Time Delay [Sidereal Day]')
ax.set_ylabel('Correlation')
ax.set_xlim([-0.01, 0.16])
# =======
# ax.set_ylabel('Correlation [Normalized to 1]')
# ax.set_xlim([-0.01, 0.05])
# >>>>>>> 4edb61abe312fcaf2c92b9b7d1a92034c9afd333

mult1010 = float(np.amax(df['mult']))
#mult1010 = float(11025)
for i in xrange(Npts):
	ix1, iy1 = int(str(df['sep'][i])[0])-1, int(str(df['sep'][i])[1])
	ix, iy = int(str(df['sep2'][i])[0])-1, int(str(df['sep2'][i])[1])
	marker = get_marker(ix1, iy1)
	ax2.scatter(np.abs(df['dT'][i]), (np.sqrt(df['mult']/mult1010)*df['peak'])[i]/peakmax, 
				s=50, 
				marker=marker, 
				c=Cx[ix], 
				edgecolors=Cy[iy],
				linewidths=2)

ax2.grid()
ax2.set_xlabel('Time Offset [Sidereal Day]')
ax2.set_ylabel('Weighted Correlation')
#plt.tight_layout()
#f.subplots_adjust(hspace=0.02,bottom=0.2)
#ax.set_xlims([-0.02, 1.16])
plt.show()


# print "========= Statistics of sensitibity contribution =========="
# df['Theta'] = df['peak']*np.sqrt(df['mult'])
# df['Theta'] /= np.amax(df['Theta'])
# def get_imp(df, Theta_min=0.0):
# 	dft = df.loc[df['Theta']>Theta_min]
# 	dfeq = dft.loc[dft['sep']==dft['sep2']]
# 	dfnq = dft.loc[dft['sep']!=dft['sep2']]
	
# 	totalsum = np.sum(dft['Theta'])
# 	eqsum = np.sum(dfeq['Theta'])
# 	totalsens = totalsum/np.sqrt(len(dft.index))
# 	eqsens = eqsum/np.sqrt(len(dfeq.index))
# 	improve = (totalsens-eqsens)/eqsens
# 	return improve
# L = np.arange(0,1,0.01)
# IL = [get_imp(df, L) for L in l]
# import IPython; IPython.embed()

