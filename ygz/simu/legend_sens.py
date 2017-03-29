import numpy as n, matplotlib.pyplot as p
from mpl_toolkits.axes_grid.axislines import SubplotZero
from random import random
# fig = p.figure(1)
# ax = SubplotZero(fig, 111)
# fig.add_subplot(ax)
# for direction in ["xzero", "yzero"]:
#     ax.axis[direction].set_axisline_style("-|>")
#     ax.axis[direction].set_visible(True)
# for direction in ["left", "right", "bottom", "top"]:
#     ax.axis[direction].set_visible(False)
def gen_color(l=1):
	colors = []
	for i in range(l): colors.append((random(),random(),random()))
	return n.array(colors)
mark = {'1':'+','-':'x','2':'d','3':'^','4':'s'}
colors_temp = gen_color(l=10)
def set_color(sep):
	a,b,c,d = int(sep[0]),int(sep[1]),int(sep[-2]),int(sep[-1])
	if sep[3] == '-': 
		dc = d+b
	else:
		dc = -d+b
	return n.sqrt(colors_temp[dc]*random())
arr = n.genfromtxt('something.csv', dtype=None,delimiter=' ',names=True)
dt = arr['dT']
##############
dt = n.abs(dt)
##############
corr = arr['peak']
sep = []
for s1, s2 in zip(arr['sep'],arr['sep2']): 
	sep.append(str(s1)+':'+str(s2))
mult = arr['mult']
colors = gen_color(len(sep))

fig = p.figure()
ax1 = fig.add_subplot(211)
for a,b,c,d in zip(dt,corr,colors,sep):
	marker = mark[d[0]]
	ax1.scatter(n.abs(a), b,color=c,label=d,marker=marker,s=40)
ax1.grid()
ax1.xaxis.set_ticks(n.arange(0, 0.225, 0.025))
ax1.set_ylabel('Raw Correlation [Normalized to 1]')
ax1.set_xlim([-0.01,0.2])
#ax1.set_title("Raw Peak Heights")
p.setp(ax1.get_xticklabels(), visible=False)
p.legend(scatterpoints=1,ncol=5,fontsize=10,loc=1,frameon=False)

ax2 = fig.add_subplot(212,sharex=ax1)
mult = mult/float(3136)   #sep10,10
for a,b,b2,c,d in zip(dt,corr,mult,colors,sep):
	marker = mark[d[0]]
	#ax2.scatter(n.abs(a), b*b2,color=set_color(d),label=d,marker=marker,s=40)
	ax2.scatter(n.abs(a), b*b2,color=c,label=d,marker=marker,s=40)
ax2.xaxis.set_ticks(n.arange(0, 0.225, 0.025))
ax2.grid()
ax2.set_xlim([-0.01,0.2])
#ax2.set_title("Corrected for Multiplicities")
ax2.set_xlabel('Time Delay [Sidereal Day]')
ax2.set_ylabel('Correlation [Weighted by Multiplicities]')
p.show()