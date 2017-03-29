import numpy as n, matplotlib.pyplot as p
from mpl_toolkits.axes_grid.axislines import SubplotZero

# fig = p.figure(1)
# ax = SubplotZero(fig, 111)
# fig.add_subplot(ax)
# for direction in ["xzero", "yzero"]:
#     ax.axis[direction].set_axisline_style("-|>")
#     ax.axis[direction].set_visible(True)
# for direction in ["left", "right", "bottom", "top"]:
#     ax.axis[direction].set_visible(False)


arr = n.genfromtxt('corr_res.csv', dtype=None,delimiter=' ',names=True)
dt = arr['dT']
corr = arr['peak']
sep = str(arr['sep'])+':'+str(arr['sep2'])
#import IPython; IPython.embed()
mult = arr['mult'].astype(n.float)
fig, ax = p.subplots()
ax = fig.add_subplot(211)
ax.scatter(n.abs(dt), corr)
for label, x, y in zip(sep, dt, corr):
	if label in ['31:32']:
		p.annotate(
	        label, 
	        xy = (x, y), xytext = (10, 25),
	        textcoords = 'offset points', ha = 'left', va = 'bottom',
	        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
	        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	elif label in ['21:22']:
		p.annotate(
	        label, 
	        xy = (x, y), xytext = (20, 20),
	        textcoords = 'offset points', ha = 'left', va = 'bottom',
	        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
	        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	elif label in ['20:21']:
		p.annotate(
	        label, 
	        xy = (x, y), xytext = (30, -15),
	        textcoords = 'offset points', ha = 'left', va = 'bottom',
	        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
	        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	else:
	    p.annotate(
	        label, 
	        xy = (x, y), xytext = (-20, 20),
	        textcoords = 'offset points', ha = 'right', va = 'bottom',
	        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
	        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
# for i, txt in enumerate(sep):
#     ax.annotate(txt, (dt[i],corr[i]))
ax.grid()
ax.set_xlabel('Time Delay [Sidereal Day]')
ax.set_ylabel('Correlation [Normalized to 1]')

ax = fig.add_subplot(212)
ax.scatter(n.abs(dt), corr)
for label, x, y,yf in zip(sep, dt, corr, mult):
	y = y*yf
	if label in ['31:32']:
		p.annotate(
	        label, 
	        xy = (x, y), xytext = (10, 25),
	        textcoords = 'offset points', ha = 'left', va = 'bottom',
	        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
	        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	elif label in ['21:22']:
		p.annotate(
	        label, 
	        xy = (x, y), xytext = (20, 20),
	        textcoords = 'offset points', ha = 'left', va = 'bottom',
	        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
	        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	elif label in ['20:21']:
		p.annotate(
	        label, 
	        xy = (x, y), xytext = (30, -15),
	        textcoords = 'offset points', ha = 'left', va = 'bottom',
	        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
	        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
	else:
	    p.annotate(
	        label, 
	        xy = (x, y), xytext = (-20, 20),
	        textcoords = 'offset points', ha = 'right', va = 'bottom',
	        bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
	        arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
# for i, txt in enumerate(sep):
#     ax.annotate(txt, (dt[i],corr[i]))
ax.grid()
ax.set_xlabel('Time Delay [Sidereal Day]')
ax.set_ylabel('Correlation [Normalized to 1]')
p.show()