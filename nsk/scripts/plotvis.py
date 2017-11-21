"""
plotvis.py
----------

casa plotting script
"""
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import argparse
import sys
import os
import subprocess
import shutil
import glob

a = argparse.ArgumentParser()
a.add_argument('--msin', default=None, type=str, help='path to CASA measurement set.', required=True)
a.add_argument('--figfile', default=None, type=str, help="output filename.", required=True)
a.add_argument('--xaxis', default='')
a.add_argument('--yaxis', default='')
a.add_argument('--timebin', default=0, type=int, help='time bin to plot'
a.add_argument('--polaxis', default=0, type=int, help='polarization axis to plot')
a.add_argument('--ant1', default=0, type=int, help='antenna1 to plot')
a.add_argument('--ant2', default=1, type=int, help='antenna2 to plot')
args = a.parse_args()

# load data
ms.open(args.msin)
data = ms.getdata(['data', 'time', 'antenna1', 'antenna2', 'axis_info'], ifraxis=True)

# parse ants
ant1 = d['antenna1']
ant2 = d['antenna2']
s = np.where((ant1==args.ant1)&(ant2==args.ant2))[0]

# make ydata
if args.yaxis == 'phs':
    y = np.angle(d['data'][args.polaxis, :, s, args.timebin])

elif args.yaxis == 'amp':
    y = np.abs(d['data'][args.polaxis, :, s, args.timebin]).T

# plot
fig, ax = plt.subplots(figsize=(10,8))
ax.grid(True)
ax.plot(y)
ax.set_xlabel('channel', fontsize=14)
ax.set_ylabel(args.yaxis, fontsize=14)
ax.set_title("({},{}) {}, {}".format(args.ant1, args.ant2, args.yaxis, args.msin))
if args.figfile is None:
    figfile = "vis_{}_{}_{}.png".format(args.yaxis, args.ant1, args.ant2)
else:
    figfile = args.figfile
fig.savefig(figfile)
ms.close()






