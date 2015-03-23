#! /usr/bin/env python
"""
Plots the output of 
uv_avg.py > nightly_avg.py
"""
from pylab import *
import sys,aipy as a
import numpy as n
import pickle
import optparse

ij2bl,bl2ij = a.miriad.ij2bl,a.miriad.bl2ij

o = optparse.OptionParser()
o.set_usage('plot_red_xtalk.py *.uv')
o.set_description(__doc__)
a.scripting.add_standard_options(o,cal=True)
o.add_option('--plotseps',default='1,0;2,0',help='The seperations to plot. [default=0,1;0,2]')
o.add_option('--verbose', action='store_true',
    help="Print a lot of stuff.")
o.add_option('--legend',action='store_true',
    help="""enable the plot legend. Its usually in the way, so its disabled by
default""")
o.add_option('--scale',type='float',
    help='Set the absolute value of the min/max scale. default = datamax, units=Jy')
o.add_option('-d',action='store_true',
    help="plot the delay spectrum")
opts, args = o.parse_args(sys.argv[1:])
if len(args)==2:
    print "plotting ",args[0],"/",args[1]
# PSA-128, JD2456240...
A_ = [49,41,47,19,29,28,34,51]
B_ = [10, 3,25,48,24,55,27,57]
C_ = [ 9,58, 1, 4,17,13,56,59]
D_ = [22,61,35,18, 5,32,30,23]
E_ = [20,63,42,37,40,14,54,50]
F_ = [43, 2,33, 6,52, 7,12,38]
G_ = [53,21,15,16,62,44, 0,26]
H_ = [31,45, 8,11,36,60,39,46]
ANTPOS_6240 = n.array([A_, B_, C_, D_,E_,F_,G_,H_])

if not opts.cal is None:
    aa = a.cal.get_aa(opts.cal,n.array([0.170]))
    ANTPOS = aa.ant_layout
else:
    ANTPOS = ANTPOS_6240

#for filename in args:
filename=args[0]
bls = {}
conj = {}
for ri in range(ANTPOS.shape[0]):
    for ci in range(ANTPOS.shape[1]):
        for rj in range(ANTPOS.shape[0]):
            for cj in range(ci,ANTPOS.shape[1]):
                if ri >= rj and ci == cj: continue # exclude repeat +/- listings of certain bls
                sep = '%d,%d' % (rj-ri, cj-ci)
                bls[sep] = bls.get(sep,[]) + [(ANTPOS[ri,ci],ANTPOS[rj,cj])]
for sep in bls.keys():
    if sep == '0,0' or len(bls[sep]) < 2: del(bls[sep])
for sep in bls:
    conj[sep] = [i>j for i,j in bls[sep]]

strbls = {}
conj_bl = {}
for sep in bls:
    strbls[sep] = []
    bl_list = []
    for (i,j),c in zip(bls[sep],conj[sep]):
        if c: i,j = j,i
        #valid = [0,1,2,3,16,17,18,19,8,9,10,11,24,25,26,27,4,5,6,7,20,21,22,23,12,13,14]
        #valid = [0,1,2,3,16,17,18,19,12,13,14,15,28,29,30,31]
        #if not i in valid or not j in valid: continue
        bl_list.append(ij2bl(i,j))
        strbls[sep].append('%d_%d' % (i,j))
        conj_bl[ij2bl(i,j)] = c
    bls[sep] = bl_list
    strbls[sep] = ','.join(strbls[sep])
    if opts.verbose: print sep, strbls[sep]
if opts.verbose: print "found %d kinds of baselines"%len(strbls)
print "loading:",filename
P = pickle.load(open(filename))
if len(args)==2:P2 = pickle.load(open(args[1]))

plotseps = opts.plotseps.split(';')
R,C = n.floor(n.sqrt(len(plotseps))),n.ceil(len(plotseps)/n.floor(n.sqrt(len(plotseps))))
print "finding global maximum for plotting homogeneity"
MAX = 0
for sep in plotseps:
    if opts.verbose:
        print "SEP = ",sep,
        print "\t BLS = ",strbls[sep]
    for bl in strbls[sep].split(','):
        if bl=='freqs':continue
        if bl.split('_')[0]==bl.split('_')[1]:continue
        if n.max(n.real(P[bl]))>MAX:MAX= n.max(n.real(P[bl]))
if not opts.scale is None:
    MAX = opts.scale
fig = figure()
for i,sep in enumerate(plotseps):
    #if not key.startswith('64'):continue
    #if key=='64_64':continue
    subplot(R,C,i+1)
    for bl in strbls[sep].split(','):
        if len(args)==2:
            myvis2 = n.ma.masked_where(P2[bl]==0,P2[bl])
            plot(P['freqs']*1e3,n.abs(P[bl])/n.abs(myvis2),'.',label=bl)
            ylim([-10,10])
            print bl,n.abs(n.mean(P[bl]/myvis2))
        else:
            if opts.d:
                W = P['counts'][bl]/n.float(n.max(P['counts'][bl]))
                weights = a.dsp.gen_window(P['counts'][bl].shape[0],'blackman-harris')
                _d,info = a.deconv.clean(n.fft.ifft(P[bl]*weights),n.fft.ifft(W),tol=1e-3)
                delays = n.fft.fftfreq(len(W),d=n.diff(P['freqs'])[0])
                plot(n.fft.fftshift(delays),n.fft.fftshift(n.abs(_d+info['res'])))
                ylim([0,MAX])
            else:
                plot(P['freqs']*1e3,P[bl],label=bl)
                ylim([-MAX,MAX])

    if opts.legend: legend()
    title(sep)
    if opts.d: xlabel('ns')
    else: xlabel('MHz')
    ylabel('real(<vis>)')
    if not opts.cal is None:
        i,j = map(int,bl.split('_'))
        uvw = aa.get_baseline(i,j)
        bl_length = n.sqrt(n.dot(uvw,uvw))
        vlines([-bl_length,bl_length],0,MAX)
mytitle = filename
if opts.d:mytitle += ' delay'
fig.canvas.set_window_title(mytitle)
show()

