from matplotlib import use
use('Agg')
from pylab import *
import numpy as n
import  sys,os,optparse

o = optparse.OptionParser()
o.set_usage('plot_wedge_singlebl.py [options] *.uv')
o.set_description(__doc__)
o.add_option('--bl',default='11,12',
    help='plot this baseline [default 0,1] units= MWA antnum units (or whatever is the standard for the original file')
opts, args = o.parse_args(sys.argv[1:])
mybl = map(int,opts.bl.split(','))
X,Y = [],[]
pols = ['xx','yy']
for filename in args:
    print "loading",filename
    D = n.load(filename)
    try:
        ant1 = D['ant1']
        ant2 = D['ant2']
        P = D['P']
        P2 = D['P2']
        C = D['C']
        bl_lengths = D['bl_lengths']
        delays = D['delays']
    except(KeyError):
        print "weird file"
        continue
    pspec = P*n.conj(P) - P2
    pspec[C>0] /= C[C>0]
    mypspec_i = n.argwhere(n.logical_and(ant1==mybl[0],ant2==mybl[1])).squeeze()
    X.append(n.fft.fftshift(pspec[mypspec_i].squeeze()[:,0]))
    Y.append(n.fft.fftshift(pspec[mypspec_i].squeeze()[:,1]))
print len(X)
figure()
subplot(121)
imshow(n.log10(n.abs(X)),aspect='auto')
subplot(122)
imshow(n.log10(n.abs(Y)),aspect='auto')
savefig('delay_spectrum_singlebl_{mybl}.png'.format(mybl=opts.bl))
