__author__ = 'yunfanzhang'

import aipy as a, numpy as n
import select_pair, export_beam, plot_pair, get_files
import time as sys_time
import optparse, sys

o = optparse.OptionParser()
#o.add_option('-t', '--lst', dest='lst', default=-1, help='Choose the time point of data')
#o.add_option('-d', '--dft', dest='dst', default=43./3600/24)
o.add_option('-s','--sys', dest='sys', default='MacPro')
#o.add_option('-d','--dir',dest='dir',default='/Users/yunfanzhang/local/simuDATA/64_UV')
opts,args = o.parse_args(sys.argv[1:])
print opts, args

if opts.sys == 'MacPro':
    dir1, dir2 = '/Users/yunfanzhang/local/simuDATA/64_UV/0_26/', '/Users/yunfanzhang/local/simuDATA/64_UV/0_38/'
sz = 200
sp = 1./sz
img = a.img.Img(200,res=0.5)   #400 by 400 image, i.e. 200 boxes, with 4 pixels per box,freq space kmax=100, dk=0.5
X,Y,Z = img.get_top(center=(200,200))
shape0 = X.shape
X,Y,Z = X.flatten(),Y.flatten(),Z.flatten()
ntop = n.array([X,Y,Z])
aa = a.cal.get_aa('psa6240_v003',n.array(0.15))
ant = aa[0]
Omp = export_beam.OmP(ant,ntop,'x')
Ompp = export_beam.OmP(ant,ntop,'x',sq=True)

dl = abs(ntop[0][2]-ntop[0][1])


print Omp,Ompp, Omp*Omp/Ompp