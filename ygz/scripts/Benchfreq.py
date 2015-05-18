__author__ = 'yunfanzhang'
import aipy as a, numpy as n
import select_pair, export_beam, plot_pair
import time as sys_time


sz = 200
d = 1./sz
img = a.img.Img(200,res=0.5)   #400 by 400 image, i.e. 200 boxes, with 4 pixels per box,freq space kmax=100, dk=0.5
X,Y,Z = img.get_top(center=(200,200))
shape0 = X.shape
X,Y,Z = X.flatten(),Y.flatten(),Z.flatten()
ntop = n.array([X,Y,Z])
list_freq = [.15]
dt = 0.001
dt_fine = 0.0004971027374267578
times_coarse = n.arange(2456249.20169,2456249.50791, dt)
times_fine = n.arange(2456249.20169,2456249.50791, dt_fine)
dist = 0.2                           #size of cells to store in dictionary.
corr_tol = 5000.                    #cutoff of minimum correlation
#aa = a.cal.get_aa('psa6240_v003',n.array(list_freq))
aa = a.cal.get_aa('psa898_v003',n.array(list_freq))
src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.18, name='test')
nants = len(aa)
bmp_list  = export_beam.beam_real(aa[0], ntop, shape0, 'x')


ni = 0
bmp = bmp_list[ni]
freq, fbmamp = export_beam.beam_fourier(bmp, d, 400)
d = select_pair.pair_coarse(aa, src,times_coarse,dist)  #coarsely determine crossings
clos_app = select_pair.alter_clos(d,freq,fbmamp)            #determine closest approach points
pairs_final = select_pair.pair_fin(clos_app,dt_fine,aa,src,freq,fbmamp,False,False,False,5000.)  #output final sorted pairs

