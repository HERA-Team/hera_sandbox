__author__ = 'yunfanzhang'
import aipy as a, numpy as n
import select_pair, export_beam, plot_pair
import time as sys_time

f1 = open('./Pairing.out', 'w')
f1.close()
f1 = open('./Pairing.out', 'a')
sz = 200
d = 1./sz
img = a.img.Img(200,res=0.5)   #400 by 400 image, i.e. 200 boxes, with 4 pixels per box,freq space kmax=100, dk=0.5
X,Y,Z = img.get_top(center=(200,200))
shape0 = X.shape
X,Y,Z = X.flatten(),Y.flatten(),Z.flatten()
ntop = n.array([X,Y,Z])
aa = a.cal.get_aa('psa6622_v001',n.array([.15]))
src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.15, name='test')
#src=a.fit.RadioSpecial("Sun")
nants = 128
dt = 0.001
dt_fine = 43./3600/24
times_coarse = n.arange(2456240.3,2456240.4, dt)
times_fine = n.arange(2456240.3,2456240.4, dt_fine)
dist = 2.                           #size of cells to store in dictionary.
corr_tol = 5000.                    #cutoff of minimum correlation
bmp  = export_beam.beam_real(aa[0], ntop, shape0, 'x')
freq, fbmamp = export_beam.beam_fourier(bmp, d, 400)
print 'Time to initialize:', sys_time.clock(), 'seconds'

d = select_pair.pair_coarse(aa, src,times_coarse,dist,2.)  #coarsely determine crossings
print 'Time after coarse selection:', sys_time.clock(), 'seconds'
pairs_sorted = select_pair.pair_sort(d,freq,fbmamp)        #sort crossings
print 'Time after sorting:', sys_time.clock(), 'seconds'
clos_app = select_pair.get_closest(pairs_sorted)           #determine closest approach points
pairs_final = select_pair.pair_fin(clos_app,5*dt,aa,src,freq,fbmamp,corr_tol)  #output final sorted pairs
for j in n.arange(len(pairs_final)):
    print pairs_final[j]
print 'Total time:', sys_time.clock(), 'seconds'


#call plotting routines
figname = './corr'+str(int(corr_tol))+'.png'
plot_pair.plot_closapp(clos_app,corr_tol,figname)
pair_xampl = select_pair.test_sample(pairs_final,dt,aa, src,freq,fbmamp,3000.)
plot_pair.plot_pair_xampl(pair_xampl)

f1.close()