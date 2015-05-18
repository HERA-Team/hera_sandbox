__author__ = 'yunfanzhang'
import aipy as a, numpy as n
import select_pair, export_beam, plot_pair
import time as sys_time


sz = 200
sp = 1./sz
img = a.img.Img(200,res=0.5)   #400 by 400 image, i.e. 200 boxes, with 4 pixels per box,freq space kmax=100, dk=0.5
X,Y,Z = img.get_top(center=(200,200))
shape0 = X.shape
X,Y,Z = X.flatten(),Y.flatten(),Z.flatten()
ntop = n.array([X,Y,Z])
#aa = a.cal.get_aa('psa6622_v001',n.array([.15]))
list_freq = [.15,.18]
dt = 0.001
#dt_fine = 43./3600/24
dt_fine = 0.0004971027374267578
times_coarse = n.arange(2456249.20169,2456249.50791, dt)
times_fine = n.arange(2456249.20169,2456249.50791, dt_fine)
dist = 1.2                           #size of cells to store in dictionary.
corr_tol = 5000.                    #cutoff of minimum correlation
aa = a.cal.get_aa('psa6240_v003',n.array(list_freq))
src = a.fit.RadioFixedBody(0, aa.lat, janskies=0., mfreq=.18, name='test')
#src=a.fit.RadioSpecial("Sun")
nants = len(aa)
bmp_list  = export_beam.beam_real(aa[0], ntop, shape0, 'x')

for ni in range(len(list_freq)):
    bmp = bmp_list[ni]
    freq, fbmamp = export_beam.beam_fourier(bmp, sp, 400)
    bm_intpl = export_beam.beam_interpol(freq,fbmamp,'cubic')
    print 'Time to initialize:', sys_time.clock(), 'seconds'

    d = select_pair.pair_coarse(aa, src,times_coarse,dist)  #coarsely determine crossings
    print 'Time after coarse selection:', sys_time.clock(), 'seconds'
    #pairs_sorted = select_pair.pair_sort(d,freq,fbmamp)        #sort crossings
    #clos_app = select_pair.get_closest(pairs_sorted)           #determine closest approach points
    clos_app = select_pair.alter_clos(d,bm_intpl)            #determine closest approach points
    print 'Found closest approach points after:', sys_time.clock(), 'seconds'
    pairs_final = select_pair.pair_fin(clos_app,dt_fine,aa,src,freq, fbmamp,multweight=True,noiseweight=False,ovlpweight=False,cutoff=5000.)
    print 'Total time:', sys_time.clock(), 'seconds'

    #write result to file and screen
    Oname = './P'+str(n.around(list_freq[ni],decimals=3))+'.out'
    f1 = open(Oname, 'w')
    f1.close()
    f1 = open(Oname, 'a')
    for j in n.arange(len(pairs_final)):
        #print pairs_final[j]
        f1.write(str(pairs_final[j])+'\n')
    f1.close()

    #call plotting routines
    figname = './corr'+str(int(corr_tol))+str(n.around(list_freq[ni],decimals=3))+'.png'
    plot_pair.plot_closapp(clos_app,corr_tol,figname)
    #pair_xampl = select_pair.test_sample(pairs_final)
    #plot_pair.plot_pair_xampl(pair_xampl)

