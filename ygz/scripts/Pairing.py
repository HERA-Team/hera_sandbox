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
o.add_option('-y','--typ',dest='typ',default='simu',help='simu or real')
opts,args = o.parse_args(sys.argv[1:])
print opts, args

dt_fine = 0.0004971027374267578
if opts.typ == 'simu':
    dir1, dir2 = '/Users/yunfanzhang/local/simuDATA/64_UV/0_26/', '/Users/yunfanzhang/local/simuDATA/64_UV/0_38/'
    dt_file = dt_fine*8
elif opts.typ == 'real':
    dir1, dir2 = '/Users/yunfanzhang/local/DATA64/runDIR/', '/Users/yunfanzhang/local/DATA64/runDIR/'
    dt_file = dt_fine*13
sz = 200
sp = 1./sz
img = a.img.Img(200,res=0.5)   #400 by 400 image, i.e. 200 boxes, with 4 pixels per box,freq space kmax=100, dk=0.5
X,Y,Z = img.get_top(center=(200,200))
shape0 = X.shape
X,Y,Z = X.flatten(),Y.flatten(),Z.flatten()
ntop = n.array([X,Y,Z])
#aa = a.cal.get_aa('psa6622_v001',n.array([.15]))
list_freq = [.15]
dt = 0.001
#dt_fine = 43./3600/24


dist = 1.5                           #size of cells to store in dictionary.
corr_tol = 5000.                    #cutoff of minimum correlation
aa = a.cal.get_aa('psa6240_v003',n.array(list_freq))
src = a.fit.RadioFixedBody(0, aa.lat)
#src=a.fit.RadioSpecial("Sun")
fdict1, fdict2 = get_files.get_fdict(dir1), get_files.get_fdict(dir2)
T_files = fdict1.keys()
T_files.sort()
T_iter = T_files[:]              #This makes a new list, rather than just a pointer, which screws up the iteration with remove
for t in T_iter:
    aa.set_jultime(t)
    src.compute(aa)
    try: uvw = aa.gen_uvw(0,26,src=src).flatten()
    except(a.phs.PointingError):                                #test if below horizon
        print 'remove',t, len(T_iter), len(T_files)
        T_files.remove(t)
#print T_files
#times_coarse = n.arange(T_files[0],T_files[-1]+dt_file, dt)
#times_fine = n.arange(T_files[0],T_files[-1]+dt_file, dt_fine)
times_coarse = n.arange(T_files[0],T_files[-1], dt)
times_fine = n.arange(T_files[0],T_files[-1], dt_fine)

nants = len(aa)
bmp_list  = export_beam.beam_real(aa[0], ntop, shape0, 'x')

for ni in range(len(list_freq)):
    bmp = bmp_list[ni]
    freq, fbmamp = export_beam.beam_fourier(bmp, sp, 400)
    bm_intpl = export_beam.beam_interpol(freq,fbmamp,'cubic')
    print 'Time to initialize:', sys_time.clock(), 'seconds'
    print 'fbmampshape, midval', fbmamp.shape, fbmamp[200][200]

    d = select_pair.pair_coarse(aa, src,times_coarse,dist,False, 0.1, northsouth=False)  #coarsely determine crossings
    print 'Time after coarse selection:', sys_time.clock(), 'seconds'
    #pairs_sorted = select_pair.pair_sort(d,freq,fbmamp)        #sort crossings
    #clos_app = select_pair.get_closest(pairs_sorted)           #determine closest approach points
    clos_app = select_pair.alter_clos(d,bm_intpl)            #determine closest approach points
    print 'Found closest approach points after:', sys_time.clock(), 'seconds'
    pairs_final = select_pair.pair_fin(clos_app,dt_fine,aa,src,freq,fbmamp,multweight=True,noiseweight=True,ovlpweight=True,puv=False)
    print 'Total time:', sys_time.clock(), 'seconds'

    #write result to file and screen
    Oname = './P'+str(n.around(list_freq[ni],decimals=3))+'.out'
    Cname = './P'+str(n.around(list_freq[ni],decimals=3))+'.cue'
    f1 = open(Oname, 'w')
    f1.close()
    f1 = open(Cname, 'w')
    f1.close()
    print "Writting ourput files", Oname, Cname
    f1 = open(Oname, 'a')
    for j in n.arange(len(pairs_final)):
        #print pairs_final[j]
        f1.write(str(pairs_final[j])+'\n')
    f1.close()

    f1 = open(Cname, 'a')
    f1.write("bl1_bl2,T2-T1,Opp \n")
    for j in n.arange(len(pairs_final)):
        #print pairs_final[j]
        T1, T2 = float(pairs_final[j][2][1]), float(pairs_final[j][3][1])
        #fn1, fn2 = get_files.get_file(T1,dt_file, fdict1), get_files.get_file(T2,dt_file, fdict2)

        blstr = str(pairs_final[j][2][0][0])+'_'+str(pairs_final[j][2][0][1])+'_'+str(pairs_final[j][3][0][0])+'_'+str(pairs_final[j][3][0][1])
        lststr = str(-pairs_final[j][2][1]+pairs_final[j][3][1])     #time lag T2-T1

        OPP = pairs_final[j][1]
        #f1.write(blstr+','+lststr+' '+str(fn1)+' '+str(fn2)+' '+str(OPP)+'\n')
        f1.write(blstr+','+lststr+','+str(OPP)+'\n')
    f1.close()
    print 'compute extension:'
    print select_pair.pair_ext(pairs_final, aa, bm_intpl, ind=1)

    #call plotting routines
   # figname = './corr'+str(int(corr_tol))+str(n.around(list_freq[ni],decimals=3))+'.png'
   # print "Saving scatterplot", figname
   # plot_pair.plot_closapp(clos_app,corr_tol,figname)

    #plot sample approach points, puv in pair_fin must be True
    #pair_xampl = select_pair.test_sample(pairs_final)
    #plot_pair.plot_pair_xampl(pair_xampl)

