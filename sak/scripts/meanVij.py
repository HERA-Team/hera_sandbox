#! /usr/bin/env python
import sys, numpy as np, aipy, optparse, capo
from matplotlib import pyplot as plt

o = optparse.OptionParser()
o.set_usage('meanVij.py [options] *.uvcRRE')
o.set_description(__doc__)
aipy.scripting.add_standard_options(o, cal=True, pol=True)
o.add_option('-c', '--chan', dest='chan', default=None,\
                     help='Channel range in form "lo_hi" to average over.')
o.add_option('-t', '--time', dest='time', default=None,\
                     help='Time range in form "lo_hi" (index) to average over.')
o.add_option('--ba', dest='badants', default='', help='comma separated list of bad antennae')
o.add_option('--autos',dest='autos',default=False,action='store_true',\
             help ='Plot autocorrelations instead of crosses.')
o.add_option('--antpos',dest='pos',default=False,action='store_true',\
             help ='Plot mean(Vij) as color on antenna positions')

opts,args = o.parse_args(sys.argv[1:])

#get antenna positions
print 'reading, %s'%opts.cal
exec("import {calfile} as cal".format(calfile=opts.cal))
antpos = cal.prms['antpos']

#get initial info from first file
uv = aipy.miriad.UV(args[0])
nants,nchan = uv['nants'],uv['nchan']
uv.select('antennae',0,1,include=True) #XXX assumes ants 0 and 1 are in uv file
_T=[]
for p,d in uv.all(): 
    _,_t,_ = p
    _T.append(_t)
ntimes = len(_T)
del(uv)
assert(nants == len(antpos.keys())) #check the cal file and uv files are going to cooperate

vis_stor = np.zeros((nants,nchan,ntimes),dtype='complex128')
flg_stor = np.zeros_like(vis_stor)

for uv in args:
    print '    Reading %s...'%uv
    if not opts.autos: times,data,flags = capo.arp.get_dict_of_uv_data([uv],'cross',opts.pol) 
    else: times,data,flags = capo.arp.get_dict_of_uv_data([uv],'all',opts.pol)
    for i in range(nants):
        for j in range(nants):
            if not opts.autos:
                if i==j: continue #neglect autos
            else:
                if i!=j: continue #neglect crosses
            try: 
                vis_stor[i,:,:]+=np.absolute(data[(i,j)][opts.pol].T)
                #flg_stor[i,:,:]+=np.logical_not(flags[(i,j)][opts.pol]).astype('complex128')
                flg_stor[i,:,:] += np.ones((nchan,ntimes),dtype='complex128')
            except KeyError:
                if i<j: print 'KeyError on (%i,%i)'%(i,j) #this should not happen
                continue
    #del(uv)
#average all abs visibilities |Vij| over j per i
mean_stor = np.absolute(vis_stor)/np.absolute(flg_stor) 

#parse options
if not opts.time is None: tlo,thi = map(int,opts.time.split('_'))
else: tlo,thi = 0,ntimes-1
if not opts.chan is None: clo,chi = map(int,opts.chan.split('_'))
else: clo,chi = 0,nchan-1
if not len(opts.badants)==0: badants=map(int,opts.badants.split(','))
else: badants=None

_x,_y,_avg = [],[],[]
for i in range(nants):
    x,y = antpos[i]['top_x'],antpos[i]['top_y']
    avg = np.nanmean(mean_stor[i,clo:chi,tlo:thi])
    _x.append(x)
    _y.append(y)
    _avg.append(avg)

#Plot channel/time average on antenna positions
print 'Plotting'
## the :112 is for the psa128 grid.
ga = range(0,112)

plt.scatter(_x[:ga[-1]],_y[:ga[-1]],s=40,c=np.log10(np.array(_avg[:ga[-1]]))) 
cbar = plt.colorbar()
cbar.ax.set_ylabel(r'$\left\langle |V_{ij}| \right\rangle_{j,t,\nu}$')

#there has to be a simpler way to do this
if not badants is None:
    for i in range(nants):
        if i in badants:
            plt.plot(_x[i],_y[i],'kx',ms=20)
plt.xlabel('E-W')
plt.ylabel('N-S')
if opts.pos: plt.show()
plt.close()

_avg = np.array(_avg)
mean_avg = np.nanmean(_avg) #I know, terrible variable name, but that's what it is!
std_avg = np.nanstd(_avg)

out = []

for i in ga:
    plt.plot(i,_avg[i],'bo',ms=8)
    if _avg[i]<=mean_avg-std_avg: out.append(i)
    if not badants is None:
        if i in badants:
            plt.plot(i,_avg[i],'kx',ms=10)
plt.axhline(mean_avg,color='k')

plt.fill_between(ga,mean_avg-std_avg,mean_avg+std_avg,color='b',alpha=0.5)
plt.fill_between(ga,mean_avg-2*std_avg,mean_avg+2*std_avg,color='b',alpha=0.2)

plt.grid()
plt.xlim(ga[0]-0.5,ga[-1]+0.5)
plt.ylabel(r'$\left\langle |V_{ij}| \right\rangle_{j,t,\nu}$',size=15)
plt.xlabel('Antenna number')
plt.show()
plt.ylim(0,0.045)
#plt.savefig('%s_meanVij.png'%uv)
print out
