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
ntimes = len(_T)*len(args)
del(uv)
assert(nants == len(antpos.keys())) #check the cal file and uv files are going to cooperate

times,data,flags = capo.arp.get_dict_of_uv_data(args,'cross',opts.pol)
vis_stor = np.zeros((nants,nchan,ntimes),dtype='complex128')
flg_stor = np.zeros_like(vis_stor)

for i in range(nants):
    for j in range(nants):
        if i==j: continue #neglect autos
        try: 
            vis_stor[i,:,:]+=data[(i,j)][opts.pol].T
            #flg_stor[i,:,:]+=np.logical_not(flags[(i,j)][opts.pol]).astype('complex128')
            flg_stor[i,:,:] += np.ones((nchan,ntimes),dtype='complex128')
        except KeyError:
            if i<j: print 'KeyError on (%i,%i)'%(i,j) #this should not happen
            continue

#average all abs visibilities |Vij| over j per i
mean_stor = np.absolute(vis_stor)/np.absolute(flg_stor) 

#parse options
if not opts.time is None: tlo,thi = map(int,opts.time.split('_'))
else: tlo,thi = 0,ntimes-1
if not opts.chan is None: clo,chi = map(int,opts.chan.split('_'))
else: clo,chi = 0,nchan-1

#Plot channel/time average on antenna positions
print 'Plotting'
_x,_y,_avg = [],[],[]
for i in range(nants):
    x,y = antpos[i]['top_x'],antpos[i]['top_y']
    avg = np.nanmean(mean_stor[i,clo:chi,tlo:thi])
    _x.append(x)
    _y.append(y)
    _avg.append(avg)
plt.scatter(_x[:112],_y[:112],s=40,c=_avg[:112]) #XXX psa128 specific - add option
plt.colorbar()
plt.show()

    


    
        

