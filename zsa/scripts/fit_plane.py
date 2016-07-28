#! /usr/bin/env python
import numpy as np, aipy as a, pylab as p, time
import capo.omni as omni, capo.miriad as miriad, capo.zsa as zsa
import scipy
from mpl_toolkits.mplot3d import Axes3D
import sys,os

def flatten_reds(reds):
    freds = []
    for r in reds:
        freds += r
    return freds

def fit_line_to_phase(phs, fqs, valid):
    '''Fits a line (slope and offset) to the phase components. 
       Returns m in inverse radians (i.e. need to multiply by 2pi) and offset in radians'''
    fqs = fqs.compress(valid)
    dly = phs.compress(valid)
    B = np.zeros((fqs.size,1)); B[:,0] = dly
    A = np.zeros((fqs.size,2)); A[:,0] = fqs*2*np.pi; A[:,1] = 1
    dt,off = np.linalg.lstsq(A,B)[0].flatten()
    return dt,off

def fit_plane(d, x, y, verbose=False):
    '''Generic function to fit a plane to d. Takes data (d) and the corresponding 
       x and y values. 
       By default, return slope in x,y direction and offset in z (units of d).
       If verbose include residual error, rank of A and singular values. 
       See linalg.lstsq'''
    A = np.c_[x, y, np.ones(len(x))]
    C,res,rank,sing = scipy.linalg.lstsq(A, d)
    if verbose:
        return C, res, rank, sing
    else:
        return C

def vector_fit_line_to_phase(phs, fqs, valid):
    '''Vectorize the fit line function with a for loop.'''
    dts, offs = [], []
    for ph in phs:
        dt, off = fit_line_to_phase(ph, fqs, valid)
        dts.append(dt); offs.append(off)
    return np.array(dts), np.array(offs)
        

filename = sys.argv[1:] #uv file with first cal solutions applied to it.
aa = a.cal.get_aa('hsa7458_v000_HH_delaytest', np.array([.150]))
fqs = np.linspace(.1,.2,1024)
valid = np.ones_like(fqs)
#only fit phase to within this frequency range.
valid = np.logical_and(valid, np.logical_and(fqs>.11,fqs<.19))

info = omni.aa_to_info(aa)
reds = flatten_reds(info.get_reds())
antstr = zsa.list2str(reds)
#wh = aa.antpos_ideal[:,2]!=-1
#ants = np.arange(len(aa.ants)) 
#ants = ants[wh]

integration = 0 # testing for a single integration
pol='xx'
times, data, flags = miriad.read_files(filename, antstr, pol)


#bls={}
bs = []
xs = []
ys = []
ds = []
os = []
for bl in data.keys():
    x,y = aa.get_baseline(*bl)[:2]
#    if x > 0.0:
    phs = np.unwrap(np.angle(data[bl][pol]))[integration:integration+1,:] # try plane fitting for a single time.
#    phs = np.unwrap(np.angle(data[bl][pol]))
#    elif x == 0 and y<0:
#        print 'conjugating'
#        phs = np.conj(np.unwrap(np.angle(data[bl][pol][time])))
#        bl = (bl[1],bl[0])
#        x,y = aa.get_baseline(*bl)[:2]
#    elif x < 0.0: 
#        print 'conjugating'
#        phs = np.conj(np.unwrap(np.angle(data[bl][pol][time])))
#        bl = (bl[1],bl[0])
#        x,y = aa.get_baseline(*bl)[:2]
#    dlys, offsets = vector_fit_line(phs, fqs, valid)
    dlys, offsets = vector_fit_line_to_phase(phs, fqs, valid)
    if False:
        fq_plot = np.c_[ [fqs for i in range(phs.shape[0])] ].T
    #    import IPython; IPython.embed()
        p.subplot(211)
        p.plot(fqs, 2*np.pi*fq_plot*dlys + offsets, lw=4)
        p.plot(fqs, 2*np.pi*fq_plot*dlys, lw=4)
        p.plot(fqs,phs.flatten(), 'k')
        p.plot(fqs,np.angle(data[bl][pol][integration]))
        p.subplot(212)
        p.semilogy(fqs,np.abs(data[bl][pol][integration]))
        p.show()
    bs.append(bl)
    xs.append(x)
    ys.append(y)
    ds.append(dlys)
    os.append(offsets)
    #bls[bl] = np.array((x,y,dlys))


#Now fit it the data points, xs,ys,ds
i,k = np.meshgrid(np.arange(-150,150,5), np.arange(-150,150,5))
ii = i.flatten()
kk = k.flatten()

#ds = np.abs(ds)


A = np.c_[xs, ys, np.ones(len(xs))]
C,residuals,rank,singulars = scipy.linalg.lstsq(A, ds)
O,oresiduals,orank,osingulars = scipy.linalg.lstsq(A, os)

zs = C[0]*i + C[1]*k + C[2]
ozs = O[0]*i + O[1]*k + O[2]

fig = p.figure()
ax1 = fig.add_subplot(121, projection='3d')
ax1.plot_surface(i,k,zs,rstride=5,cstride=5,alpha=.2)
ax1.scatter(xs,ys,ds,c='r',s=20)
ax1.set_xlabel('x-length (ns)')
ax1.set_ylabel('y-length (ns)')
ax1.set_zlabel('dly')

ax2 = fig.add_subplot(122, projection='3d')
ax2.plot_surface(i,k,ozs,rstride=5,cstride=5,alpha=.2)
ax2.scatter(xs,ys,os,c='r',s=20)
ax2.set_xlabel('x-length (ns)')
ax2.set_ylabel('y-length (ns)')
ax2.set_zlabel('offsets')

p.show()

import IPython; IPython.embed()


def save_gains(s,f,name='fcgains',verbose=False):
    s2 = {}
    for k,i in s.iteritems():
        try:
            len(i)
            s2[str(k)] = omni.get_phase(f,i,offset=True)
            s2['d'+str(k)] = i[0]
            if verbose:
                print 'ant=%d dly=%f , off=%f '%(k,i[0],i[1])
        except(TypeError):
            s2[str(k)] = omni.get_phase(f,i)
            s2['d'+str(k)] = i
            if verbose:
                print 'ant=%d dly=%f  '%(k,i)
    import sys
    cmd = sys.argv
    s2['cmd'] = ' '.join(cmd)
    np.savez('%s.npz'%name,**s2)

update_delay = {}
#update_offset = {}
cx,cy = C[:2]
ox,oy = O[:2]
for i in info.subsetant:
    bl = aa.get_baseline(80,i)[:2]
    update_delay[i] = np.sum(bl*[cx,cy])#, np.sum(bl*[ox,oy]))
#    update_offset[i] = np.sum(bl*[ox,oy])

save_gains(update_delay,fqs,name='updates')

#fcalnpz = filename[0][:-1] + '.npz'
#npz = np.load(fcalnpz)
#for key in npz.keys():
#    if key,startswith('d'):
#        update_delay[int(key[1:])] += npz[key]
    





