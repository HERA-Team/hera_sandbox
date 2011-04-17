#mixes AIPY and CASA to grid down uv samples using central freq to compute 
#uvw
# DCJ 2011
import aipy as a
#!/home/djacobs/src/casa/bin/casapy --logfile test.log -c
# NB: this script works as an executable on linux kernels 
#                        V>2.6.27.p [1]
#   otherwise run with
#     casapy --logfile test.log --nologger -c <scriptname> \
#                  <options> <args>
#
# [1] http://www.in-ulm.de/~mascheck/various/shebang/#interpreter-script
import sys,os,optparse,inspect
from glob import glob
####
#PAPER imports
from casachecksum import checkmscolumn
from IPython.Magic import Magic
import IPython.ipapi
ip = IPython.ipapi.get()
ip.magic('%colors linux')
import numpy as n
from cosmo_units import *
import time
################################################################
##     Parse inputs
####################################
o = optparse.OptionParser()
o.set_usage('bash_selfcal_test.py [options] *.ms')
o.set_description(__doc__)
a.scripting.add_standard_options(o,chan=True)
o.add_option('--uvsize',default=200,type='float',
    help='size of the uv grid in wavelengths [200]')
o.add_option('--uvres',default=4,type='float',
    help='size of the uv cells in wavelengths [4]')
o.add_option('--fconfig',default='120_180_6',
    help='Start_stop_step for output uveta in MHz.[120_180_6]')
o.add_option('--fspace',default=1,type='float',
    help='spacing between band centers = bw/fspace. [1]')
o.add_option('--dodiff',action='store_true',
    help='Also output power spectra of time diffed data.')
o.add_option('--clean', dest='clean', type='float', default=1e-3,
    help='Deconvolve delay-domain data by the response that results from flagged data.  Specify a tolerance for termination (usually 1e-2 or 1e-3).')
#o.add_option('--ubins_range',default='45_200',
#    help="""Power spectrum will be integrated in logarithmic radial bins. 
#    Specify radial limits in wavelengths. [45_200]""")
#o.add_option('--nu_bins',default=6,
#    help="Number of radial bins in which to bin [6]")
#clean off the casa args:
for i in range(len(sys.argv)):
    if sys.argv[i]==inspect.getfile( inspect.currentframe()):break
opts, args = o.parse_args(sys.argv[i+1:])
###################################################


uvsize=opts.uvsize
uvres =opts.uvres
#umin = 45
#umax = 200
#nu = 6
DIM = int(uvsize/uvres)
ulim = [1.,2000.]

pm = {'L': None, 'taps':None, 'fwidth':None, 'window':None, 'window_name':None, 'sinx_x':None, 'window_sinx_x':None}

def __set_pm__(L, window, taps, fwidth):
    global pm
    if pm['L'] == L and pm['taps'] == taps and pm['fwidth'] == fwidth:
        if type(window) == str and pm['window_name'] == window: return
        elif window is pm['window']: return
    else:
        pm['L'] = L
        pm['taps'] = taps
        pm['fwidth'] = fwidth
        def sinx_x(x):
            t = n.pi * taps * fwidth * (x/float(L) - .5)
            v = n.where(t != 0, t, 1)
            return n.where(t != 0, n.sin(v) / v, 1)
        pm['sinx_x'] = n.fromfunction(sinx_x, (L,))
    if type(window) == str:
        wf = {}
        wf['hamming'] = lambda x: .54 + .46 * cos(2*n.pi*x/L - n.pi)
        wf['hanning'] = lambda x: .5 + .5 * n.cos(2*n.pi*x/(L+1) - n.pi)
        wf['none'] = lambda x: 1
        pm['window'] = n.fromfunction(wf[window], (L,))
        pm['window_name'] = window
    else:
        pm['window'] = window
        pm['window_name'] = None
    pm['window_sinx_x'] = pm['window'] * pm['sinx_x']

def __pfb_fir__(data, window='hamming', taps=8, fwidth=1):
    L = data.shape[-1]
    __set_pm__(L, window, taps, fwidth)
    d = data * pm['window_sinx_x']
    print d.shape,taps,L/taps
    try: d.shape = d.shape[:-1] + (taps, L/taps)
    except: raise ValueError("More taps than samples")
    return n.sum(d, axis=len(d.shape) - 2)

def pfb(data, window='hamming', taps=8, fwidth=1, fft=n.fft.fft,axis=-1):
    """Perform PFB on last dimension of 'data' for multi-dimensional arrays.
    'window' may be a name (e.g. 'hamming') or an array with length of the
    last dimension of 'data'.  'taps' is the number of PFB taps to use.  The
    number of channels out of the PFB will be length out of the last 
    dimension divided by the number of taps. 'fwidth' scales the width of 
    each channel bandpass (to create overlapping filters, for example)."""
    return fft(__pfb_fir__(data, window, taps, fwidth))
NTAPS = 2
#highchan =int(opts.chan.split('_')[1])
#lowchan = int(opts.chan.split('_')[0])
fstart = int(opts.fconfig.split('_')[0])*1e6
fstop = int(opts.fconfig.split('_')[1])*1e6
bw = int(opts.fconfig.split('_')[2])*1e6
print opts.fconfig,fstart,fstop,bw
c = 3e8 #m/s
from numpy import array,int32
csystemplate = {'linear0': {'axes': array(['UU', 'VV'], 
      dtype='|S3'),
             'cdelt': array([-2.29183118,  2.29183118]),
             'crpix': array([ 250.,  250.]),
             'crval': array([ 0.,  0.]),
             'pc': array([[ 1.,  0.],
       [ 0.,  1.]]),
             'units': array(['lambda', 'lambda'], 
      dtype='|S7')},
 'linear2': {'axes': array(['Time'], 
      dtype='|S5'),
             'cdelt': array([  2.04834084e-07]),
             'crpix': array([ 50.]),
             'crval': array([ 0.]),
             'pc': array([[ 1.]]),
             'units': array(['s'], 
      dtype='|S2')},
 'obsdate': {'m0': {'unit': 'd', 'value': 55455.173560246825},
             'refer': 'UTC',
             'type': 'epoch'},
 'observer': '',
 'pixelmap0': array([0, 1], dtype=int32),
 'pixelmap1': array([2], dtype=int32),
 'pixelmap2': array([3], dtype=int32),
 'pixelreplace0': array([ 0.,  0.]),
 'pixelreplace1': array([ 0.]),
 'pixelreplace2': array([ 0.]),
 'pointingcenter': {'initial': False,
                    'value': array([ 1.3962    , -0.53619133])},
 'stokes1': {'axes': array(['Stokes'], 
      dtype='|S7'),
             'cdelt': array([ 1.]),
             'crpix': array([ 0.]),
             'crval': array([ 1.]),
             'pc': array([[ 1.]]),
             'stokes': array(['I'], 
      dtype='|S2')},
 'telescope': 'PAPER',
 'worldmap0': array([0, 1], dtype=int32),
 'worldmap1': array([2], dtype=int32),
 'worldmap2': array([3], dtype=int32),
 'worldreplace0': array([ 0.,  0.]),
 'worldreplace1': array([ 1.]),
 'worldreplace2': array([ 0.])
}



def grid_it(im,us,vs,ws,ds,wgts):
    #print 'Gridding %d integrations' % n_ints
    sys.stdout.write('|'); sys.stdout.flush()
    if len(ds) == 0: raise ValueError('No data to use.')
    ds,wgts = n.concatenate(ds), n.concatenate(wgts).flatten()
    us,vs,ws = n.concatenate(us), n.concatenate(vs), n.concatenate(ws)
    # Grid data into UV matrix
    (us,vs,ws),ds,wgts = im.append_hermitian((us,vs,ws),ds,wgts)
    im.put((us,vs,ws), ds, wgts)
    #im.put((us,vs,ws), ds, wgts, invker2=bm_im)
def plop(X,x):
    return n.argwhere(n.min(n.abs(X-x))==n.abs(X-x)).squeeze()
def length(X): return n.sqrt(n.dot(X,X))
def MJDs2JD(t):
    return t/86400 + 2400000.5
def JD2MJDs(t):
    return (t-2400000.5)*86400
def difftime(D,times):
    times = list(set(times))
    inshape = D.shape
    D.shape = D.shape[:-1] + (D.shape[-1]/len(times),len(times))
    dD = n.diff(D,axis=(D.ndim-1))
    dD = n.concatenate((dD,n.zeros(dD.shape[:-1]+(1,))),axis=(dD.ndim-1))
    return n.reshape(dD,inshape)
def diffchan(D):
    dD = n.diff(D,axis=1)
    return n.concatenate((dD,n.zeros((D.shape[0],1,D.shape[2]))),axis=1)
def difftc(D,times):
    dD = n.diff(D,axis=1)
    dD = difftime(dD,times)
    return n.concatenate((dD,n.zeros((D.shape[0],1,D.shape[2]))),axis=1)
def CLEAN_uvs(uvs):
    inshape = uvs.shape
    uvs.shape = (uvs.shape[0]*uvs.shape[1],uvs.shape[2])
    M = n.logical_not(uvs.mask).astype(n.float)
    uvs_c = n.zeros_like(uvs)
    for i in range(M.shape[0]): #iterate over uv
        val = M[i]
        if not val.sum()>0:continue
        d = uvs[i]
        ker = n.fft.ifft(val)
        _d = n.fft.ifft(d)
        _d, info = a.deconv.clean(_d, ker, tol=opts.clean)
        # Not sure if dividing by gain here is the right thing... once clean components are removed, want
        # residuals to be in original data units.
        #if True: _d = info['res'] / gain
        if True: _d = info['res']
        else: _d += info['res'] / gain
        d = n.fft.fft(_d) * val
        uvs_c[i] = d
    uvs_c = n.reshape(n.ma.array(uvs_c,mask=uvs.mask),inshape)
    return uvs_c

flush = sys.stdout.flush
for vis in args:
#    outfile = vis+'_%d_%d.uveta'%(lowchan,highchan)
    outfile = vis +'.uveta'
    if opts.dodiff: outfile +='.diff'
    print vis,'->',outfile
    print "Analyzing %s "%(vis)
    flush()
    ms.open(vis)
    rec = ms.getdata(['axis_info'])
    F = rec['axis_info']['freq_axis']['chan_freq'].squeeze()
    df = rec['axis_info']['freq_axis']['resolution'].squeeze()[0]
    ms.close()
    fs = n.arange(fstart,fstop,bw/opts.fspace)
    nchan = int(bw/df) + int(bw/df)%2
    print fstart,fstop,bw
    uvetastack = []
    #start channel loop
    for fmin in fs:
        tstart = time.time()
        lowchan = plop(F,fmin)
        highchan = lowchan + nchan
        if highchan>(len(F)-1):continue
        ms.open(vis)
        ms.selectinit()
        ms.select({'uvdist':ulim})
        #ms.iterinit(columns=['TIME'])
        print lowchan,highchan,nchan
        ms.selectchannel(highchan-lowchan,lowchan,1,1)
        rec = ms.getdata(['axis_info'])
        f = n.median(rec['axis_info']['freq_axis']['chan_freq'].squeeze())
        z = f212z(f)
                #moredata=True
        #while(moredata):
        rec = ms.getdata(['u','v','w','data','flag','antenna1','antenna2','time'])
        t = n.median(rec['time'])
        I,J = rec['antenna1'],rec['antenna2']
        D = n.ma.array(rec['data'],mask=rec['flag'])
        if opts.dodiff: D = difftc(D,rec['time'])
        D = D.squeeze()
#        print type(D),D.mask.sum(),D.size
        U,V,W = rec['u']*f/c,rec['v']*f/c,rec['w']*f/c
        del(rec)
        ms.close()
        tb.open(vis+'/FIELD')
        direction = tb.getcol('PHASE_DIR').squeeze()
        tb.close()
        
        
        
        print " done"
        flush()
        
        print "match uvw with ijs",;flush()
        Nant = n.max(n.vstack((I,J)))
        UVWs = n.zeros((Nant,Nant,3))
        bls = {}
        #for i,j in n.indices(UVWs.shape):
        for (i,j,u,v,w) in zip(I,J,U,V,W):
            bl = '%d&&%d'%(i,j)
            if not bls.has_key(bl): bls[bl] = []
            bls[bl].append((u,v,w))
        for bl in bls:
            BL = n.array(bls[bl])
            bls[bl] = (n.median(BL[:,0]),n.median(BL[:,1]),n.median(BL[:,2]))
        print ". done";flush()
        #TODO grid the single set of uvws, use get_indices()
        #TODO Then grid by the incoherent uvs radially.  
        ci = [] #coherent indices (into Im a uv grid)
        ici = [] #incoherent indices (into Ps, a list of radial us)
        
        Im = a.img.Img(uvsize, uvres, mf_order=0)
    #    Ps = n.logspace(n.log10(umin),n.log10(umax),num=(nu+2))  #we'll eventually throw out 0 & -1
        print "badly written, but exact, gridding",;flush()
        for i,j in zip(I,J):
            bl = '%d&&%d'%(i,j)
            u,v = bls[bl][0],bls[bl][1]
            uv = n.array([u,v])
            ci.append(Im.get_indices(u,v))
        uvs = n.ma.zeros(Im.uv.shape+(D.shape[0],))    #uvf cube
        uvin = n.zeros_like(uvs).astype(n.int)
        for l,(ui,vi) in enumerate(ci):
            if ui<0:continue
            uvs[ui,vi,:] += D[:,l]
#            print uvs.sum(),D[:,l].sum(),':',
            uvin[ui,vi,:] += n.logical_not(D[:,l].mask)
        print ".. done";flush()
        uvs[uvin>0] /= uvin[uvin>0]
#        print uvs.sum()
#        uvs = n.ma.array(uvs,mask=n.logical_not(uvin).astype(n.bool))
        #STOP! Save the uvgrid and exit.

        nchan = uvs.shape[2]
        print "CLEANing before max=%f"%(uvs.max());flush()
        if opts.clean<1: uvs = CLEAN_uvs(uvs)
        print "after max=%f"%(uvs.max());flush()
        print "FFT";flush()
        print uvs.shape
#        uveta = n.abs(n.fft.ifft(uvs,axis=2))[:,:,:nchan/2]
        uveta = n.abs(pfb(uvs,taps=NTAPS,window='hanning', fft=n.fft.ifft))
        uveta = a.img.recenter(uveta,(uvs.shape[0]/2,uvs.shape[1]/2,0))
        uvetastack.append(uveta)
        print uveta.shape
        print "t = %6.1fs"%(time.time()-tstart,)
    #end spectral loop
    uvetastack = n.array(uvetastack)
    print uvetastack.shape
    uvetastack = n.transpose(uvetastack,axes=[1,2,0,3])
    print uvetastack.shape
    dly = n.fft.fftfreq(uveta.shape[2],df)[:nchan/2]
    csupdate = {'linear0':{'cdelt':n.array([-uvres, uvres]).astype(n.float),
                            'crpix':n.array((uveta.shape[0]/2,uveta.shape[1]/2)).astype(n.float)},
                'obsdate':{'m0':{'unit':'d','value':t/86400},
                            'refer':'UTC','type':'epoch'},
                'pointingcenter':{'initial':False,
                                    'value':direction
                           },
                'linear2':{'cdelt':array([n.diff(dly)[0]]),
                            'crpix':array([0.]),
                            'crval':array([dly[0]])},    
                  }
    csysrecord = csystemplate
    for k in csupdate:
        csysrecord[k].update(csupdate[k])
    csysrecord['observer']=str(z)
#    csys = cs
#    csys.fromrecord(csysrecord)
#    print csys.coordinatetype()
    newspec = cs.newcoordsys(spectral=T)
    newspec.setreferencelocation(pixel=0,world=n.round(n.median(fs.min()),decimals=1))
    newspec.setincrement(value=bw/opts.fspace,type='spectral')
    csys = cs
    csys.fromrecord(csysrecord)
    csys.replace(newspec.torecord(),0,1)
    print "writing ",outfile;flush()
    uveta.shape = (uveta.shape[0],uveta.shape[1],1,uveta.shape[2])
    ia.fromarray(outfile=outfile,pixels=uvetastack,csys=csys.torecord(),overwrite=True)
#    ia.open(outfile)
#    csys = ia.coordsys()
#    csys.replace(newspec.torecord(),0,1)
#    ia.setcoordsys(csys.torecord())
#    ia.close()
