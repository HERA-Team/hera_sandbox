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
from capo.pfb import *
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
o.add_option('--data_type',default='corrected',type='str',
    help="Can be: raw,corrected,model, or residual")
o.add_option('--name',default=None,type='str',
    help="optional name to insert in filename")
o.add_option('--window',default='hamming',
    help="can be any of: "+','.join(WINDOWS))
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

#def __set_pm__(L, window, taps, fwidth):
#    global pm
#    if pm['L'] == L and pm['taps'] == taps and pm['fwidth'] == fwidth:
#        if type(window) == str and pm['window_name'] == window: return
#        elif window is pm['window']: return
#    else:
#        pm['L'] = L
#        pm['taps'] = taps
#        pm['fwidth'] = fwidth
#        def sinx_x(x):
#            t = n.pi * taps * fwidth * (x/float(L) - .5)
#            v = n.where(t != 0, t, 1)
#            return n.where(t != 0, n.sin(v) / v, 1)
#        pm['sinx_x'] = n.fromfunction(sinx_x, (L,))
#    if type(window) == str:
#        wf = {}
#        wf['hamming'] = lambda x: .54 + .46 * cos(2*n.pi*x/L - n.pi)
#        wf['hanning'] = lambda x: .5 + .5 * n.cos(2*n.pi*x/(L+1) - n.pi)
#        wf['none'] = lambda x: 1
#        pm['window'] = n.fromfunction(wf[window], (L,))
#        pm['window_name'] = window
#    else:
#        pm['window'] = window
#        pm['window_name'] = None
#    pm['window_sinx_x'] = pm['window'] * pm['sinx_x']
#
#def __pfb_fir__(data, window='hamming', taps=8, fwidth=1):
#    L = data.shape[-1]
#    __set_pm__(L, window, taps, fwidth)
#    d = data * pm['window_sinx_x']
##    print d.shape,taps,L/taps
#    try: d.shape = d.shape[:-1] + (taps, L/taps)
#    except: raise ValueError("More taps than samples")
#    return n.sum(d, axis=len(d.shape) - 2)
#
#def pfb(data, window='hamming', taps=8, fwidth=1, fft=n.fft.fft,axis=-1):
#    """Perform PFB on last dimension of 'data' for multi-dimensional arrays.
#    'window' may be a name (e.g. 'hamming') or an array with length of the
#    last dimension of 'data'.  'taps' is the number of PFB taps to use.  The
#    number of channels out of the PFB will be length out of the last 
#    dimension divided by the number of taps. 'fwidth' scales the width of 
#    each channel bandpass (to create overlapping filters, for example)."""
#    return fft(__pfb_fir__(data, window, taps, fwidth))
NTAPS = 2
#highchan =int(opts.chan.split('_')[1])
#lowchan = int(opts.chan.split('_')[0])
fstart = int(opts.fconfig.split('_')[0])*1e6
fstop = int(opts.fconfig.split('_')[1])*1e6
bw = int(opts.fconfig.split('_')[2])*1e6
print opts.fconfig,fstart,fstop,bw
#form up the data request word
if opts.data_type=='raw':
    datastr = 'data'
else:
    datastr = opts.data_type+'_'+'data'
print "Getting column: %s"%(datastr,)

nbins = 2
favg = 4
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
def squish(A,factor):
    out = n.zeros(A.size/factor)
    for l in range(len(out)-1):
        out[l] = n.median(A[l*factor:l*factor+factor])
    out[-1] = n.median(A[-factor:])
    return out
        
flush = sys.stdout.flush
for vis in args:
#    outfile = vis+'_%d_%d.uveta'%(lowchan,highchan)
    if not opts.name is None:
        outfile = vis +'.'+opts.name+'.uveta'
    else:
        outfile = vis+'.uveta'
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
    print fs
    nchan = int(bw/df)
    print fstart,fstop,bw
    uvetastack = []
    uveta_powerstack = []
    #start channel loop
    for fmed in fs:
        #calculate channels, start the timer clock
        tstart = time.time()
        medchan = plop(F,fmed)
#        lowchan = medchan - (nchan + nchan*NTAPS)/2
#        highchan = medchan + (nchan + nchan*NTAPS)/2
        lowchan = medchan - int(1./3*(NTAPS*nchan)) - int(1./3*NTAPS*nchan)%NTAPS
        highchan = medchan + int(2./3*(NTAPS*nchan)) + int(2./3*NTAPS*nchan)%NTAPS
        if highchan>(len(F)-1):continue

        #get data
        print "loading data",;flush()
        ms.open(vis)
        ms.selectinit()
        ms.select({'uvdist':ulim})
        #ms.iterinit(columns=['TIME'])
        print lowchan,highchan,highchan-lowchan,nchan,df
        ms.selectchannel(highchan-lowchan,lowchan,1,1)
        rec = ms.getdata(['axis_info'])
        f = n.median(rec['axis_info']['freq_axis']['chan_freq'].squeeze())
        z = f212z(f)
                #moredata=True
        #while(moredata):
        rec = ms.getdata(['u','v','w',datastr,'flag','antenna1','antenna2','time'])
        t = n.median(rec['time'])
        I,J = rec['antenna1'],rec['antenna2']
        D = n.ma.array(rec[datastr],mask=rec['flag'])
        if opts.dodiff: D = difftime(D,rec['time'])
        D = D.squeeze()
#        print type(D),D.mask.sum(),D.size
        U,V,W = rec['u']*f/c,rec['v']*f/c,rec['w']*f/c
        del(rec)
        ms.close()

        #get the field center
        tb.open(vis+'/FIELD')
        direction = tb.getcol('PHASE_DIR').squeeze()
        tb.close()       
        print " done"
        flush()
        
#        print "match uvw with ijs",;flush()
#        Nant = n.max(n.vstack((I,J)))
#        UVWs = n.zeros((Nant,Nant,3))
#        bls = {}
#        #for i,j in n.indices(UVWs.shape):
#        #find the median uvws over freq and time
#        for (i,j,u,v,w) in zip(I,J,U,V,W):
#            bl = '%d&&%d'%(i,j)
#            if not bls.has_key(bl): bls[bl] = []
#            bls[bl].append((u,v,w))
#        for bl in bls:
#            BL = n.array(bls[bl])
#            bls[bl] = (n.median(BL[:,0]),n.median(BL[:,1]),n.median(BL[:,2]))
#        print ". done";flush()
#
#
#        ci = [] #coherent indices (into Im a uv grid)
#        Im = a.img.Img(uvsize, uvres, mf_order=0)
#        print "badly written, but exact, grid association",;flush()
#        for i,j in zip(I,J):
#            bl = '%d&&%d'%(i,j)
#            u,v = bls[bl][0],bls[bl][1]
#            uv = n.array([u,v])
#            ci.append(Im.get_indices(u,v))
#        print "done"

        print "PFB [%s]. "%(opts.window,),;flush()
        FD = n.zeros((D.shape[0]/NTAPS,D.shape[1])).astype(n.complex64)
        for l in range(D.shape[1]):
            FD[:,l] = pfb(D[:,l],taps=NTAPS,window=opts.window, 
                        fft=n.fft.ifft)
        print "done"
        Im = a.img.Img(uvsize, uvres, mf_order=0)
        uveta = n.zeros(Im.uv.shape+(D.shape[0]/NTAPS,)).astype(n.complex64)    #uvf cube
#        uvetas = [uveta.copy() for i in range(nbins)]
        uvetan = n.zeros_like(uveta)
#        uvs = n.ma.zeros(Im.uv.shape+(D.shape[0],)).astype(n.complex64)
#        uveta_all = n.zeros_like(uveta).astype(n.complex64)
#        uveta_buff = n.zeros_like(uveta).astype(n.complex64)
#        uvin = n.zeros_like(uvs).astype(n.int)
        del(D)
        #grid by cross multiplying FTd samples
        print 'Gridding.',;flush()
#        print "gridding into %d bins"%(nbins,)
#        for l,(ui,vi) in enumerate(ci):
        for l,(u,v) in enumerate(zip(U,V)):
            ui,vi = Im.get_indices(u,v)
#            if ui<0:continue
            uveta[ui,vi,:]  +=  FD[:,l]
            uvetan[ui,vi,:] +=  1
        print ".. done";flush()
        uveta[uvetan>0] /= uvetan[uvetan>0]
        print uveta.shape
        uveta = a.img.recenter(uveta,
            (uveta.shape[0]/2,uveta.shape[1]/2,0))
        uvetastack.append(n.fft.fftshift(uveta,axes=(2,)))

        uveta = uveta*n.conj(uveta)
        #average +k and -k parts, take take the abs
        print uveta.shape[2]/2
        uveta_power = n.abs((uveta[:,:,:uveta.shape[2]/2] + \
            uveta[:,:,:-uveta.shape[2]/2-(uveta.shape[2]-1)%2:-1])/2)

        print "uveta power min,max",uveta_power.min(),uveta_power.max()
        if uveta_power.min()<0: print "WARNING: ",uveta_power.min()
        uveta_powerstack.append(uveta_power)

        dly = n.fft.fftfreq(nchan,df)[:uveta.shape[2]/2]
        print nchan,df,dly[:3]
        print "t = %6.1fs"%(time.time()-tstart,)


    #end spectral loop
    print "uveta0 min",uvetastack[0].min()
    uvetastack = n.array(uvetastack).astype(n.complex64)
    uveta_powerstack = n.array(uveta_powerstack)
    print uvetastack.shape
    uvetastack = n.transpose(uvetastack,axes=[1,2,0,3])
    uveta_powerstack = n.transpose(uveta_powerstack,axes=[1,2,0,3])
    print uvetastack.shape
    print "full uveta data type",uvetastack.dtype
    print "uveta min",uvetastack.min()

    #save the power
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

    newspec = cs.newcoordsys(spectral=T)
    newspec.setreferencelocation(pixel=0,world=n.round(n.median(fs.min()),decimals=1))
    newspec.setincrement(value=bw/opts.fspace,type='spectral')
    csys = cs
    csys.fromrecord(csysrecord)
    csys.replace(newspec.torecord(),0,1)
    print "writing ",outfile;flush()
#    uveta.shape = (uveta.shape[0],uveta.shape[1],1,uveta.shape[2])
    ia.fromarray(outfile=outfile,pixels=uveta_powerstack,csys=csys.torecord(),overwrite=True)

    ########################################################
    #save the complex uveta
    csupdate = {'linear0':{'cdelt':n.array([-uvres, uvres]).astype(n.float),
                            'crpix':n.array((uveta.shape[0]/2,uveta.shape[1]/2)).astype(n.float)},
                'obsdate':{'m0':{'unit':'d','value':t/86400},
                            'refer':'UTC','type':'epoch'},
                'pointingcenter':{'initial':False,
                                    'value':direction
                           },
                'linear2':{'cdelt':array([n.diff(dly)[0]]),
                            'crpix':array([0.]),
                            'crval':array([dly.min()])},    
                  }
    csysrecord = csystemplate
    for k in csupdate:
        csysrecord[k].update(csupdate[k])
    csysrecord['observer']=str(z)

    newspec = cs.newcoordsys(spectral=T)
    newspec.setreferencelocation(pixel=0,world=n.round(n.median(fs.min()),decimals=1))
    newspec.setincrement(value=bw/opts.fspace,type='spectral')
    csys = cs
    csys.fromrecord(csysrecord)
    csys.replace(newspec.torecord(),0,1)
    print "writing ",outfile+'.im,.re';flush()
#    uveta.shape = (uveta.shape[0],uveta.shape[1],1,uveta.shape[2])
    ia.fromarray(outfile=outfile+'.re',pixels=n.real(uvetastack),csys=csys.torecord(),overwrite=True)
    ia.fromarray(outfile=outfile+'.im',pixels=n.imag(uvetastack),csys=csys.torecord(),overwrite=True)



print "\a"
#    ia.open(outfile)
#    csys = ia.coordsys()
#    csys.replace(newspec.torecord(),0,1)
#    ia.setcoordsys(csys.torecord())
#    ia.close()
