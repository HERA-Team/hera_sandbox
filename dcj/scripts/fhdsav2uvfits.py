import time
from astropy.io import fits
import os,sys
import numpy as n
import progressbar
from scipy.io import readsav
import ipdb
"""
Converts the residual uv file output by FHD into an uvfits file 
"""

"""
fields we care about
data.file('DATE') = times  FHD equiv = obs.baseline_info.jdate
data.field('BASELINE') = ant2*256 + ant1 
data.field('DATA') = dimensions[blt,:,:,Nfreqs,Npol,[Re,Im,W]] # not sure what those two middle dims are
data.field('UU') etc bl vector in ns cal['cal']['uu'][0]*1e9

"""
template_file = 'template.uvfits'
templatehdulist = fits.open(template_file)
templatehdu = templatehdulist[0]
for savfile in sys.argv[1:]:
    #get the other poles
    pol = savfile.split('_')[-1][:2]
    if pol=='xx':
        uvfile_altpol = savfile.replace('xx','yy')
        pols = [0,1]
    elif pol=='yy':
        uvfile_altpol = savfile.replace('yy','xx')
        pols = [1,0]
    else:
        print "polarization not found in filename. skipping"
        continue
    if not os.path.exists(uvfile_altpol):
        print "pol file",uvfile_altpol,"not found. please find"
        continue
    paramfile = savfile.split('_')[0]+'_params.sav'
    if not os.path.exists(paramfile):
        print "error: paramfile=",paramfile,"not found. please find"
        continue

    weightfile = savfile.split('_')[0]+'_flags.sav'
    if not os.path.exists(weightfile):
        print "error: weightfile",weightfile,"not found, please find"
        continue
    print "loading:",savfile
    uvfile = readsav(savfile)
    ant1 =  uvfile['obs']['baseline_info'][0]['tile_a'][0]
    ant2 =  uvfile['obs']['baseline_info'][0]['tile_b'][0]
    data =  uvfile['vis_ptr']
    #times = uvfile['obs']['baseline_info'][0]['jdate'][0]
    baseline = ant2*256+ant1
    freqs = uvfile['obs']['baseline_info'][0]['freq'][0]

    print "loading alternate polarization",uvfile_altpol
    uv_altpol = readsav(uvfile_altpol)
    data_altpol = uv_altpol['vis_ptr']

    print "loading baselines from params file:",paramfile
    params = readsav(paramfile)
    U = params['params']['uu'][0]*1e9
    V = params['params']['vv'][0]*1e9
    W = params['params']['ww'][0]*1e9
    times = params['params']['time'][0]
    print "loading weights from :",weightfile
    flags = readsav(weightfile)
    weights = n.dstack([flags['flag_arr'][0],flags['flag_arr'][1]]) 
    #create the new fits file
    outdata = n.zeros((data.shape[0],1,1,data.shape[1],2,3))
    outdata[:,0,0,:,pols[0],0] = n.real(data)
    outdata[:,0,0,:,pols[0],1] = n.imag(data)
    outdata[:,0,0,:,pols[1],0] = n.real(data_altpol)
    outdata[:,0,0,:,pols[1],1] = n.imag(data_altpol)
    outdata[:,0,0,:,:,2] = weights
    
    dtype =[('UU', '>f4'), ('VV', '>f4'), ('WW', '>f4'), ('BASELINE', '>f4'),('DATE', '>f4'), ('DATA', '>f4',(1,1, 1, data.shape[1], 2, 3))]
    names = ['UU','VV','WW','BASELINE','DATE','DATA']
    ### Possibly do with fits.ColDefs
#    templatehdu.data.field('UU') = U
#    templatehdu.data.field('VV') = V
#    templatehdu.data.field('WW') = W
#    templatehdu.data.field('BASELINE') = baseline
#    templatehdu.data.field('DATE') = times
#    templatehdu.data.field('DATA') = outdata
    cols = []
    cols.append(fits.hdu.groups.Column(name='UU',format='E',array=U))
    cols.append(fits.hdu.groups.Column(name='VV',format='E',array=V))
    cols.append(fits.hdu.groups.Column(name='WW',format='E',array=W))
    cols.append(fits.hdu.groups.Column(name='BASELINE',format='E',array=baseline))
    cols.append(fits.hdu.groups.Column(name='DATE',format='E',array=times))
    cols.append(fits.hdu.groups.Column(name='DATA',format='EP',array=outdata))
    try:
        tbhdu=fits.new_table(fits.hdu.groups.ColDefs(cols))
    except StandardError as e:
        print e
        ipdb.set_trace()

    
    outfile = savfile.split('_')[0]+'_FHDres.uvfits'
    print "saving ",outfile

    tbhdu.writeto(outfile,clobber=True)
    sys.exit()
#    try:
#        #O = n.rec.fromarrays([U,V,W,baseline,times,outdata],names=names)
#        O = n.array([U,V,W,baseline,times,outdata],dtype=dtype) 
#        hdu = fits.PrimaryHDU(O)
#    except StandardError as e:
#        print e
#        ipdb.set_trace()
    #print O.dtype
    #add the spectrum info
    templatehdu.header['NAXIS4'] = len(freqs)
    templatehdu.header['CRVAL4'] = freqs[0]
    templatehdu.header['CDELT4'] = n.diff(freqs)[0]
    templatehdu.header['CRPIX4'] = 0
    outfile = savfile.split('_')[0]+'_FHDres.uvfits'
    print "saving ",outfile
    templatehdulist[0] = templatehdu
    templatehdulist.writeto(outfile)


