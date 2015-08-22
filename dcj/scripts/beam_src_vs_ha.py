#! /usr/bin/env python
import matplotlib
matplotlib.use('Agg')
import aipy as a, numpy as n
import optparse, sys, scipy.optimize
import capo as C
import pylab as p
import ipdb
from pyslalib import slalib
def sky_sep(A,B):
    #compute distance on sphere
    #using input vectors in xyz
    (theta1,phi1) = slalib.sla_cc2s(A)
    (theta2,phi2) = slalib.sla_cc2s(B)
    return slalib.sla_sep(theta1,phi1,theta2,phi2)
o = optparse.OptionParser()
a.scripting.add_standard_options(o, cal=True, chan=True, pol=True)
o.add_option('-s', '--src', dest='src', type='str', 
    help='Source to use for calibration.')
o.add_option('--cat', dest='cat', type='str', default='helm,misc',
    help='A comma-delimited list of catalogs from which sources are to be drawn.  Default is "helm,misc".  Other available catalogs are listed under aipy._src.  Some catalogs may require a separate data file to be downloaded and installed.')
o.add_option('--sep_min',default=2,type='float',
    help="maximum allowed calibrator distance in degrees default=2")
o.add_option('--plot_src_track',action='store_true',
    help='Does what it says.')
opts,args = o.parse_args(sys.argv[1:])

def filename2src(f):
    #return f.split('_')[-1]
    return f.split('_')[-1].split('.')[0]

uv = a.miriad.UV(args[0])
chans = a.scripting.parse_chans(opts.chan, uv['nchan'])
aa = a.cal.get_aa(opts.cal, uv['sdf'], uv['sfreq'], uv['nchan'])
aa.select_chans(chans)

if opts.src != None:
    srclist,cutoff,catalogs, = a.scripting.parse_srcs(opts.src, opts.cat)
    calsrc = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)
    calsrc = calsrc.values()[0]
else:
    calsrc = None
srclist = [filename2src(f) for f in args]

#load a catalog of the sources found in the file
if not calsrc is None: assert(calsrc.src_name in srclist)
srclist,cutoff,catalogs, = a.scripting.parse_srcs(','.join(srclist), opts.cat)
cat = a.cal.get_catalog(opts.cal, srclist, cutoff, catalogs=catalogs)

tracks = {}

for filename in args:
    print filename
    srcname = filename2src(filename)
    if not tracks.has_key(srcname):
        tracks[srcname] = {}
        tracks[srcname]['top'] = []
        tracks[srcname]['bm'] = []
        tracks[srcname]['dat'] = []
        tracks[srcname]['wgt'] = []
        tracks[srcname]['azalt'] = []
    src = cat[srcname]
    uv = a.miriad.UV(filename)
    a.scripting.uv_selector(uv, 'cross', opts.pol)
    for (uvw,t,(i,j)),d,f in uv.all(raw=True):
        d,f = d[chans], f[chans]
        f = n.where(d == 0, 1, f)
        f = n.where(n.isnan(d), 1, f)
        d = n.where(f, 0, d)
        w = n.logical_not(f)
        #if False and w.sum() > 0: # Do delay filtering
        #    DLYWID = 5
        #    dmdl,dres = C.dspec.wideband_dspec(d, w, DLYWID+1, -DLYWID, tol=1e-4, window='none')
        #    dres[DLYWID+1:-DLYWID] = 0
        #    d = (dmdl+dres) * w
        aa.set_jultime(t)
        src.compute(aa)
        src_xyz = src.get_crds('top')
        tracks[srcname]['top'].append(src_xyz)
        src_azalt = src.get_crds('top',ncrd=2)
        tracks[srcname]['azalt'].append(src_azalt)
        bmi = aa[i].bm_response(src_xyz, pol=opts.pol[0])
        bmj = aa[j].bm_response(src_xyz, pol=opts.pol[-1])
        tracks[srcname]['bm'].append(bmi * n.conj(bmj))
        #print src.src_name, src.ra, src.dec, src_xyz#, bm[-1]
        tracks[srcname]['dat'].append(d)
        tracks[srcname]['wgt'].append(w)



calsrc = calsrc.src_name

calsrc_x   = n.array(tracks[calsrc]['top'])[:,0]
calsrc_bm  = n.array(tracks[calsrc]['bm'])[...,0]
calsrc_dat = n.array(tracks[calsrc]['dat'])
calsrc_wgt = n.array(tracks[calsrc]['wgt'])

print calsrc_x.shape, calsrc_x.max()
print calsrc_bm.shape, calsrc_bm.max()
print calsrc_dat.shape, calsrc_dat.max()
print calsrc_wgt.shape, calsrc_wgt.max()

def pwrlaw(fqs, flx, ind, mfreq=.150):
    return flx * (fqs/mfreq)**ind
    
def fit_pwrlaw(fqs, dat, err, flx, ind, mfreq=.150):
    spec = pwrlaw(fqs, flx, ind, mfreq=mfreq)
    scr = n.abs(dat-spec)**2 / err**4
    wgt = 1. / err**4
    if False: # favor middle of band
        scr = n.where(n.abs(fqs - .150) > .25, 0, scr)
        wgt = n.where(n.abs(fqs - .150) > .25, 0, wgt)
    rv = n.sqrt(n.sum(scr)/n.sum(wgt))
    #print flx, ind, rv
    return rv

calxyz = n.array(tracks[calsrc]['top'])
badsrc = []
import pylab as p
for cnt,srcname in enumerate(tracks.keys()):
    color = 'krbgcm'[cnt % 6]
    x = n.array(tracks[srcname]['top'])[:,0]
    alt = n.array(tracks[srcname]['azalt'])[:,1]
    xyz = n.array(tracks[srcname]['top'])
    ha = n.array(tracks[srcname]['azalt'])[:,0]
    _cbm, _cd, _cw = [], [], []
    _cdist = []
    _too_far = []
    print "working on ",srcname
    print "finding overlap with calibrator",calsrc
    for i,_x in enumerate(x): # Find closest point bin the beam between two tracks
        skyseps = n.array([sky_sep(xyz[i],cxyz) for cxyz in calxyz])
        m = n.argmin(skyseps)
        _m = n.argmin(n.abs(calsrc_x - _x))
        _cbm.append(calsrc_bm[m])
        _cd.append(calsrc_dat[m])
        _cw.append(calsrc_wgt[m])
        _cdist.append(skyseps[m])
        #if n.abs(calsrc_x[m] - _x) >= .005 or \
        if skyseps[m]> (opts.sep_min*n.pi/180) or \
            alt[i]<0: # XXX this depends on time step
            _cw[-1] *= 0
        #if skyseps[m]>(opts.sep_min*n.pi/180):_too_far.append(skyseps[m])
        if n.abs(calsrc_x[_m] - _x) >= 0.005: _too_far.append(skyseps[m])
    if n.sum(_cw)==0:
        print "zero(!) calibrator overlap for %s"%(srcname)
        badsrc.append(srcname)
        continue
    mask = n.array(_cw) * n.array(tracks[srcname]['wgt']) # flag data that are flagged for either source
    _cdist = n.ma.masked_where(n.mean(_cw,axis=1)>0,_cdist)*180/n.pi
    print "hour angles:",set(n.ceil(ha[n.mean(_cw,axis=1)>0]*12/n.pi)-12)
    print "range of distances",_cdist.min(),_cdist.max()
    print "threw out %d data points for not having cal overlap"%len(_too_far)
    print _too_far
    #print "_cw.shape",n.array(_cw).shape
    tracks[srcname]['include'] = n.mean(_cw,axis=1)
    cbm = n.array(_cbm) * mask
    cd = n.array(_cd) * mask
    bm = n.array(tracks[srcname]['bm'])[...,0] 
    d = n.array(tracks[srcname]['dat']) * mask
    sum_w2 = n.sum((mask*bm)**2, axis=0); sum_w2 = n.where(sum_w2 == 0, 1, sum_w2)
    sum_cw2 = n.sum(cbm**2, axis=0); sum_cw2 = n.where(sum_cw2 == 0, 1, sum_cw2)
    srcest = n.sum(d * bm, axis=0) / sum_w2
    calest = n.sum(cd * cbm, axis=0) / sum_cw2
    #perform a beam-weighted rms on the residuals
    #errest = n.sqrt(n.sum(bm * ( bm * srcest - d ) * n.conj(bm * srcest -d),axis=0)/n.sum(bm,axis=0))
    _cc = n.argwhere(n.abs(ha-n.pi)<(10*n.pi/180))
    p.figure()
    p.plot(ha,(bm[:,100]*srcest[100] - d[:,100])/n.abs(d[:,100]))
    p.ylim([-1,1])
    p.savefig('beamslice_residual.png')
    errest = n.sqrt(n.sum(bm*(bm[_cc,:] * srcest - d[_cc,:])*n.conj(bm[_cc,:] * srcest - d[_cc,:]),axis=0)/sum_w2)
    mask_vs_fq = n.sum(mask, axis=0)
    valid = n.where(mask_vs_fq > mask_vs_fq.max() / n.sqrt(2), 1, 0)
    print "cleaning calibrator"
    print n.sum(n.isnan(calest)),n.sum(n.isnan(valid))
    if True: # Do delay filtering on calsrc
        DLYWID = 10
        dmdl,dres = C.dspec.wideband_dspec(calest*valid, valid, DLYWID+1, -DLYWID, tol=1e-10, window='none')
        dres[DLYWID+1:-DLYWID] = 0
        sm_calest = dmdl + dres
    
    srcest = srcest.compress(valid)
    calest = calest.compress(valid)
    sm_calest = sm_calest.compress(valid)
    errest = errest.compress(valid)
    fq = aa.get_afreqs().compress(valid)

    cal = {}
    cal['1932-464'] = 93.7 * (fq/.150)**-0.82 
    #cal['2331-416'] = 31.9 * (fq/.150)**-0.69
    cal['2331-416'] = 33.9 * (fq/.150)**-0.76
    calsrc_spec = cal[calsrc]
    #gain = calsrc_spec / calest
    sm_gain = calsrc_spec / sm_calest
    #sm_gain = n.polyfit(fq, gain, deg=7)
    #sm_gain = n.polyfit(fq, gain, deg=10)
    #sm_gain = n.polyval(sm_gain, fq)

    jy_src = srcest * sm_gain
    jy_cal = calest * sm_gain
    #err = .1 * n.abs(jy_src)
    #err = n.ones_like(jy_src)

    print 'Writing', srcname+'_spec.npz'
    n.savez(srcname+'_spec.npz', freq=fq, spec=jy_src,res=errest,gain=sm_gain)

        
    p.subplot(211)
    #p.plot(x, dat/gain, label=srcname)
    p.plot(fq, jy_src, color+'.', label=srcname)
    p.plot(fq, jy_cal, color+',', label=srcname)

    def fitfunc(prms): return fit_pwrlaw(fq, jy_src, 1, prms[0], prms[1])
    rv = scipy.optimize.fmin(fitfunc, n.array([0.,0]), full_output=1, disp=0, ftol=1e-4, xtol=1e-4)
    prms,score = rv[:2]
    flx,ind = prms
    print srcname, flx, ind, score
    
    spec = pwrlaw(fq, flx, ind)
    p.plot(fq, spec, color+'-', label=srcname)

    p.subplot(212)
    #p.plot(x, bm, label=srcname)
    if True: # plot gain gainection
        #p.plot(fq, gain, color+'.', label=srcname)
        p.plot(fq, sm_gain, color+'-', label=srcname)
        p.plot(fq, calsrc_spec / calest, color+'.', label=srcname)
    else: # plot beam pattern
        p.plot(100*mask_vs_fq, color+'-')
        #p.plot(fq, srcest, color+'.')
        p.plot(n.sum(d * bm, axis=0), color+'.')
        p.plot(100*sum_w2, color+'x')

        #src_amp_vs_fq = spec * sm_gain
        
if opts.plot_src_track:
    try:
        from mpl_toolkits.basemap import Basemap
    except(ImportError):
        from matplotlib.toolkits.basemap import Basemap
        #if its not here, I die
    try:
        map = Basemap(projection='ortho',lat_0=90,lon_0=180)
        radec_map = Basemap(projection='ortho',lat_0=-30,lon_0=0)
        p.figure()
        map.drawmapboundary()
        radec_map.drawparallels(n.arange(-90,90,30)[1:])
        radec_map.drawmeridians(n.arange(-180, 180, 30))
        print "plotting source tracks"
        for srcname in tracks:
            if srcname in badsrc:continue
            print srcname,
            azalt = n.array(tracks[srcname]['azalt'])*a.img.rad2deg
            include = n.array(tracks[srcname]['include'])
            print "masking %d/%d non-overlapping points"%(n.sum(include==0),include.size)
            lon = azalt[:,0]
            lat = azalt[:,1]
            sx,sy = map(lon,lat)
            if srcname == calsrc:
                map.plot(sx[include>0],sy[include>0],',k',label=srcname+'*')
            else:
                map.plot(sx[include>0],sy[include>0],',',label=srcname)
        if len(tracks)<10: p.legend(loc='upper right',numpoints=1,ncol=2)
    except:
        ipdb.set_trace()
    #ipdb.set_trace()
    try:
        p.savefig("%s_srctrack.png"%'_'.join(tracks.keys()))
    except(IOError):
        p.savefig("many_sources_srctrack.png")
        
    
p.subplot(212)
p.legend()
p.savefig(srcname+'_spec.png')
#p.show()
sys.exit(0)

bp_cal = {}
for src in [calsrc] + cat.values():
    if src is None: continue
    srcfiles = [f for f in args if f.startswith(src.src_name)]
    for cnt, filename in enumerate(srcfiles):
        filetype = filename[len(src.src_name):]
        color = colors[filetype]
        print 'Reading', filename
        _f = open(filename)
        f = n.load(_f)
        spec,afreqs = f['spec'].flatten(), f['freq'].flatten()
        _f.close()
        #dspec = spec - a.rfi.remove_spikes(spec, order=8)
        #sig = n.std(dspec)
        if True:
            valid = n.where(n.logical_and(afreqs > .120, afreqs < 0.170), 1, 0)
            spec = spec.compress(valid)
            afreqs = afreqs.compress(valid)
        src.update_jys(afreqs)
        bp = n.sqrt(spec / src.jys)
        bp_poly = n.polyfit(afreqs, bp, deg=opts.deg)
        if not calsrc is None and src.src_name == calsrc.src_name:
            print 'Calibrating to', src.src_name
            bp_cal[filetype] = bp_poly
        if opts.src is None:
            bp_fit = 1.
        else: bp_fit = n.polyval(bp_cal[filetype], afreqs).clip(.1,10)**2
        spec /= bp_fit
        
        src_poly = n.polyfit(n.log10(afreqs/src.mfreq), n.log10(spec), deg=1)
        n.set_printoptions(threshold=n.nan)
        print 'bp =', list(bp_poly)
        print "'%s':" % src.src_name + "{ 'jys':10**%f, 'index':  %f , }," % (src_poly[-1], src_poly[-2])
        print 'RMS residual:', n.sqrt(n.average((spec - 10**n.polyval(src_poly, n.log10(afreqs/src.mfreq)))**2))
        if n.all(spec <= 0): continue

        if not opts.quiet:
            p.loglog(afreqs, spec, color+'.', label='Measured')
            p.loglog(afreqs, 10**n.polyval(src_poly, n.log10(afreqs/src.mfreq)), color+'-', 
                label='Fit Power Law')
    if not opts.quiet:
        #p.loglog(afreqs, src.jys, color+':', label='%f, %s' % (src._jys, str(src.index)))
        p.xticks(n.arange(.1,.2,.02), ['100','120','140','160','180'])
        p.xlim(afreqs[0], afreqs[-1])
        p.ylim(3,3e3)
        p.grid()
        p.title(src.src_name)
        #p.legend()
p.show()
