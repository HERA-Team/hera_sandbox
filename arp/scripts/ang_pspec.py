#! /usr/bin/env python
import aipy as a, numpy as n, pylab as p
import capo as C
import sys, optparse
#from matplotlib import rc
#rc('font',**{'family':'serif','serif':['Times']})
#rc('text', usetex=True)

o = optparse.OptionParser()
o.add_option('-d', '--dat', dest='datfile',
    help='Output file from previous run of this script, to regenerate plots without redoing analysis')
opts,args = o.parse_args(sys.argv[1:])

# Generating C_l of EoR from Santos 2005
x = n.array([1e2,2e2,4e2,1e3,3e3,1e4])
y = n.array([1.,3,10,30,60,50])
C_l_sim = n.polyfit(n.log10(x), n.log10(y), deg=3)

# Generating C_l of synchrotron
def C_l_sync(L, freq):
    return L*(L+1)/(2*n.pi) * 700. * (1e3/L)**2.4 * (130./freq)**(2*2.8)

#f_re = re.compile(r'(.*_c(\d+)_(\d+)_).*')
#def ch2freq(chan, sfreq=0.1212890625, sdf=0.00029296875):
#    return (sfreq + chan*sdf) * 1e9

"""
# Sort files according to channel
if not opts.datfile is None:
    exec(open(opts.datfile).read())
else:
    pspec = {}
    chs = {}
    for filename in args:
        print '#', filename
        prefix, ch1, ch2 = f_re.match(filename).groups()
        ch1,ch2 = float(ch1),float(ch2)
        ch = (ch1 + ch2) / 2.
        chs[ch] = chs.get(ch, []) + [prefix]

    chs_list = chs.keys()
    chs_list.sort()
    
    # Loop over all channels:
    print 'pspec = {'
    for cnt, ch in enumerate(chs_list):
        print '# Working on CH=%f, FREQ=%f MHz' % (ch, ch2freq(ch)/1e6)

        def jy2t(chan=ch, sfreq=0.1212890625, sdf=0.00029296875):
            '''Convert Jy to K for this channel.'''
            freq = ch2freq(chan, sfreq=sfreq, sdf=sdf)
            lam = a.const.c / freq
            pb = n.polyval(pb_poly, freq*1e-9)
            return 1e-23 * lam**2 / (2 * a.const.k * pb)

        print '# Jy to T conversion factor:', jy2t()

        # Get the data
        uv1sq,uv2sq,bmsq = [],[],[]
        fov = None
        for pfx in chs[ch]:
            print '#', ch, pfx
            dima, kwds = a.img.from_fits(pfx+'00a.dim.fits')
            dimb, kwds = a.img.from_fits(pfx+'00b.dim.fits')
            dbma, kwds = a.img.from_fits(pfx+'00a.dbm.fits')
            dbmb, kwds = a.img.from_fits(pfx+'00b.dbm.fits')
            # Some parameters determining size of pixels/image
            dpx_dg = kwds['d_ra']
            dpx_rd = dpx_dg * a.img.deg2rad
            sh = dima.shape[:2]
            DIM = sh[0]
            # Generate coherent (1) and incoherent (2) maps 
            dim1,dim2 = (dima + dimb) / 2., (dima - dimb) / 2
            dbm = (dbma + dbmb) / 2.
            dim1,dim2,dbm = dim1.squeeze(), dim2.squeeze(), dbm.squeeze()
            # Limit the FoV used in PS measurement
            if fov is None:
                rpx = 10. / dpx_dg
                print '# FoV: rpx=%f' % (rpx)
                #fov = n.fromfunction(lambda x,y: gen_circ(x,y,r=rpx), sh)
                fov = circ(DIM,r=rpx)
                #fov = a.img.gaussian_beam(rpx, sh, center=(DIM/2,DIM/2))
            dim1 *= fov; dim2 *= fov; dbm *= fov
            # Transform into UV domain and square visibilities (now Jy^2)
            uv1 = a.img.recenter(n.fft.fft2(dim1), (DIM/2,DIM/2))
            uv2 = a.img.recenter(n.fft.fft2(dim2), (DIM/2,DIM/2))
            bm = a.img.recenter(n.fft.fft2(dbm), (DIM/2,DIM/2))
            uv1sq.append(n.abs(uv1)**2)
            uv2sq.append(n.abs(uv2)**2)
            bmsq.append(n.abs(bm)**2)

        uv1sq,uv2sq,bmsq = n.array(uv1sq),n.array(uv2sq),n.array(bmsq)

        ps1,ps2 = [],[]
        wgts,sigs,ells = [],[],[]
        
        # Draw rings of geometrically increasing radius
        #logstep = .04
        logstep = .08
        rpxs = 10**n.arange(n.log10(5.), n.log10(385.), logstep)
        for rpx in rpxs:
            # Choose a ring width (for averaging) ~1/2 step size between radii
            w = rpx * (10**logstep - 10**(-logstep)) / 4
            ell = r2ell(rpx, dpx_rd, DIM)
            rng = ring(DIM, rpx, w)
            while len(rng.shape) < len(uv1sq.shape): rng.shape = (1,) + rng.shape
            print '# R_px = %f, width = %f, ell = %f' % (rpx, w, ell)
            # Select a ring out of square visibilities
            uv1sq_r = uv1sq * rng
            uv2sq_r = uv2sq * rng
            bmsq_r = bmsq * rng
            # Sum square visibilities around ring and divide by their weighting
            # Giving a weighted average for a sample on the ring
            wgt1 = bmsq_r.sum()
            if wgt1 == 0: continue
            wgt2 = (bmsq_r**2).sum()
            ps1.append(uv1sq_r.sum())
            ps2.append(uv2sq_r.sum())
            # Estimate "sample variance" around ring by removing average
            # vis^2 and finding how variance is beating down
            sigsq_r = (uv1sq_r - ps1[-1] * bmsq_r)**2
            # Compute stddev of measurements of C_l
            sig = n.sqrt(sigsq_r.sum() / wgt2)
            # Degrade resolution to get # of actual independent samples
            degres = 10
            bmsq_rd = bmsq_r.copy()
            bmsq_rd.shape = (bmsq_rd.shape[0], DIM/degres,degres,DIM/degres,degres)
            bmsq_rd = bmsq_rd.sum(axis=4)
            bmsq_rd = bmsq_rd.sum(axis=2)
            wgt2 = (bmsq_rd**2).sum()
            # Beat down stddev by number of independent samples of C_l
            # Divide by 2 b/c only half of ring are independent measurements
            sig = sig * n.sqrt(wgt2 / wgt1**2 / 2)
            sigs.append(sig); wgts.append(wgt1); ells.append(ell)
        ps1,ps2 = n.array(ps1), n.array(ps2)
        ells,wgts,sigs = n.array(ells), n.array(wgts), n.array(sigs)

        # Convert power spectra in Jy^2 to K^2
        c_l1 = ps1 / wgts * jy2t()**2
        c_l2 = ps2 / wgts * jy2t()**2  
        sigs = sigs * jy2t()**2

        # Compute a "total error" that is combo of sample variance and
        # the computed power spectrum of thermal noise.  Do 2sig for 95% confidence
        errs = 2*n.sqrt(sigs**2 + c_l2**2)

        pspec[ch] = {
            'c_l':list(c_l1),
            'dc_l':list(c_l2),
            'errs':list(errs),
            'ells':list(ells),
        }
        print '\t%.1f : {' % (ch)
        print '\t\t"c_l": %s,' % (str(list(c_l1)))
        print '\t\t"dc_l": %s,' % (str(list(c_l2)))
        print '\t\t"errs": %s,' % (str(list(errs)))
        print '\t\t"ells": %s,' % (str(list(ells)))
        print '\t},'
    
    # Fit a line (in log space) to c_l to determine approx. power law
    #def pwrlaw(x):
    #    A, ind = x
    #    wgt = 1 / errs**2
    #    wgt /= wgt.sum()
    #    dif = (c_l1 - A*ells**ind) * wgt
    #    return n.sqrt((dif**2).sum())
    #A,ind = a.optimize(pwrlaw, n.array([0., 0.]))
"""
        
print 'Reading', args[0]
im,kwds = a.img.from_fits(args[0])
print kwds
print 'Reading', args[1]
bm,kwds = a.img.from_fits(args[1])

umag_px = 1. / (kwds['d_ra'] * a.img.deg2rad * im.shape[0])
print 'Resolution:', umag_px, 'per px'

print 'Currently assuming center frequency of 150MHz...'
uvs = n.fft.fft2(im.squeeze()) * C.pspec.jy2T(.150)
uvs = a.img.recenter(uvs, n.array(uvs.shape)/2)
bms = n.fft.fft2(bm.squeeze())
bms = a.img.recenter(bms, n.array(bms.shape)/2)

umag, Trms2, errs, wgts = C.pspec.Trms2_vs_umag(uvs, bms, umag_px, umin=1, logstep=.05)
        
#color = .8 * float(cnt) / len(chs_list)
#color = (color, color, color)
color = (0,0,0)
if True:
    errs_dn = n.where(errs >= Trms2-1e-10, Trms2-1e-10, errs)
    p.errorbar(umag, Trms2, [errs_dn, errs], fmt='-', color=color)
else:
    p.plot(umag, Trms2, '-', color=color)
p.plot(umag, wgts*1e10, 'k:')

ax = p.gca(); ax.set_xscale('log'); ax.set_yscale('log')
#p.xlim(60,1e3)
#p.xlim(1e2,1e3)
#p.ylim(1e-1,1e11)
p.xlabel(r'$u$', fontsize=20)
p.ylabel(r'$T_{\rm rms}^2\ \ [mK^2]$', fontsize=20)
p.show()
