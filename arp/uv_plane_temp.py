#! /usr/bin/env python
"""
This is a general-purpose script for plotting simple FITS images.
"""

import aipy as a, sys, optparse, os, re
import numpy as n, pylab as p, ephem, math

o = optparse.OptionParser()
o.set_usage('uv_plane_temp.py [options] *.fits')
o.set_description(__doc__)
a.scripting.add_standard_options(o, chan=True, cmap=True, max=True, drng=True)
o.add_option('-m', '--mode', dest='mode', default='log',
    help='Plot mode can be log (logrithmic), lin (linear), phs (phase), real, or imag.')
o.add_option('-o', '--outfile', dest='outfile', default='',
    help='If provided, will save the figure to the specified file instead of popping up a window.')
o.add_option('-p', '--pol', dest='pol', type='int', default=0, 
    help='Polarization index if FITS file has multiple polarizations.  Default 0.')
o.add_option('--batch', dest='batch', action='store_true',
    help='Process files in batch mode (one plot each) and output to a <input file>.png file')
o.add_option('--fov', dest='fov', type='float', 
    help='Radius of field of View (in degrees) to limit images to before doing power-spectrum analysis.  Used to enforce the flat-sky approximation inherent to mapping the UV-plane to the angular power spectrum.  Default none')
opts, args = o.parse_args(sys.argv[1:])

sqrt2 = n.sqrt(2)
def circ(dim, r, thresh=.4):
    '''Generate a circle of specified radius (r) in pixel
    units.  Determines sub-pixel weighting using adaptive mesh refinement.
    Mesh refinement terminates at pixels whose side length is <= the specified
    threshold (thresh).'''
    x,y = n.indices((dim,dim), dtype=n.float)
    x -= dim/2 ; y -= dim/2
    rin,rout = int(r/sqrt2)-1, int(r)+1
    d1,d2,d3,d4 = dim/2-rout,dim/2-rin,dim/2+rin,dim/2+rout
    # If big circle, start as 1 and set a bounding box to 0.  
    # If small, start as 0 and set a bounded box to 1.
    if r > dim/2:
        rv = n.ones((dim,dim), dtype=n.float)
        rv[d1:d4,d1:d4] = 0
    else:
        rv = n.zeros((dim,dim), dtype=n.float)
        rv[d2:d3,d2:d3] = 1
    # Select 4 rects that contain boundary areas and evaluate them in detail
    for a1,a2,a3,a4 in ((d1,d2,d1,d4), (d3,d4,d1,d4),
            (d2,d3,d1,d2), (d2,d3,d3,d4)):
        x_, y_ = x[a1:a2,a3:a4], y[a1:a2,a3:a4]
        rs = n.sqrt(x_**2 + y_**2)
        # Get rough answer
        rv_ = (rs <= r).astype(n.float)
        # Fine-tune the answer
        brd = n.argwhere(n.abs(rs.flatten() - r) < 1 / sqrt2).flatten()
        rv_.flat[brd] = _circ(x_.flat[brd], y_.flat[brd], r, 1., thresh)
        # Set rectangle in the actual matrix
        rv[a1:a2,a3:a4] = rv_
    return rv
def _circ(x, y, r, p, thresh):
    # Subdivide into 9 pixels
    p /= 3.
    x0,x1,x2 = x, x+p, x-p
    y0,y1,y2 = y, y+p, y-p
    x = n.array([x0,x0,x0,x1,x1,x1,x2,x2,x2]).flatten()
    y = n.array([y0,y1,y2,y0,y1,y2,y0,y1,y2]).flatten()
    r2 = x**2 + y**2
    # Get the rough answer
    rv = (r2 <= r**2).astype(n.float) * p**2
    # Fine-tune the answer
    if p > thresh:
        brd = n.argwhere(n.abs(n.sqrt(r2) - r) < p / sqrt2).flatten()
        rv[brd] = _circ(x[brd], y[brd], r, p, thresh)
    rv.shape = (9, rv.size / 9)
    rv = rv.sum(axis=0)
    return rv


# Pair dim/dbm files
pairs = {}
for arg in args:
    name = arg[:-len('dim.fits')]
    pairs[name] = pairs.get(name, []) + [arg]
cmap = p.get_cmap(opts.cmap)
if opts.batch: m1,m2 = 1,1
else:
    m2 = int(math.sqrt(len(pairs)))
    m1 = int(math.ceil(float(len(pairs)) / m2))

#chan_regx = re.compile(r'_c(\d+)_(\d+)\.')
chan_regx = re.compile(r'_c(\d+)\.')
sfreq = 0.121142578125
sdf = 7.32421875e-05
pb_poly = [1.02854332e+09, -9.49707493e+08, 3.64775002e+08, -7.46038156e+07, 
           8.56951433e+06, -5.24246222e+05, 1.33464786e+04]

for cnt, filename in enumerate(pairs):
    print filename
    #chans = map(int, chan_regx.search(filename).groups())
    #freq = sfreq + sdf * .5 * (chans[0] + chans[1])
    #freq = sfreq + sdf * chans[0]
    print 'Assuming freq = 150 MHz'
    freq = .150
    wavelen = a.const.c / (freq * 1e9)
    pbeam = n.polyval(pb_poly, freq)
    dim,dbm = pairs[filename]
    if dim.find('dbm.fits') != -1: dim,dbm = dbm,dim
    if opts.batch:
        cnt = 0
        outfile = filename+'.png'
        if os.path.exists(outfile):
            print 'Output file exists... skipping.'
            continue
    # Gather data
    dim, kwds = a.img.from_fits(dim)
    dbm, kwds = a.img.from_fits(dbm)
    print dbm.shape
    print kwds
    print '-----------------------------------------------------------'

    # Some parameters determining size of pixels/image
    dpx_dg = kwds['d_ra']
    dpx_rd = dpx_dg * a.img.deg2rad
    sh = dim.shape[:2]
    DIM = sh[0]

    # Parse command-line options
    compress_axes = []
    ra_ax,dec_ax = (0,1)
    try:
        for i,ax in enumerate(kwds['axes']):
            if ax.startswith('ra'): ra_ax = i
            elif ax.startswith('dec'): dec_ax = i
            elif ax.startswith('freq'):
                chans = a.scripting.parse_chans(opts.chan, dim.shape[i])
                dim = dim.take(chans, axis=i)
                dbm = dbm.take(chans, axis=i)
                compress_axes.append(i)
            elif ax.startswith('stokes'):
                dim = dim.take([opts.pol], axis=i)
                dbm = dbm.take([opts.pol], axis=i)
                compress_axes.append(i)
    except(KeyError): pass

    compress_axes.reverse()
    for ax in compress_axes:
        dim = n.average(dim, axis=ax)
        dbm = n.average(dbm, axis=ax)

    # Put array in (ra,dec) order for plotting
    dim = dim.transpose((ra_ax,dec_ax))
    dbm = dbm.transpose((ra_ax,dec_ax))

    if not opts.fov is None:
        rpx = opts.fov / dpx_dg
        print 'FoV: rpx=%f' % (rpx)
        #fov = n.fromfunction(lambda x,y: gen_circ(x,y,r=rpx), sh)
        fov = circ(DIM,r=rpx)
        #fov = a.img.gaussian_beam(rpx, sh, center=(DIM/2,DIM/2))
        dim *= fov; dbm *= fov


    # Generate plots
    dim = n.fft.fft2(dim)
    dbm = n.fft.fft2(dbm)
    d = dim / n.where(dbm == 0, 1, dbm)
    # Generate T scaling
    jy_to_T_scaling = 1e-23 * wavelen**2 / (2 * a.const.k * pbeam)
    #jy_to_T_scaling = 1
    print 'Temp scaling using \\nu=%f \\lambda=%f, \\Omega=%f: %f' % (freq, wavelen, pbeam, jy_to_T_scaling)
    d *= jy_to_T_scaling
    d = a.img.recenter(d, (d.shape[0]/2, d.shape[1]/2))
    if opts.mode.startswith('phs'): d = n.angle(d)
    elif opts.mode.startswith('lin'): d = n.ma.absolute(d)
    elif opts.mode.startswith('real'): d = d.real
    elif opts.mode.startswith('imag'): d = d.imag
    elif opts.mode.startswith('log'):
        d = n.ma.absolute(d)
        d = n.ma.masked_less_equal(d, 0)
        d = n.ma.log10(d)

    if not opts.max is None: max = opts.max
    else: max = d.max()
    if not opts.drng is None: min = max - opts.drng
    else: min = d.min()

    p.subplot(m2, m1, cnt+1)
    p.imshow(d, vmin=min, vmax=max, origin='lower', cmap=cmap, interpolation='nearest')
    p.colorbar(shrink=.5, fraction=.05)
    p.title(filename)

    if opts.batch:
        print 'Saving to', outfile
        p.savefig(outfile)
        p.clf()
        

# Add right-click functionality for finding locations/strengths in map.
cnt = 1
def click(event):
    global cnt
    if not event.button in [2,3]: return
    xpx = n.around(event.xdata)
    ypx = n.around(event.ydata)
    flx = d[ypx,xpx]
    if opts.mode.startswith('log'): flx = 10**flx
    print '#%d PX: (%d,%d) Jy: %f' % (cnt, xpx, ypx, flx)
    cnt += 1

#register this function with the event handler
p.connect('button_press_event', click)

if not opts.batch:
    if opts.outfile != '':
        print 'Saving to', opts.outfile
        p.savefig(opts.outfile)
    else: p.show()
