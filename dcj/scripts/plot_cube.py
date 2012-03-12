#! /usr/bin/env python
"""
This is a general-purpose script for plotting simple FITS images.

Author: Aaron Parsons
Cribbed and mashed by: Danny Jacobs
"""

import aipy as a, sys, optparse, os
import numpy as n, pylab as p, ephem, math
from mpl_toolkits.basemap import Basemap

o = optparse.OptionParser()
o.set_usage('plot_cube.py [options] *.fits')
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
o.add_option('--nogrid', dest='nogrid', action='store_true',
    help='Do not display RA/DEC grid.')
o.add_option('-f', '--fft', dest='fft', action='store_true',
    help='Perform 2D FFT of image.')
o.add_option('--moments',dest='moments',
    help="Select from any of the first three moments: m0,m1,m2.")
opts, args = o.parse_args(sys.argv[1:])

moments = opts.moments.split(',')
cmap = p.get_cmap(opts.cmap)
if opts.batch: m1,m2 = 1,1
elif not opts.moments is None: m1,m2=len(moments),1
else:
    m2 = int(math.sqrt(len(args)))
    m1 = int(math.ceil(float(len(args)) / m2))

for cnt, filename in enumerate(args):
    print filename
    if opts.batch:
        cnt = 0
        outfile = filename+'.png'
        if os.path.exists(outfile):
            print 'Output file exists... skipping.'
            continue
    # Gather data
    d, kwds = a.img.from_fits(filename)
    print d.shape
    print kwds
    print '-----------------------------------------------------------'
    
    # Parse command-line options
    compress_axes = []
    ra_ax,dec_ax = (0,1)
    try:
        for i,ax in enumerate(kwds['axes']):
            if ax.startswith('ra'): ra_ax = i
            elif ax.startswith('dec'): dec_ax = i
            elif ax.startswith('freq'):
                chans = a.scripting.parse_chans(opts.chan, d.shape[i])
                d = d.take(chans, axis=i)
                freq_ax = i
               # compress_axes.append(i)
            elif ax.startswith('stokes'):
                d = d.take([opts.pol], axis=i)
                compress_axes.append(i)
                pol_ax = i
    except(KeyError): pass
    cube = d
    for mcnt,mom in enumerate(moments):
        compress_axes.reverse()
        for ax in compress_axes: d = n.average(cube, axis=ax)
        print d.shape
        if mom=='m0': d = n.average(d,axis=freq_ax)
        elif mom=='m1': d=n.average(n.diff(n.log10(d),axis=freq_ax),axis=freq_ax)
        elif mom=='m2': d=n.average(n.diff(n.diff(n.log10(d),axis=freq_ax),axis=freq_ax),axis=freq_ax)
        # Put array in (ra,dec) order for plotting
        print d.shape
        d = d.transpose((ra_ax,dec_ax))

        # Generate plots
        if opts.fft:
            d = n.fft.fft2(d)
            d = a.img.recenter(d, (d.shape[0]/2, d.shape[1]/2))
        if not mom in ['m1','m2']:
            if opts.mode.startswith('phs'): d = n.angle(d.filled(0))
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
        
        p.subplot(m2, m1, mcnt+1)
        p.title(mom)
        if not opts.nogrid:
            xpx,ypx = d.shape
            dx1 = -(xpx/2 + .5) * kwds['d_ra'] * a.img.deg2rad
            dx2 = (xpx/2 - .5) * kwds['d_ra'] * a.img.deg2rad
            dy1 = -(ypx/2 + .5) * kwds['d_dec'] * a.img.deg2rad
            dy2 = (ypx/2 - .5) * kwds['d_dec'] * a.img.deg2rad
            map = Basemap(projection='ortho', lon_0=180, lat_0=kwds['dec'],
                rsphere=1, llcrnrx=dx1, llcrnry=dy1, urcrnrx=dx2,urcrnry=dy2)
            map.drawmeridians(n.arange(kwds['ra']-180,kwds['ra']+180,30))
            map.drawparallels(n.arange(-90,120,30))
            map.drawmapboundary()
            map.imshow(d, vmin=min, vmax=max, cmap=cmap)
        else: p.imshow(d, vmin=min, vmax=max, origin='lower', cmap=cmap)
    p.colorbar(shrink=.5, fraction=.05)
    p.suptitle(filename)

    if opts.batch:
        print 'Saving to', outfile
        p.savefig(outfile)
        p.clf()
        

# Add right-click functionality for finding locations/strengths in map.
cnt = 1
def click(event):
    global cnt
    if event.button == 3: 
        lon,lat = map(event.xdata, event.ydata, inverse=True)
        lon = (180 + kwds['ra'] - lon) % 360
        lon *= a.img.deg2rad; lat *= a.img.deg2rad
        ra,dec = ephem.hours(lon), ephem.degrees(lat)
        xpx = n.around(event.xdata / (kwds['d_ra'] * a.img.deg2rad))
        ypx = n.around(event.ydata / (kwds['d_dec'] * a.img.deg2rad))
        flx = d[ypx,xpx]
        if opts.mode.startswith('log'): flx = 10**flx
        print '#%d (RA,DEC): (%s, %s), PX: (%d,%d) Jy: %f' % \
            (cnt, ra, dec, xpx, ypx, flx)
        cnt += 1
    elif event.button == 2: 
        p.figure(200)
        lon,lat = map(event.xdata, event.ydata, inverse=True)
        p.plot(cube[lat,lon,:])
        lon = (180 + kwds['ra'] - lon) % 360
        lon *= a.img.deg2rad; lat *= a.img.deg2rad
        ra,dec = ephem.hours(lon), ephem.degrees(lat)
        xpx = n.around(event.xdata / (kwds['d_ra'] * a.img.deg2rad))
        ypx = n.around(event.ydata / (kwds['d_dec'] * a.img.deg2rad))
        flx = d[ypx,xpx]
        if opts.mode.startswith('log'): flx = 10**flx
        p.title('#%d (RA,DEC): (%s, %s), PX: (%d,%d) Jy: %f' % \
            (cnt, ra, dec, xpx, ypx, flx))
        cnt += 1
    else: return

    


#register this function with the event handler
p.connect('button_press_event', click)

if not opts.batch:
    if opts.outfile != '':
        print 'Saving to', opts.outfile
        p.savefig(opts.outfile)
    else: p.show()
