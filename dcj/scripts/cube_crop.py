#!/usr/bin/env python
#
#  cube_crop.py
#  
#
#  Created by Danny Jacobs on 1/10/10.
#  PAPER Project
#

import aipy as a, numpy as n, pylab as p,math as m
import sys, optparse,ephem,datetime

o = optparse.OptionParser()
o.set_usage('cube_crop.py --crop <degrees> *.fits')
#a.scripting.add_standard_options(o, cal=True)
o.add_option('--crop',dest='crop',default=70,type=float,
    help='Width of square central crop region, 0 to skip, default=70.')
opts, args = o.parse_args(sys.argv[1:])

def hname(ra):
	one = str((int(str(ephem.hours(ra*n.pi/180)).split(':')[0])+100000))[-2::]
	two =  str((int(str(ephem.hours(ra*n.pi/180)).split(':')[1])+100000))[-2::]
	ret = one + two
	return ret
def dname(dec):
	one = str((int(str(ephem.degrees(dec*n.pi/180)).split(':')[0])+100000))[-2::]
	two = str((int(str(ephem.degrees(dec*n.pi/180)).split(':')[1])+100000))[-2::]
	if dec < 0: add = '-'
	else: add = '+'
	ret = add + one + two
	return ret

for file in args:
    s = file.split('.')
    out_file = '.'.join(['.'.join(s[:-1]),'c',s[-1]])
    print file + ' > ' + out_file
    c,kwrd = a.img.from_fits(file)
    pic_res = (n.abs(kwrd['d_ra']) + n.abs(kwrd['d_dec']))/2 
    cen = n.round(n.array(c.shape).astype(float)[:2]/2)
    print "image center at pixels",kwrd
    hw = round(opts.crop/pic_res/2)
    print "extracting central %d pixels"%hw
    cnrs = n.vstack([cen,cen])                              #corner matrix 
    cnrs += n.vstack([[-hw,-hw],[hw,hw]])                   #crnrs[0,:]=blc=[y,x]
    print "with corners ",cnrs                              #crnrs[1,:]=trc=[y,x]
    print "input shape:",c.shape
    cc = c[cnrs[0,0]:cnrs[1,0],cnrs[0,1]:cnrs[1,1]]     #crop
    print "output shape:",cc.shape
    history = "Cropped to central %3.1f degrees"%(opts.crop,)
    kwrd['object'] = 'J'+hname(kwrd['ra'])+dname(kwrd['dec'])
    a.img.from_fits_to_fits(file,out_file,cc,kwrd)

#    try: 
#        print kwrd['history']
#        history = kwrd['history'] + history
#    except(KeyError): pass
#
#    
#    a.img.to_fits(out_file, 
#        cc, 
#        clobber=True, 
#        axes=('ra---sin', 'dec--sin','freq    ','stokes  '), 
#        object='J'+hname(kwrd['ra'])+dname(kwrd['dec']), 
#        telescope='PAPER', 
#        instrument='', 
#        observer='', 
#        origin='AIPY', 
#        obs_date=kwrd['obs_date'], 
#        cur_date=str(datetime.date.today()), 
#        ra=kwrd['ra'], dec=kwrd['dec'], 
#        d_ra=kwrd['d_ra'], d_dec=kwrd['d_dec'], epoch=2000.0,
#        freq=kwrd['freq'],d_freq=kwrd['d_freq'],
#        history=history)
