#!/usr/bin/env python
#
#  mkcubes.py
#  
#
#  Created by Charles Epstein on 8/4/09.
#  Copyright (c) 2009 __MyCompanyName__. All rights reserved.
#
import matplotlib ; matplotlib.use('Agg')
import glob, aipy as a, pyfits, operator as op, numpy as n, os, ephem
from pylab import *

inputdir = '/data1/paper/arp/pgb015/lst_v004/map009/'
outputdir = '/data1/paper/summer09/newcubes/'

files = glob.glob(inputdir + '*.bim.fits')

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


cubes={}
flist = []
for file in files:
	flist.append(file.split('/')[-1])
for file in flist:
	name = file.split('_')
	num = name[4].split('.')
	cubes[num[0]]=[]
for file in flist:
	name = file.split('_')
	num = name[4].split('.')
	cubes[num[0]].append(file)
mlist = cubes.items()
flist=[]
for m in mlist:
	flist.append(m[0])
for i in range(0,len(flist)):
	flist[i]=[flist[i]]
	flist[i].append(mlist[i][1])
for x in range(0,len(flist)):
	value = sorted(flist[x][1])
	flist[x][1] = value[15:23]+value[0:15]

for cube in flist:
	cubedata = []
	for i in range(0,len(cube[1])):
		data, kwds = a.img.from_fits(inputdir+cube[1][i])
		cubedata.append(data)
	cdata = array(cubedata[0], dtype = 'float')
	for i in range(0,len(cubedata)-1):
		cdata = n.dstack([cdata,array(cubedata[i+1], dtype = 'float')])

	filename = outputdir +'J'+hname(kwds['ra'])+dname(kwds['dec']) +'.fits'
		
	a.img.to_fits(filename, cdata, clobber=False, axes=('DEC--SIN', 
		'RA---SIN','freq'), object='', telescope='', instrument='', 
		observer='', origin='AIPY', obs_date=kwds['obs_date'], cur_date='', 
		ra=kwds['ra'], dec=kwds['dec'], d_ra=-kwds['d_ra'], d_dec=kwds['d_dec'], epoch=2000.0, freq=60e6, d_freq=5e6, bscale=0, 
		bzero=1)
	file = pyfits.open(filename, mode='update')
	file[0].header['CRPIX1'] = 500
	file[0].header['CRPIX2'] = 500
	file[0].header['CRPIX3'] = 1
	file.flush()
	print 'Done a cube'